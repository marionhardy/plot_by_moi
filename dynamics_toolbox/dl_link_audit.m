function L = dl_link_audit(varargin)
% dl_link_audit  Is each live track paired with the right CyCIF cell? Per-well audit of
%   track integrity, live -> endpoint linking and mask alignment.
%   A broken pairing (track identity swaps, mis-registered or mis-aligned masks, cells
%   moving before fixation) erases any dynamics -> endpoint relation while the live data
%   and the stain each still look fine on their own. This audit looks for that, per well,
%   with model-free metrics, and (optionally) relates them to how well the dynamics
%   model predicts each well.
%
% USAGE
%   P = load(pipeFile, 'imputeMask','expPaths','hasEndpoint','winAlign');
%   L = dl_link_audit('S',S, 'data',data, 'imputeMask',P.imputeMask, 'rep',repCol, 'xy',origXYcol, ...
%                     'cell_id',origCellID, 'rows',isBioAll, 'has_endpoint',P.hasEndpoint, ...
%                     'endpoint_ok',endpointOK, 'channels',channels, 'dt_min',dtMin, ...
%                     'win_align',P.winAlign, 'exp_paths',P.expPaths, 'pred',R3, 'table',A);
%   L = dl_link_audit('S',S, 'data',data, 'rep',repCol, 'xy',origXYcol, 'channels',channels, ...
%                     'dt_min',dtMin, 'plots','none');                          % model-free, no dataloc
%   L = dl_link_audit(..., 'highlight',["3|12","3|14"]);   % wells you re-aligned (rep|xy), drawn in red
%   L = dl_link_audit(..., 'contrast',{'CD68','auto'; 'DAPI','nuc'});   % chosen compartment checks
%
% NAME-VALUE (paired)
%  Required
%   'S'            master object with .cycifT (row-aligned to data; <marker>_nuc / _cyt columns)
%   'data'         [n x T x C] live traces in the encoder window (as saved by the pipeline)
%   'rep', 'xy'    [n x 1] replicate (index into exp_paths) and original XY number
%  Optional, per cell (row-aligned to data)
%   'imputeMask'   [n x T x C] true = imputed point                                   ([])
%   'cell_id'      original cell id (dataloc.d{xy}.cellindex), needed for full traces  ([])
%   'rows'         logical subset to audit (e.g. isBioAll)                            ([] = all)
%   'has_endpoint' any CyCIF value                                       ([] = from S.cycifT)
%   'endpoint_ok'  endpoint not below secondary-only background (pipeline 3a)  ([] = all true)
%  Track steps (identity swaps change the level of EVERY sensor at the same frame)
%   'channels' ({})  'dt_min' (10)  'win_align' ('first' | 'last')                    ('first')
%   'step_k'       level change between the median of the 'step_win' frames before
%                  and after a frame, in units of the cell's frame noise              (8)
%   'step_win'     frames per side for the step medians                               (5)
%   'late_min'     "late" = this many minutes before the end of the movie (fixation)  (90)
%   'sync_frac'    frames where more than this fraction of a well's cells co-step are
%                  well-wide events (stimulus addition, focus/stage jumps), not identity
%                  swaps (and >= 3 cells): removed from the swap counts, reported as
%                  n_sync_frames                                                       (0.05)
%  Endpoint / mask alignment
%   'contrast'     'all' = every marker with _nuc and _cyt columns, or {marker, 'cyt'|'nuc'|'auto'; ...}.
%                  Per cell (cyt - nuc)/(|cyt| + |nuc|), signed so that POSITIVE = the expected
%                  compartment ('auto': the sign of the pooled median, i.e. whichever compartment
%                  is brighter overall). Masks off their cells sample a mix -> every marker's
%                  contrast shrinks toward 0. contrast_mean = mean over markers per well.
%                  Only ONE marker shifting = stain/biology; ALL shifting = mask alignment   ('all')
%  Registration and full traces (loads each dataloc once)
%   'exp_paths'    cellstr, one experiment path per replicate (dl_dataloc_v3 'exp')   ({})
%  Model link (optional)
%   'pred'         result of dl_dynamics_predict (uses .pred and .keep_idx)          ([])
%   'table'        the table that was passed to dl_dynamics_predict                  ([])
%   'rep_col' ('rep')  'unit_col' ('xy')  'min_cells' (10, per well for correlations)
%  Output
%   'highlight' (string "rep|xy")  'flag_z' (3)  'plots' 'summary' | 'none' ('summary')
%   'verbose' (true)  'save' (true)  'savepath' ('')  'closefigs' (true)
%
% OUTPUT  L struct
%   .wells   one row per well: well | rep | xy | n | pct_endpoint | pct_ok | pct_imputed |
%            pct_costep (>= 1 simultaneous step in all sensors, window) | pct_late_costep
%            (last late_min of the window) | n_sync_frames (well-wide jumps: stimulus,
%            focus/stage) | contrast_<marker> | [reg_post_px | reg_improve_px |
%            reg_failed | pct_fix_costep (co-step in the last late_min before fixation, full
%            trace) | pct_absent_last (endpoint present but track absent in the last frame) |
%            pct_id_missing] | [r_pred (mean within-well corr(observed, predicted))] | flags
%   .reps    per-replicate medians over wells + frames_full, gap_h (window end -> fixation)
%   .rep_tests  per replicate and metric: median vs the other replicates, rank-sum p
%   .r_pred  [wells x targets] within-well corr(observed, out-of-fold prediction)
%   .cells   per audited cell: well | imputed | n_costep | late_costep | fix_costep | absent_last
%   .figs
%
% NOTES
%   New in DatalocHandler v5 (no v4 counterpart).
%   Steps are counted on running medians, so single-frame spikes are ignored; a step in
%   ONE sensor can be biology (a sustained pulse onset), a step in ALL sensors at the same
%   frame is the signature of a track jumping to another cell.
%   Contrast assumes the live endpoint masks were warped onto the CyCIF images
%   (ct_cycif_link): a mask that lands off its cell samples neighbours/background, and the
%   expected compartment difference shrinks toward 0.
%   Registration QC is read from dataloc.d{xy}.cycif_reg (worst round per well).
%   Requires Statistics and Machine Learning Toolbox (mad, ranksum).

%% ---- Parse ----------------------------------------------------------------
p = inputParser;  p.FunctionName = mfilename;
isTxt = @(x) ischar(x) || isstring(x);
addParameter(p, 'S',            [],  @isstruct);
addParameter(p, 'data',         [],  @isnumeric);
addParameter(p, 'imputeMask',   []);
addParameter(p, 'rep',          []);
addParameter(p, 'xy',           []);
addParameter(p, 'cell_id',      []);
addParameter(p, 'rows',         []);
addParameter(p, 'has_endpoint', []);
addParameter(p, 'endpoint_ok',  []);
addParameter(p, 'channels',     {},  @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(p, 'dt_min',       10,  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'win_align',    'first', @(x) any(strcmpi(x, {'first','last'})));
addParameter(p, 'step_k',       8);
addParameter(p, 'step_win',     5);
addParameter(p, 'late_min',     90);
addParameter(p, 'sync_frac',    0.05);
addParameter(p, 'contrast',     'all', @(x) iscell(x) || (isTxt(x) && strcmpi(x, 'all')));
addParameter(p, 'exp_paths',    {});
addParameter(p, 'pred',         []);
addParameter(p, 'table',        []);
addParameter(p, 'rep_col',      'rep', isTxt);
addParameter(p, 'unit_col',     'xy',  isTxt);
addParameter(p, 'min_cells',    10);
addParameter(p, 'highlight',    strings(0, 1));
addParameter(p, 'flag_z',       3);
addParameter(p, 'plots',        'summary', @(x) any(strcmpi(x, {'summary','none'})));
addParameter(p, 'verbose',      true);
addParameter(p, 'save',         true);
addParameter(p, 'savepath',     '',  isTxt);
addParameter(p, 'closefigs',    true);
parse(p, varargin{:});  o = p.Results;
tag = '[dl_link_audit]';

X = double(o.data);
assert(~isempty(X) && ~isempty(o.S) && isfield(o.S, 'cycifT'), '%s ''S'' (with cycifT) and ''data'' are required.', tag);
[n, T, C] = size(X);
ch = cellstr(o.channels);  if isempty(ch), ch = compose('ch%d', 1:C); end
assert(numel(ch) == C, '%s %d channel names for %d channels.', tag, numel(ch), C);
assert(numel(o.rep) == n && numel(o.xy) == n, '%s ''rep'' and ''xy'' need one entry per row of data.', tag);
rep = double(o.rep(:));  xy = double(o.xy(:));
key = string(rep) + "|" + string(xy);
rows = true(n, 1);
if ~isempty(o.rows)
    if islogical(o.rows), rows = o.rows(:); else, rows = false(n, 1);  rows(o.rows) = true; end
end
CT = o.S.cycifT;
assert(height(CT) == n, '%s S.cycifT has %d rows, data %d.', tag, height(CT), n);
isNum = varfun(@isnumeric, CT, 'OutputFormat', 'uniform');
hasEP = o.has_endpoint;  if isempty(hasEP), hasEP = any(~isnan(double(CT{:, isNum})), 2); end
hasEP = logical(hasEP(:));
epOK  = o.endpoint_ok;   if isempty(epOK), epOK = true(n, 1); end
epOK  = logical(epOK(:));
w     = o.step_win;  lateF = max(1, round(o.late_min / o.dt_min));
cid   = o.cell_id;   if ~isempty(cid), cid = double(cid(:)); end

%% ---- Per cell: imputation and steps inside the window ---------------------------
imp = nan(n, 1);
if ~isempty(o.imputeMask)
    assert(size(o.imputeMask, 1) == n, '%s imputeMask has %d rows, data %d.', tag, size(o.imputeMask, 1), n);
    imp = 100 * mean(reshape(double(o.imputeMask), n, []), 2);
end
ia = find(rows);
[wells, ~, gw] = unique(key(ia));  nW = numel(wells);
CO = i_steps(X, o.step_k, w);
[CO, nSync] = i_desync(CO, ia, gw, nW, o.sync_frac);             % drop well-wide (synchronous) frames
nCo    = i_count(CO);
lateCo = any(CO(:, max(1, end - lateF + 1):end), 2);
if strcmpi(o.win_align, 'first') && o.verbose
    fprintf('%s window = first %d frames: window steps do not cover the frames before fixation (use ''exp_paths'').\n', tag, T);
end

%% ---- Per cell: compartment contrast (mask on its cell?) -------------------------
vn = CT.Properties.VariableNames;
if ~iscell(o.contrast)                                            % 'all'
    base = erase(vn(endsWith(vn, '_nuc')), '_nuc');
    base = base(ismember(strcat(base, '_cyt'), vn));
    cs = [base(:), repmat({'auto'}, numel(base), 1)];
else
    cs = o.contrast;  if ~isempty(cs) && size(cs, 2) ~= 2, cs = reshape(cs, 2, [])'; end
end
CTR = nan(n, 0);  ctrNames = strings(1, 0);
for j = 1:size(cs, 1)
    m = char(cs{j, 1});
    if ~all(ismember({[m '_nuc'], [m '_cyt']}, vn))
        warning('%s contrast %s: needs %s_nuc and %s_cyt in S.cycifT; skipped.', tag, m, m, m);  continue
    end
    nu = double(CT.([m '_nuc']));  cy = double(CT.([m '_cyt']));
    v  = (cy - nu) ./ (abs(cy) + abs(nu));
    sg = 1;
    if strcmpi(cs{j, 2}, 'nuc'), sg = -1;
    elseif strcmpi(cs{j, 2}, 'auto'), sg = sign(median(v(rows & hasEP), 'omitnan'));  if sg == 0 || isnan(sg), sg = 1; end
    end
    v = sg * v;
    CTR(:, end+1) = v;  ctrNames(end+1) = "contrast_" + string(m); %#ok<AGROW>
end

%% ---- Per well -----------------------------------------------------------------
first = accumarray(gw, ia, [], @min);
Wt = table(wells, rep(first), xy(first), accumarray(gw, 1), 'VariableNames', {'well','rep','xy','n'});
mw = @(v) accumarray(gw, v, [nW 1], @(z) mean(z, 'omitnan'));
Wt.pct_endpoint    = 100 * mw(double(hasEP(ia)));
okv = double(epOK(ia));  okv(~hasEP(ia)) = NaN;
Wt.pct_ok          = 100 * mw(okv);
Wt.pct_imputed     = mw(imp(ia));
Wt.pct_costep      = 100 * mw(double(nCo(ia) > 0));
Wt.pct_late_costep = 100 * mw(double(lateCo(ia)));
Wt.n_sync_frames   = nSync;
for j = 1:numel(ctrNames)
    v = CTR(ia, j);  v(~hasEP(ia)) = NaN;
    Wt.(char(ctrNames(j))) = accumarray(gw, v, [nW 1], @(z) median(z, 'omitnan'));
end
if numel(ctrNames) > 1
    Wt.contrast_mean = mean(Wt{:, cellstr(ctrNames)}, 2, 'omitnan');
    ctrNames = [ctrNames, "contrast_mean"];
end

%% ---- Registration QC + full traces from each dataloc -------------------------------
fixCo = nan(n, 1);  absLast = nan(n, 1);  idMiss = nan(n, 1);
repInfo = table(unique(rep(ia)), 'VariableNames', {'rep'});
repInfo.frames_full = nan(height(repInfo), 1);
if ~isempty(o.exp_paths)
    [Wt.reg_post_px, Wt.reg_improve_px, Wt.reg_failed] = deal(nan(nW, 1));
    for r = repInfo.rep'
        if r > numel(o.exp_paths), warning('%s rep %d has no entry in exp_paths.', tag, r);  continue; end
        if o.verbose, fprintf('%s loading rep %d dataloc ...\n', tag, r); end
        dl = dl_dataloc_v3('exp', o.exp_paths{r});
        for k = find(Wt.rep == r)'
            x0 = Wt.xy(k);
            if x0 > numel(dl.d) || isempty(dl.d{x0}), continue; end
            dd = dl.d{x0};
            if isfield(dd, 'cycif_reg') && ~isempty(dd.cycif_reg)
                q = dd.cycif_reg;
                Wt.reg_post_px(k)    = max(i_num(q, 'median_post_px'));       % worst round
                Wt.reg_improve_px(k) = min(i_num(q, 'improvement_px'));
                Wt.reg_failed(k)     = nnz(i_num(q, 'failed') == 1);
            end
            if isempty(cid) || ~isfield(dd, 'cellindex') || ~all(isfield(dd.data, ch)), continue; end
            ci = ia(gw == k);
            [tf, loc] = ismember(cid(ci), double(dd.cellindex(:)));
            Tf = size(dd.data.(ch{1}), 2);
            repInfo.frames_full(repInfo.rep == r) = Tf;
            Xf = nan(numel(ci), Tf, C);
            for c = 1:C, Xf(tf, :, c) = double(dd.data.(ch{c})(loc(tf), :)); end
            COf = i_steps(Xf, o.step_k, w);
            COf(:, mean(COf, 1) > o.sync_frac & sum(COf, 1) >= 3) = false;  % well-wide frames
            fixCo(ci)   = double(any(COf(:, max(1, end - lateF + 1):end), 2));
            absLast(ci) = double(any(isnan(Xf(:, end, :)), 3));
            idMiss(ci)  = double(~tf);
            fixCo(ci(~tf)) = NaN;  absLast(ci(~tf)) = NaN;
        end
        clear dl
    end
    Wt.pct_fix_costep  = 100 * mw(fixCo(ia));
    a = absLast(ia);  a(~hasEP(ia)) = NaN;                         % only cells that HAVE an endpoint
    Wt.pct_absent_last = 100 * mw(a);
    Wt.pct_id_missing  = 100 * mw(idMiss(ia));
end
repInfo.gap_h = zeros(height(repInfo), 1);
if strcmpi(o.win_align, 'first'), repInfo.gap_h = (repInfo.frames_full - T) * o.dt_min / 60; end

%% ---- Model link: within-well corr(observed, out-of-fold prediction) ----------------
RP = nan(nW, 0);  rpNames = strings(1, 0);
if ~isempty(o.pred) && ~isempty(o.table)
    R  = o.pred;  Tk = o.table(R.keep_idx, :);
    kk = string(Tk.(o.rep_col)) + "|" + string(Tk.(o.unit_col));
    [~, wl] = ismember(kk, Wt.well);
    for t = 1:numel(R.pred)
        col = i_col(Tk, char(R.pred(t).target));
        if isempty(col), continue; end
        rt = ~isnan(double(Tk.(col)));  y = R.pred(t).y(:);  yh = R.pred(t).yhat(:);
        if numel(y) ~= nnz(rt)
            warning('%s %s: %d predictions for %d measured rows; skipped.', tag, R.pred(t).target, numel(y), nnz(rt));  continue
        end
        g = wl(rt);  v = nan(nW, 1);
        for k = 1:nW
            i = g == k & ~isnan(y) & ~isnan(yh);
            if nnz(i) >= o.min_cells && std(yh(i)) > 0, v(k) = corr(y(i), yh(i)); end
        end
        RP(:, end+1) = v;  rpNames(end+1) = string(R.pred(t).target); %#ok<AGROW>
    end
    if ~isempty(RP), Wt.r_pred = mean(RP, 2, 'omitnan'); end
end

%% ---- Flags (robust z over all wells, in the BAD direction) -------------------------
bad = {'pct_endpoint',-1; 'pct_ok',-1; 'pct_imputed',1; 'pct_costep',1; 'pct_late_costep',1; 'n_sync_frames',1; ...
       'reg_post_px',1; 'reg_improve_px',-1; 'pct_fix_costep',1; 'pct_absent_last',1; 'pct_id_missing',1; 'r_pred',-1};
bad = [bad; [cellstr(ctrNames)', num2cell(-ones(numel(ctrNames), 1))]];
bad = bad(ismember(bad(:, 1), Wt.Properties.VariableNames), :);
Wt.flags = strings(nW, 1);
for j = 1:size(bad, 1)
    v = Wt.(bad{j, 1});  s = 1.4826 * mad(v, 1);  if ~(s > 0), continue; end
    z = bad{j, 2} * (v - median(v, 'omitnan')) / s;
    f = z > o.flag_z;  Wt.flags(f) = strip(Wt.flags(f) + " " + string(bad{j, 1}));
end
if ismember('reg_failed', Wt.Properties.VariableNames)
    f = Wt.reg_failed > 0;  Wt.flags(f) = strip(Wt.flags(f) + " reg_failed");
end
Wt.highlight = ismember(Wt.well, string(o.highlight(:)));

%% ---- Per replicate -----------------------------------------------------------------------
mets = [bad(:, 1); {'reg_failed'}];  mets = mets(ismember(mets, Wt.Properties.VariableNames));
Rt = repInfo;
Rt.n_wells = arrayfun(@(r) nnz(Wt.rep == r), Rt.rep);
Rt.n_cells = arrayfun(@(r) sum(Wt.n(Wt.rep == r)), Rt.rep);
Rt.n_flagged = arrayfun(@(r) nnz(Wt.rep == r & strlength(Wt.flags) > 0), Rt.rep);
for j = 1:numel(mets)
    Rt.(mets{j}) = arrayfun(@(r) median(Wt.(mets{j})(Wt.rep == r), 'omitnan'), Rt.rep);
end
tests = {};
for r = Rt.rep'
    for j = 1:numel(mets)
        a = Wt.(mets{j})(Wt.rep == r);  b = Wt.(mets{j})(Wt.rep ~= r);
        a = a(~isnan(a));  b = b(~isnan(b));
        pv = NaN;  if numel(a) >= 3 && numel(b) >= 3, pv = ranksum(a, b); end
        dir = 0;  k = find(strcmp(bad(:, 1), mets{j}));  if ~isempty(k), dir = bad{k, 2}; end
        worse = dir * (median(a) - median(b)) > 0;
        tests(end+1, :) = {r, string(mets{j}), median(a), median(b), pv, worse}; %#ok<AGROW>
    end
end
RT = cell2table(tests, 'VariableNames', {'rep','metric','median_rep','median_others','p','worse'});

L = struct('wells', Wt, 'reps', Rt, 'rep_tests', RT, 'r_pred', array2table(RP, 'VariableNames', cellstr(rpNames)), ...
           'cells', table(key, imp, nCo, lateCo, fixCo, absLast, 'VariableNames', ...
                          {'well','pct_imputed','n_costep','late_costep','fix_costep','absent_last'}), ...
           'params', o, 'figs', gobjects(0));
L.r_pred.well = Wt.well;
% Compartment contrast per marker x replicate (median over wells)
L.contrast_by_rep = table();
cn = ctrNames(ctrNames ~= "contrast_mean");
if ~isempty(cn)
    M = nan(numel(cn), height(Rt));
    for j = 1:numel(cn), M(j, :) = arrayfun(@(r) median(Wt.(char(cn(j)))(Wt.rep == r), 'omitnan'), Rt.rep'); end
    L.contrast_by_rep = array2table(M, 'VariableNames', cellstr(compose('rep%d', Rt.rep')), 'RowNames', cellstr(erase(cn, "contrast_")));
end
% Within-replicate association r_pred ~ contrast (separates alignment from a replicate batch effect)
L.within_rep = table();
cm = '';  if ismember('contrast_mean', Wt.Properties.VariableNames), cm = 'contrast_mean'; elseif ~isempty(cn), cm = char(cn(1)); end
if ismember('r_pred', Wt.Properties.VariableNames) && ~isempty(cm)
    rr = [];
    for r = Rt.rep'
        i = Wt.rep == r & ~isnan(Wt.r_pred) & ~isnan(Wt.(cm));
        [rho, pv] = deal(NaN);  if nnz(i) >= 5, [rho, pv] = corr(Wt.(cm)(i), Wt.r_pred(i), 'Type', 'Spearman'); end
        rr = [rr; r, nnz(i), rho, pv]; %#ok<AGROW>
    end
    L.within_rep = array2table(rr, 'VariableNames', {'rep','n_wells','rho_rpred_contrast','p'});
end

%% ---- Console --------------------------------------------------------------------------
if o.verbose
    fprintf('\n%s %d cells in %d wells | co-step = level change > %g x frame noise in ALL %d sensors at the same frame\n', ...
            tag, numel(ia), nW, o.step_k, C);
    disp(Rt)
    for r = Rt.rep'
        s = RT(RT.rep == r & RT.worse & RT.p < 0.01, :);
        if isempty(s), fprintf('  rep %d: no metric worse than the other reps (rank-sum p < 0.01)\n', r);
        else
            fprintf('  rep %d WORSE: %s\n', r, strjoin(compose('%s %.3g vs %.3g (p %.2g)', s.metric, s.median_rep, ...
                    s.median_others, s.p), ' | '));
        end
    end
    fl = Wt(strlength(Wt.flags) > 0, :);
    if ~isempty(fl)
        [~, oi] = sort(arrayfun(@(s) numel(split(s)), fl.flags), 'descend');  fl = fl(oi, :);
        fprintf('  flagged wells (robust z > %g vs all wells): %d\n', o.flag_z, height(fl));
        for i = 1:min(15, height(fl)), fprintf('    %-8s n %4d  %s\n', fl.well(i), fl.n(i), fl.flags(i)); end
    end
    if ~isempty(L.contrast_by_rep)
        fprintf('  compartment contrast per marker (median over wells; signed so the pooled median > 0; ~0 = no compartment difference):\n');
        disp(L.contrast_by_rep)
    end
    if ~isempty(L.within_rep)
        fprintf('  within replicate: Spearman r_pred ~ %s across wells (positive = better-aligned wells are better predicted)\n', cm);
        disp(L.within_rep)
    end
    if any(Wt.highlight)
        h = Wt(Wt.highlight, :);
        fprintf('  highlighted wells: %s\n', strjoin(compose('%s [%s]', h.well, h.flags), ', '));
    end
end

%% ---- Figure -----------------------------------------------------------------------------
if strcmpi(o.plots, 'summary')
    L.figs = i_fig(Wt, mets, bad, o, tag);
    i_save(L.figs, o);
    if o.closefigs && o.save && ~isempty(o.savepath), close(L.figs);  L.figs = gobjects(0); end
end
end


%% =====================================================================================
function CO = i_steps(X, k, w)
% Level steps per frame: |median of the w frames after - median of the w frames before|
% in units of the cell's frame noise; CO = step in ALL channels at the same frame.
[n, T, C] = size(X);
S = false(n, T - 1, C);
for c = 1:C
    x  = X(:, :, c);
    sg = 1.4826 * mad(diff(x, 1, 2), 1, 2) / sqrt(2);  sg(~(sg > 0)) = NaN;
    mb = movmedian(x, [w-1 0], 2, 'omitnan', 'Endpoints', 'shrink');
    ma = movmedian(x, [0 w-1], 2, 'omitnan', 'Endpoints', 'shrink');
    S(:, :, c) = abs(ma(:, 2:end) - mb(:, 1:end-1)) ./ sg > k;
end
CO = all(S, 3);
end

function nCo = i_count(CO)
% Number of co-step EVENTS per cell (runs of flagged frames)
nCo = sum(diff(double([false(size(CO, 1), 1), CO]), 1, 2) == 1, 2);
end

function [CO, nSync] = i_desync(CO, ia, gw, nW, frac)
% Frames where > frac of a well's cells co-step are well-wide events: unflag them
nSync = zeros(nW, 1);
for k = 1:nW
    iw = ia(gw == k);  sy = mean(CO(iw, :), 1) > frac & sum(CO(iw, :), 1) >= 3;   % >= 3 cells: not 1 swap in a sparse well
    nSync(k) = nnz(sy);  CO(iw, sy) = false;
end
end

function v = i_num(q, fld)
% Numeric field of a struct array, NaN where missing / empty
v = nan(1, numel(q));
if ~isfield(q, fld), return; end
for i = 1:numel(q)
    x = q(i).(fld);
    if ~isempty(x) && (isnumeric(x) || islogical(x)), v(i) = double(x(1)); end
end
end

function col = i_col(T, nm)
% Column of T for target nm: exact name, else first '<nm>_<compartment>'
vn = T.Properties.VariableNames;  col = '';
if ismember(nm, vn), col = nm; return; end
k = find(startsWith(vn, [nm '_']), 1);  if ~isempty(k), col = vn{k}; end
end

function f = i_fig(Wt, mets, bad, o, tag) %#ok<INUSD>
% One strip per metric: wells as dots (size ~ cells), x = replicate; red = highlighted
mets = mets(cellfun(@(m) any(~isnan(Wt.(m))), mets));
if any(strcmp(mets, 'contrast_mean')), mets = mets(~startsWith(mets, 'contrast_') | strcmp(mets, 'contrast_mean')); end
hasR = ismember('r_pred', Wt.Properties.VariableNames);
ctr  = mets(startsWith(mets, 'contrast_'));
np   = numel(mets) + double(hasR && ~isempty(ctr));
nc   = min(4, np);  nr = ceil(np / nc);
f = figure('Color', 'w', 'Name', 'link audit', 'Units', 'pixels');
f.Position(3:4) = [330 * nc, 260 * nr + 60];
tl = tiledlayout(f, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {'Link audit: is each live track paired with the right CyCIF cell?', ...
           'dots = wells (size ~ cells); red = highlighted wells; arrow = direction of a problem'}, ...
      'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
rs = RandStream('mt19937ar', 'Seed', 0);
ur = unique(Wt.rep)';  cr = lines(numel(ur));
for j = 1:numel(mets)
    ax = nexttile(tl);  hold(ax, 'on');
    v = Wt.(mets{j});  xj = Wt.rep + 0.18 * (rand(rs, height(Wt), 1) - 0.5);
    sz = 8 + 40 * sqrt(Wt.n / max(Wt.n));
    for r = 1:numel(ur)
        i = Wt.rep == ur(r) & ~Wt.highlight;
        scatter(ax, xj(i), v(i), sz(i), cr(r, :), 'filled', 'MarkerFaceAlpha', 0.6);
        md = median(v(Wt.rep == ur(r)), 'omitnan');
        plot(ax, ur(r) + [-0.3 0.3], [md md], 'k-', 'LineWidth', 1.5);
    end
    i = Wt.highlight;
    if any(i), scatter(ax, xj(i), v(i), sz(i), [0.85 0.1 0.1], 'filled', 'MarkerEdgeColor', 'k'); end
    k = find(strcmp(bad(:, 1), mets{j}));  arr = '';
    if ~isempty(k), if bad{k, 2} > 0, arr = ' (problem: high)'; else, arr = ' (problem: low)'; end, end
    title(ax, [strrep(mets{j}, '_', ' ') arr], 'FontSize', 8.5, 'Interpreter', 'none');
    set(ax, 'XTick', ur, 'XTickLabel', compose('rep %d', ur), 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
    xlim(ax, [min(ur) - 0.6, max(ur) + 0.6]);
end
if hasR && ~isempty(ctr)
    ax = nexttile(tl);  hold(ax, 'on');
    for r = 1:numel(ur)
        i = Wt.rep == ur(r);
        scatter(ax, Wt.(ctr{1})(i), Wt.r_pred(i), 20, cr(r, :), 'filled', 'MarkerFaceAlpha', 0.7, ...
                'DisplayName', sprintf('rep %d', ur(r)));
    end
    if any(Wt.highlight)
        i = Wt.highlight;
        scatter(ax, Wt.(ctr{1})(i), Wt.r_pred(i), 30, [0.85 0.1 0.1], 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'highlighted');
    end
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    xlabel(ax, strrep(ctr{1}, '_', ' '));  ylabel(ax, 'r pred (within well)');
    legend(ax, 'Box', 'off', 'FontSize', 7, 'Location', 'best');
    title(ax, 'Does model failure follow mask alignment?', 'FontSize', 8.5);
    set(ax, 'FontSize', 8, 'Box', 'off', 'TickDir', 'out');
end
end

function i_save(f, o)
% <yyyy-mm-dd>__link-audit_wells.svg ; overwrite on re-run
if ~o.save || isempty(o.savepath) || isempty(f), return; end
d = char(o.savepath);  if ~exist(d, 'dir'), mkdir(d); end
fn = fullfile(d, sprintf('%s__link-audit_wells.svg', char(datetime('today', 'Format', 'yyyy-MM-dd'))));
print(f, fn, '-dsvg', '-painters');
fprintf('  saved %s\n', fn);
end