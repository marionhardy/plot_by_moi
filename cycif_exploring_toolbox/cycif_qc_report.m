function QC = cycif_qc_report(ip, op, dataloc, cycif_info, xy_list, varargin)
% CYCIF_QC_REPORT
% Quality-control report for the fixed-mode CycIF pipeline. Verifies, per XY and
% aggregated per replicate:
%   (1) MASK QC        — live-frame reference mask nucleus/cell counts & shape.
%   (2) REGISTRATION QC — per stain round: pre/post NN distance, improvement,
%                         transform magnitude, pass/fail.
%   (3) FOLLOWABILITY   — which reference cells are successfully measured in each
%                         round after registration, i.e. how many cells can be
%                         followed from the live last frame through to the last
%                         CycIF round. Views retained in the figures:
%                           - spatial map   (reference cells colored by
%                                            #well-registered rounds survived)
%                           - survival histogram (measured vs well-registered)
%                         (The presence matrix and followability curve panels
%                          were dropped by request; the underlying metrics are
%                          still computed and returned in QC.xy.)
%
% TOPOLOGY REMINDER: every round registers its DAPI to the SAME live-frame
% reference masks (hub-and-spoke). "Same cell across rounds" therefore means the
% same reference-mask cell (fixed cellindex) measured in each round — a cell is
% "followed through round r" iff it has a non-NaN measurement in round r.
%
% INPUTS
%   ip, op       celltracer structs (op unused here but kept for signature parity)
%   dataloc      fixed-mode dataloc AFTER ct_cycif_link_fixed
%   cycif_info   struct array (per round: .chan_map maps marker->channel). Used to
%                map measured antibody columns back to their round.
%   xy_list      [1×N] XY indices (default: all populated in dataloc.d)
%
% OPTIONS (name/value)
%   'save_dir'    char    directory to save figures (default '' = don't save)
%   'per_xy'      logical render one detail figure per XY (default true)
%   'summary'     logical render the replicate-level summary figure (default true)
%   'measured_marker' char  which compartment defines "measured" (default '_nuc')
%   'min_nuc_area' scalar  shape-QC area threshold for context (default 150)
%   'reg_nn_thresh' scalar well-registered gate on post-NN residual (default 10)
%   'crop_xys'    [1×N]  XYs to render the raw-DAPI registration overlay for
%                        (reopens each round's ND2; needs ip + staininfo).
%                        Empty -> no overlay.                     (default [])
%
% OUTPUT
%   QC  struct with per-XY and aggregate metrics:
%     QC.xy(i).xy, .n_masknuc, .n_cell, .n_linkable, .rounds_acquired,
%              .presence [nCells×nRounds logical], .survived [nCells×1],
%              .reg [nRounds struct: pre,post,improvement,tx,ty,rot,failed]
%     QC.agg   replicate-level rollups.

% ---- options ----
p = struct('save_dir','', 'per_xy',true, 'summary',true, ...
           'measured_marker','_nuc', 'min_nuc_area',150, 'reg_nn_thresh',10, ...
           'crop_xys',[]);   % XYs for the raw-DAPI registration overlay ([]=none)
for k = 1:2:numel(varargin)
    key = lower(varargin{k});
    if ~isfield(p, key); error('cycif_qc_report: unknown option ''%s''', varargin{k}); end
    p.(key) = varargin{k+1};
end
if nargin < 5 || isempty(xy_list)
    xy_list = find(cellfun(@(d) ~isempty(d) && isfield(d,'cellindex'), dataloc.d));
end
nRounds = numel(cycif_info);

% ---- map each round -> the antibody field STEMS actually stored ----------
% The linker names stored columns from the platemap marker strings (via
% local_parse_markers), which may NOT equal cycif_info.chan_map keys. So we
% derive the mapping from staininfo(r).markers (same source the linker used),
% falling back to chan_map only if staininfo is unavailable. Each stem is the
% makeValidName-sanitized antibody, matching how the linker built the field.
si = [];
if isfield(dataloc,'platemapd') && isfield(dataloc.platemapd,'staininfo')
    si = dataloc.platemapd.staininfo;
end
round_markers = cell(nRounds,1);
for r = 1:nRounds
    stems = {};
    if ~isempty(si) && r <= numel(si) && isfield(si(r),'markers')
        % union of antibody names across all wells in this round's marker grid
        cells = si(r).markers(:);
        for cix = 1:numel(cells)
            s = cells{cix};
            if isempty(s) || ~(ischar(s)||isstring(s)); continue; end
            % split co-markers on newline / CR / ';', take antibody (pre last '-')
            parts = strsplit(strtrim(char(s)), {char(10),char(13),';'});
            for pp = 1:numel(parts)
                e = strtrim(parts{pp});
                if isempty(e); continue; end
                hy = strsplit(e,'-');
                if numel(hy) < 2; continue; end
                ab = strtrim(strjoin(hy(1:end-1),'-'));
                stems{end+1} = matlab.lang.makeValidName(ab); %#ok<AGROW>
            end
        end
    elseif isfield(cycif_info(r),'chan_map') && isstruct(cycif_info(r).chan_map)
        stems = cellfun(@(k) matlab.lang.makeValidName(k), ...
                        fieldnames(cycif_info(r).chan_map), 'uni',0);
    end
    round_markers{r} = unique(stems);
    % Drop background/control pseudo-markers (not real antibodies; no stored
    % column). These are blank/background wells labeled e.g. 'bkgd', 'ctrl',
    % 'background' in the platemap.
    if ~isempty(round_markers{r})
        is_bkgd = ~cellfun('isempty', regexpi(round_markers{r}, ...
                            '^(bkgd|bkg|background|ctrl|control|blank)', 'once'));
        round_markers{r}(is_bkgd) = [];
    end
end

QC = struct('xy', struct([]), 'agg', struct());
comp = p.measured_marker;

% =====================================================================
% PER-XY PASS
% =====================================================================
for ii = 1:numel(xy_list)
    xy = xy_list(ii);
    d  = dataloc.d{xy};
    if isempty(d) || ~isfield(d,'data'); continue; end
    D  = d.data;
    ci = d.cellindex(:);
    nCells = numel(ci);

    % --- reference mask centroids + segmentation method (for spatial map) ---
    [cxy, n_masknuc, n_cell, seg_method] = local_mask_centroids(dataloc, xy, nCells);

    % --- registration QC pull (needed before followability gating) ---
    reg = local_pull_reg(d, nRounds);

    % Registration-quality gate: a round is "well-registered" if it improved
    % alignment AND the residual is within tolerance (~half a cell).
    well_reg = false(1, nRounds);
    for r = 1:nRounds
        well_reg(r) = (reg(r).improvement > 0) && (reg(r).post <= p.reg_nn_thresh);
    end

    % --- presence matrix: cell measured in round r? ---
    % A round counts as "acquired" only if at least one of its marker columns
    % is present AND carries at least one non-NaN value. A column that exists
    % but is entirely NaN means the linker wrote a scaffold for a round that
    % was never actually acquired (e.g. a pending stain round whose chan_map
    % was empty -- the platemap still lists a marker name so a stem resolves,
    % but no real measurement exists). Requiring real data keeps those pending
    % rounds OUT of acq_rounds so they don't pollute the registration/transform
    % panels (all-NaN tp) or count as "poor" rounds.
    presence = false(nCells, nRounds);
    acquired = false(1, nRounds);
    postmat  = nan(nCells, nRounds);   % per-cell round post-NN (for gradient)
    for r = 1:nRounds
        got = false(nCells,1);
        for mk = 1:numel(round_markers{r})
            fld = [round_markers{r}{mk} comp];   % stems already sanitized
            if isfield(D, fld)
                col_ok = ~isnan(D.(fld));
                if any(col_ok)
                    acquired(r) = true;
                    got = got | col_ok;
                end
            end
        end
        presence(:,r) = got;
        postmat(got, r) = reg(r).post;   % measured cells carry the round's post-NN
    end
    acq_rounds = find(acquired);

    % --- followability ---
    % PRIMARY (gap-tolerant): TOTAL number of acquired / well-registered rounds
    % a cell is measured in. One isolated bad round (e.g. a single failed
    % registration) no longer zeroes a cell that is fine in every other round,
    % because each round registers INDEPENDENTLY to the reference (hub-and-spoke)
    % -- a gap at round k does not corrupt rounds k+1..N.
    % SECONDARY (strict, consecutive from round 1): kept for the "unbroken chain"
    % view; stored but not plotted.
    only_acq        = presence(:, acq_rounds);
    goodmask        = presence & well_reg;
    survived        = sum(only_acq, 2);                      % total acquired rounds measured
    survived_good   = sum(goodmask(:, acq_rounds), 2);       % total well-registered rounds
    survived_consec      = local_run_length(presence,            acq_rounds);
    survived_good_consec = local_run_length(presence & well_reg, acq_rounds);

    % --- store metrics ---
    QC.xy(ii).xy              = xy;
    QC.xy(ii).n_masknuc       = n_masknuc;
    QC.xy(ii).n_cell          = n_cell;
    QC.xy(ii).n_linkable      = sum(ci>0);
    QC.xy(ii).rounds_acquired = acq_rounds;
    QC.xy(ii).presence        = presence;
    QC.xy(ii).survived        = survived;
    QC.xy(ii).survived_good   = survived_good;
    QC.xy(ii).survived_consec      = survived_consec;
    QC.xy(ii).survived_good_consec = survived_good_consec;
    QC.xy(ii).well_reg        = well_reg;
    QC.xy(ii).postmat         = postmat;
    QC.xy(ii).reg             = reg;
    QC.xy(ii).seg_method      = seg_method;

    if p.per_xy
        f = local_plot_xy(xy, presence, survived, survived_good, postmat, ...
                          cxy, reg, well_reg, acq_rounds, ...
                          n_masknuc, n_cell, sum(ci>0), p.reg_nn_thresh, seg_method);
        % filename encodes the segmentation method (Voronoi vs Cellpose)
        local_save(f, p.save_dir, sprintf('qc_xy%02d_%s', xy, lower(seg_method)));
    end

    % --- raw-DAPI registration overlay (only for requested crop_xys) ---
    % Separate figure that reopens each round's ND2 and draws the reference
    % centroids on the raw DAPI: green = pre-transform (live_xy), blue =
    % after registration (live_xformed). Blue landing on the round's nuclei
    % = good alignment. Needs ip (live frame) + staininfo (round ND2 paths).
    if ismember(xy, p.crop_xys)
        fo = local_align_overlay(ip, op, dataloc, xy, reg, acq_rounds);
        if ~isempty(fo)
            local_save(fo, p.save_dir, sprintf('qc_align_overlay_xy%02d', xy));
        end
    end
end

% =====================================================================
% REPLICATE-LEVEL SUMMARY
% =====================================================================
if ~isempty(QC.xy)
    % Always aggregate (cheap, and carries the registration-attrition flags
    % that the caller may want even when summary plotting is off).
    QC.agg = local_aggregate(QC.xy, nRounds);
    % Convenience top-level list of flagged XYs (registration attrition):
    % many rounds attempted, <50% passed -- the XY28 signature.
    if isfield(QC.agg,'pass_frac')
        fmask_flag = QC.agg.pass_frac < 0.5 & QC.agg.n_registered >= 3;
        QC.flagged_xys = QC.agg.xy_ids(fmask_flag);
    else
        QC.flagged_xys = [];
    end
    % dominant segmentation method across XYs (for title + filename)
    methods = {QC.xy.seg_method};
    umeth   = unique(methods);
    if isscalar(umeth), seg_tag = umeth{1};
    else,               seg_tag = strjoin(umeth, '+'); end   % mixed run
    if p.summary
        f = local_plot_summary(QC.agg, QC.xy, nRounds, seg_tag);
        local_save(f, p.save_dir, sprintf('qc_summary_%s', lower(seg_tag)));
    end
end

fprintf('[cycif_qc_report] Done: %d XYs.\n', numel(QC.xy));
end


% =====================================================================
% Helpers
% =====================================================================
function [cxy, n_masknuc, n_cell, seg_method] = local_mask_centroids(dataloc, xy, nCells)
% Reference nucleus centroids from the saved endpoint mask file, indexed by
% cellindex label so row i of cxy = centroid of cell i. Also returns the
% segmentation method (from the file's what_seg provenance). Falls back to NaN
% coords / 'unknown' method if the mask file is unavailable.
%
% cxy is ALWAYS [nCells x 2]: max(nuc_labels(:)) is not guaranteed to equal
% nCells (cellindex can have more entries than there are labeled nuclei -- e.g.
% label IDs with no surviving region after reconciliation), so any label ID
% beyond the mask's max up to nCells simply stays NaN. The spatial-survival
% panel filters NaN rows, so it renders partial coverage rather than treating a
% partial fill as a total failure.
cxy = nan(nCells,2); n_masknuc = NaN; n_cell = NaN; seg_method = 'unknown';
try
    fmask = fullfile(dataloc.exp.proc_dir, ...
        sprintf('%s_xy%02d_endpoint_masks.mat', dataloc.exp.name, xy));
    if ~exist(fmask,'file'); return; end   % no masks on disk; leave defaults
    M = load(fmask, 'nuc_labels', 'cyt_labels', 'what_seg');
    if isfield(M,'nuc_labels') && ~isempty(M.nuc_labels)
        st = regionprops(M.nuc_labels, 'Centroid');
        n_masknuc = numel(unique(M.nuc_labels(:))) - 1;
        for lbl = 1:min(numel(st), nCells)
            if ~isempty(st(lbl).Centroid), cxy(lbl,:) = st(lbl).Centroid; end
        end
    end
    if isfield(M,'cyt_labels') && ~isempty(M.cyt_labels)
        n_cell = numel(unique(M.cyt_labels(:))) - 1;
    end
    if isfield(M,'what_seg') && ~isempty(M.what_seg)
        seg_method = local_seg_label(M.what_seg);
    end
catch
    % leave defaults
end
end

function lbl = local_seg_label(what_seg)
% Map what_seg provenance to a concise segmentation-method label.
switch lower(char(what_seg))
    case {'cyt','nuc'}
        lbl = 'Voronoi';                 % legacy nuclei + Voronoi cytoplasm
    case {'whole','cytwhole'}
        lbl = 'Cellpose';                % cellpose whole-cell body
    otherwise
        lbl = char(what_seg);
end
end

function s = local_run_length(mask, acq_rounds)
% Per-cell count of consecutive TRUE rounds from the first acquired round.
nCells = size(mask,1);
s = zeros(nCells,1);
if isempty(acq_rounds); return; end
for cc = 1:nCells
    run = 0;
    for r = acq_rounds
        if mask(cc,r), run = run + 1; else, break; end
    end
    s(cc) = run;
end
end

function reg = local_pull_reg(d, nRounds)
% Extract per-round registration QC into a compact struct array.
% Recompute the post-registration NN distance in the DENSITY-ROBUST direction
% (moving DAPI -> nearest reference) from the stored point clouds, so existing
% dataloc objects (whose median_post_px was stored in the inflated ref->dapi
% direction) get corrected without re-running the linker. Falls back to the
% stored value if the clouds are unavailable.
reg = repmat(struct('pre',NaN,'post',NaN,'improvement',NaN, ...
                    'tx',NaN,'ty',NaN,'rot',NaN,'failed',true, ...
                    'live_xformed',[],'dapi_xy',[],'live_xy',[]), 1, nRounds);
if ~isfield(d,'cycif_reg'); return; end
R = d.cycif_reg;
for r = 1:min(numel(R), nRounds)
    q = R(r);
    pre = NaN; post = NaN;
    if all(isfield(q,{'live_xformed','dapi_xy','live_xy'})) ...
            && ~isempty(q.dapi_xy) && ~isempty(q.live_xformed)
        % density-robust: for each DAPI point, distance to nearest reference
        [~, dpost] = knnsearch(q.live_xformed, q.dapi_xy);
        post = median(dpost);
        if ~isempty(q.live_xy)
            [~, dpre] = knnsearch(q.live_xy, q.dapi_xy);
            pre = median(dpre);
        end
    end
    if isnan(post) && isfield(q,'median_post_px'), post = q.median_post_px; end
    if isnan(pre)  && isfield(q,'median_pre_px'),  pre  = q.median_pre_px;  end
    reg(r).pre  = pre;
    reg(r).post = post;
    reg(r).improvement = pre - post;
    if isfield(q,'tp') && numel(q.tp)>=3
        reg(r).tx = q.tp(1); reg(r).ty = q.tp(2); reg(r).rot = rad2deg(q.tp(3));
    end
    if isfield(q,'failed'), reg(r).failed = logical(q.failed); end
    % carry clouds for the overlay panels (empty if unavailable)
    if isfield(q,'live_xformed'), reg(r).live_xformed = q.live_xformed; else, reg(r).live_xformed = []; end
    if isfield(q,'dapi_xy'),      reg(r).dapi_xy      = q.dapi_xy;      else, reg(r).dapi_xy = []; end
    if isfield(q,'live_xy'),      reg(r).live_xy      = q.live_xy;      else, reg(r).live_xy = []; end
end
end

function f = local_align_overlay(ip, op, dataloc, xy, reg, acq_rounds)
% RAW-DAPI registration overlay for a single XY. Reopens each round's ND2 and
% overlays the reference nucleus centroids on the raw DAPI frame:
%   green = live_xy      (reference cloud, pre-transform)
%   blue  = live_xformed (reference cloud, after registration)
% Blue landing on the round's nuclei = good alignment; blue off the nuclei =
% that round mis-registered. The leftmost tile is the LIVE last frame with the
% reference centroids, so a nucleus the rounds are "missing" can be checked
% against the reference it registers to. Returns [] (no figure) if the inputs
% needed to load frames are unavailable.
f = [];
if isempty(op) || isempty(ip) || ~isfield(ip,'fname')
    warning('cycif_qc_report: align overlay needs ip.fname + op -- skipping XY%02d.', xy);
    return
end
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'staininfo')
    warning('cycif_qc_report: no staininfo -- cannot resolve round ND2s for overlay (XY%02d).', xy);
    return
end
si = dataloc.platemapd.staininfo;
dapi_chan = 1;

% Live nuclear channel resolver (the reference is the LIVE last frame).
if isfield(op,'seg') && isfield(op.seg,'provider') && strcmpi(op.seg.provider,'cellpose') ...
        && isfield(op.seg,'cp') && isfield(op.seg.cp,'chan_nuc') && ~isempty(op.seg.cp.chan_nuc)
    live_chan = op.seg.cp.chan_nuc;
elseif isfield(op,'seg') && isfield(op.seg,'chan') && ~isempty(op.seg.chan)
    plane = op.seg.chan(1);
    if isfield(op,'cind') && ~isempty(op.cind) && plane<=numel(op.cind); live_chan=op.cind(plane); else; live_chan=plane; end
else
    live_chan = 1;
end

% Choose up to 4 acquired rounds that actually have stored clouds, evenly spaced.
have_cloud = acq_rounds(arrayfun(@(r) ~isempty(reg(r).live_xformed), acq_rounds));
if isempty(have_cloud)
    warning('cycif_qc_report: no registration clouds stored for XY%02d -- skipping overlay.', xy);
    return
end
sel = have_cloud(round(linspace(1, numel(have_cloud), min(4,numel(have_cloud)))));

f = figure('Color','w','Name',sprintf('QC align overlay XY%02d',xy), ...
           'Position',[40 40 1450 340]);
tl = tiledlayout(f, 1, numel(sel)+1, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf(['XY%02d reg overlay: leftmost = LIVE last frame (reference) with its ' ...
                   'centroids; then per round green=reference blue=after-registration'], xy));

% --- reference tile: LIVE last frame with the reference centroids ---
axr = nexttile(tl);
imgr = []; lt = 1;
fmask = fullfile(dataloc.exp.proc_dir, ...
    sprintf('%s_xy%02d_endpoint_masks.mat', dataloc.exp.name, xy));
if exist(fmask,'file')
    Lt = load(fmask,'last_t');
    if isfield(Lt,'last_t') && ~isempty(Lt.last_t); lt = Lt.last_t; end
end
try
    dao_live = ImageAccess(ip.fname, ip.indsz);
    imgr = double(dao_live.GetFrame([lt, live_chan, xy, 1]));
catch
end
if isempty(imgr)
    axis(axr,'off'); text(0.1,0.5,'live frame n/a','Parent',axr);
else
    lo=prctile(imgr(:),2); hi=prctile(imgr(:),99.5);
    imagesc(axr,imgr,[lo hi]); colormap(axr,gray); axis(axr,'image','off'); hold(axr,'on');
    lxy0 = reg(sel(1)).live_xy;
    if ~isempty(lxy0)
        plot(axr,lxy0(:,1),lxy0(:,2),'o','Color',[0.2 0.9 0.3],'MarkerSize',5,'LineWidth',0.7);
    end
    title(axr,'REFERENCE (live last frame)','FontSize',9);
end

% --- per-round tiles: raw DAPI + green (pre) / blue (post) reference clouds ---
for k = 1:numel(sel)
    r = sel(k);
    ax = nexttile(tl);
    img = [];
    if r <= numel(si) && ~isempty(si(r).nd2_path)
        try
            dao = iman_imageaccess(si(r).nd2_path);
            img = double(iman_getframe(dao,[1,dapi_chan,xy,1]));
        catch
        end
    end
    if isempty(img)
        axis(ax,'off'); text(0.1,0.5,sprintf('round %d n/a',r),'Parent',ax); continue
    end
    lo=prctile(img(:),2); hi=prctile(img(:),99.5);
    imagesc(ax,img,[lo hi]); colormap(ax,gray); axis(ax,'image','off'); hold(ax,'on');
    if ~isempty(reg(r).live_xy)
        plot(ax,reg(r).live_xy(:,1),reg(r).live_xy(:,2),'o', ...
             'Color',[0.2 0.9 0.3],'MarkerSize',5,'LineWidth',0.7);
    end
    if ~isempty(reg(r).live_xformed)
        plot(ax,reg(r).live_xformed(:,1),reg(r).live_xformed(:,2),'o', ...
             'Color',[0.2 0.5 1],'MarkerSize',5,'LineWidth',0.7);
    end
    title(ax,sprintf('round %d (%.0fpx)', r, reg(r).post),'FontSize',9);
end
end

function f = local_plot_xy(xy, presence, survived, survived_good, postmat, ...
                           cxy, reg, well_reg, acq_rounds, ...
                           n_masknuc, n_cell, n_linkable, nn_thresh, seg_method)
% Per-XY detail (2×2, 4 panels): spatial survival map, survival histogram,
% flagged registration QC, transform magnitude. Title reports the segmentation
% method (Voronoi/Cellpose). (Presence matrix and followability curve dropped
% by request; a fifth registration-overlay panel is also omitted here.)
if nargin < 14 || isempty(seg_method), seg_method = 'unknown'; end
FS = 10;  % base font size
f = figure('Color','w','Position',[50 50 1500 900], ...
           'Name',sprintf('CyCIF QC — XY%02d [%s]', xy, seg_method));
tl = tiledlayout(f, 2, 2, 'Padding','compact','TileSpacing','compact');
title(tl, sprintf('CyCIF QC — XY%02d   ·   seg: %s   ·   %d mask nuclei · %d linkable', ...
      xy, seg_method, n_masknuc, n_linkable), 'FontWeight','bold','FontSize',13);

rr = acq_rounds; nAcq = numel(rr);

% (1) SPATIAL MAP — reference cells colored by well-registered survival.
%     cxy is [nCells x 2] with NaN rows for cellindex slots lacking a labeled
%     nucleus; those are filtered here so the panel renders whatever coverage
%     exists rather than going all-or-nothing on a size mismatch.
ax = nexttile(tl);
valid = ~any(isnan(cxy),2);
if any(valid)
    scatter(ax, cxy(valid,1), cxy(valid,2), 16, survived_good(valid), 'filled');
    colormap(ax, parula); ax.CLim = [0 max(nAcq,1)];
    cb = colorbar(ax); cb.Label.String = '#well-reg rounds survived'; cb.FontSize = FS-2;
    axis(ax,'image'); ax.YDir='reverse';
else
    text(ax,0.5,0.5,'no mask centroids','Units','normalized','HorizontalAlignment','center');
end
xlabel(ax,'X (px)','FontSize',FS); ylabel(ax,'Y (px)','FontSize',FS);
ax.FontSize = FS-1; title(ax,'Spatial survival (well-registered)','FontSize',FS);

% (2) SURVIVAL HISTOGRAM — measured vs well-registered side by side.
%     NOTE: this uses a plain grouped bar() over integer survival counts (0..smax)
%     -- NOT histogram(). MATLAB 2022's histogram() throws a recurring internal
%     error ("Operands to && must be convertible to logical scalar" at its
%     factor==-2||... auto-bin-count heuristic) regardless of calling form, so
%     any histogram of integer/discrete data in this file is built with bar().
ax = nexttile(tl);
smax = max([survived(:); survived_good(:); 1]);
ctrs = 0:double(smax);
c1 = arrayfun(@(v) sum(survived==v),      ctrs);
c2 = arrayfun(@(v) sum(survived_good==v), ctrs);
bar(ax, ctrs, [c1(:) c2(:)], 1.0, 'grouped');
if numel(ax.Children) >= 2
    ax.Children(2).FaceColor = [0.55 0.6 0.75];   % measured
    ax.Children(1).FaceColor = [0.15 0.45 0.30];  % well-reg
end
legend(ax, {'measured','well-registered'}, 'FontSize',FS-2,'Box','off','Location','north');
xlabel(ax,'#rounds survived per cell','FontSize',FS); ylabel(ax,'# cells','FontSize',FS);
xlim(ax,[-0.5 double(smax)+0.5]); ax.FontSize = FS-1;
title(ax, sprintf('Mask: %d nuclei / %d cell bodies', n_masknuc, n_cell),'FontSize',FS);
grid(ax,'on');

% (3) REGISTRATION QC — post-NN per round with threshold line; poor rounds
%     (failed OR not well-registered) marked with a red x.
ax = nexttile(tl);
post = arrayfun(@(r) reg(r).post, rr);
impr = arrayfun(@(r) reg(r).improvement, rr);
poor = ~well_reg(rr);
if nAcq >= 1
    yyaxis(ax,'left');
    plot(ax, 1:nAcq, post, '-o','LineWidth',1.5,'Color',[0.2 0.4 0.7], ...
         'MarkerFaceColor',[0.2 0.4 0.7]); hold(ax,'on');
    yline(ax, nn_thresh, '--','Color',[0.85 0.3 0.3], ...
          'Label',sprintf('%gpx gate',nn_thresh),'FontSize',FS-3);
    bp = find(poor);
    if ~isempty(bp)
        plot(ax, bp, post(bp), 'x','Color',[0.85 0.15 0.15],'MarkerSize',11,'LineWidth',2.2);
    end
    ylabel(ax,'post-reg NN dist (px)','FontSize',FS); ax.YColor = [0.2 0.4 0.7];
    yyaxis(ax,'right');
    bar(ax, 1:nAcq, impr, 0.4, 'FaceColor',[0.9 0.5 0.2],'FaceAlpha',0.35,'EdgeColor','none');
    ylabel(ax,'improvement (px)','FontSize',FS); ax.YColor = [0.7 0.4 0.1];
    ax.XTick = 1:nAcq; ax.XTickLabel = rr; ax.FontSize = FS-1;
else
    text(ax,0.5,0.5,'no acquired rounds','Units','normalized','HorizontalAlignment','center');
end
xlabel(ax,'stain round','FontSize',FS);
title(ax, sprintf('Registration QC — %d/%d rounds poor', sum(poor), nAcq),'FontSize',FS);

% (4) TRANSFORM MAGNITUDE per round.
ax = nexttile(tl);
txv = arrayfun(@(r) hypot(reg(r).tx,reg(r).ty), rr);
rotv= arrayfun(@(r) abs(reg(r).rot), rr);
if nAcq >= 1 && any(~isnan(txv))
    yyaxis(ax,'left');  plot(ax,1:nAcq,txv,'-o','LineWidth',1.5); ylabel(ax,'|translation| (px)','FontSize',FS);
    yyaxis(ax,'right'); plot(ax,1:nAcq,rotv,'-s','LineWidth',1.5); ylabel(ax,'|rotation| (deg)','FontSize',FS);
    ax.XTick = 1:nAcq; ax.XTickLabel = rr; ax.FontSize = FS-1;
else
    axis(ax,'off'); text(ax,0.5,0.5,'transform magnitude n/a','Units','normalized','HorizontalAlignment','center');
end
xlabel(ax,'stain round','FontSize',FS); title(ax,'Transform magnitude','FontSize',FS);
end

function agg = local_aggregate(xyQC, nRounds)
% Roll up per-XY metrics to replicate level.
agg = struct();
agg.n_xy        = numel(xyQC);
agg.masknuc     = arrayfun(@(x) x.n_masknuc, xyQC);
agg.linkable    = arrayfun(@(x) x.n_linkable, xyQC);
% per-round measured fraction + post-NN (mean across XYs)
frac = nan(numel(xyQC), nRounds);
postd = nan(numel(xyQC), nRounds);
poorrate = zeros(1, nRounds); poorden = zeros(1, nRounds);
% Per-XY registration attrition: how many rounds were REGISTERED (tp populated,
% i.e. the linker ran registration for that round) vs how many PASSED (became
% acquired = registered well enough that measurement was written). A large gap
% means the linker attempted many rounds for this XY but failed most of them --
% the XY28 signature (huge stage drift -> geometric-outlier auto-fail -> no
% column written -> round silently absent from the per-XY figure's x-axis).
% Rounds whose tp is all-NaN are genuinely pending/never-run (empty chan_map),
% NOT failures, so they are excluded from the "attempted" count.
n_registered = zeros(1, numel(xyQC));   % rounds with a real tp (attempted)
n_passed     = zeros(1, numel(xyQC));   % rounds that became acquired
n_failed     = zeros(1, numel(xyQC));   % registered but failed (attempted - passed-among-registered)
for i = 1:numel(xyQC)
    P = xyQC(i).presence;
    frac(i,:) = mean(P,1);
    for r = 1:nRounds
        postd(i,r) = xyQC(i).reg(r).post;
        if ismember(r, xyQC(i).rounds_acquired)
            poorden(r)  = poorden(r) + 1;
            poorrate(r) = poorrate(r) + double(~xyQC(i).well_reg(r));
        end
    end
    % attempted = tp populated (tx non-NaN); failed flag comes straight from
    % the linker's own geometric-outlier / improvement gate.
    tp_run = arrayfun(@(rr) ~isnan(xyQC(i).reg(rr).tx), 1:nRounds);
    fl     = arrayfun(@(rr) logical(xyQC(i).reg(rr).failed), 1:nRounds);
    n_registered(i) = sum(tp_run);
    n_passed(i)     = numel(xyQC(i).rounds_acquired);
    n_failed(i)     = sum(tp_run & fl);
end
agg.round_measured_frac = frac;
agg.round_post_nn       = postd;
agg.round_poor_frac     = poorrate ./ max(poorden,1);   % frac of XYs poorly-reg per round
agg.xy_ids           = arrayfun(@(x) x.xy, xyQC);
agg.n_registered     = n_registered;   % [1 x nXY] rounds with real tp
agg.n_passed         = n_passed;        % [1 x nXY] rounds acquired
agg.n_failed         = n_failed;        % [1 x nXY] registered-but-failed rounds
agg.pass_frac        = n_passed ./ max(n_registered,1);   % fraction of attempted rounds that passed
% overall followability per XY: measured vs well-registered
agg.frac_survive_all = arrayfun(@(x) ...
    mean(x.survived      >= numel(x.rounds_acquired)), xyQC);
agg.frac_survive_good = arrayfun(@(x) ...
    mean(x.survived_good >= numel(x.rounds_acquired)), xyQC);

% Console flag for problem XYs: many rounds attempted, few passed. Surfaces
% XY28-type failures without needing to eyeball each per-XY figure. Threshold
% is deliberately loose (pass_frac < 0.5 AND at least 3 rounds attempted) so it
% catches genuine attrition without flagging XYs that only ever had 1-2 rounds.
flagged = find(agg.pass_frac < 0.5 & agg.n_registered >= 3);
if ~isempty(flagged)
    fprintf('\n[cycif_qc_report] REGISTRATION-ATTRITION WARNING: %d XY(s) passed <50%% of attempted rounds:\n', numel(flagged));
    for j = flagged
        fprintf('   XY%02d: %d/%d rounds passed (%d failed registration), %d mask nuclei\n', ...
            agg.xy_ids(j), agg.n_passed(j), agg.n_registered(j), ...
            agg.n_failed(j), agg.masknuc(j));
    end
    fprintf('   -> large stage drift or sparse/low-quality fields; consider excluding.\n');
end
end

function f = local_plot_summary(agg, xyQC, nRounds, seg_tag)
% Replicate-level summary (2×3, 5 panels): per-round registration quality, mask
% nucleus counts, whole-panel followability by XY, linkable vs mask nuclei, and
% per-XY registration attrition (rounds passed vs attempted -- the panel that
% surfaces XY28-type failures where many rounds registered but few passed).
if nargin < 4 || isempty(seg_tag), seg_tag = 'unknown'; end
FS = 10;
f = figure('Color','w','Position',[60 60 1500 880], ...
           'Name',sprintf('CyCIF QC — replicate summary [%s]', seg_tag));
tl = tiledlayout(f,2,3,'Padding','compact','TileSpacing','compact');
title(tl, sprintf('CyCIF QC — replicate summary (%d XYs)   ·   seg: %s', ...
      agg.n_xy, seg_tag), 'FontWeight','bold','FontSize',13);

% (1) per-round registration quality: post-NN + fraction of XYs poorly-reg
ax = nexttile(tl);
mp = mean(agg.round_post_nn,1,'omitnan');
sp = std(agg.round_post_nn,0,1,'omitnan');
yyaxis(ax,'left');
errorbar(ax, 1:nRounds, mp, sp, '-o','LineWidth',1.5,'Color',[0.7 0.3 0.3], ...
         'MarkerFaceColor',[0.7 0.3 0.3]);
ylabel(ax,'post-reg NN dist (px)','FontSize',FS); ax.YColor=[0.7 0.3 0.3];
yyaxis(ax,'right');
bar(ax, 1:nRounds, agg.round_poor_frac, 0.5, 'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor','none');
ylabel(ax,'frac XYs poorly-registered','FontSize',FS); ylim(ax,[0 1]); ax.YColor=[0.4 0.4 0.4];
xlabel(ax,'stain round','FontSize',FS); ax.FontSize=FS-1;
title(ax,'Per-round registration quality','FontSize',FS);

% (2) mask nucleus count distribution
ax = nexttile(tl);
mn = agg.masknuc(~isnan(agg.masknuc));
if isempty(mn)
    text(ax,0.5,0.5,'no data','Units','normalized','HorizontalAlignment','center');
else
    local_hist_bars(ax, mn, [0.4 0.5 0.7]);
end
xlabel(ax,'mask nuclei per XY','FontSize',FS); ylabel(ax,'# XYs','FontSize',FS);
ax.FontSize=FS-1; title(ax,'Mask nucleus counts','FontSize',FS); grid(ax,'on');

% (3) whole-panel followability per XY — measured vs well-registered
ax = nexttile(tl);
xyids = arrayfun(@(x)x.xy,xyQC);
bar(ax, xyids, [agg.frac_survive_all(:) agg.frac_survive_good(:)], 1.0, 'grouped');
if numel(ax.Children) >= 2
    ax.Children(2).FaceColor = [0.55 0.6 0.75];
    ax.Children(1).FaceColor = [0.15 0.45 0.30];
end
legend(ax,{'measured','well-registered'},'FontSize',FS-2,'Box','off','Location','southoutside','Orientation','horizontal');
xlabel(ax,'XY','FontSize',FS); ylabel(ax,'frac surviving all rounds','FontSize',FS);
ylim(ax,[0 1]); ax.FontSize=FS-1; title(ax,'Whole-panel followability by XY','FontSize',FS); grid(ax,'on');

% (4) linkable vs mask nuclei
ax = nexttile(tl);
scatter(ax, agg.masknuc, agg.linkable, 32, [0.2 0.4 0.7], 'filled','MarkerFaceAlpha',0.7);
hold(ax,'on'); mx = max(agg.masknuc); plot(ax,[0 mx],[0 mx],'k--');
xlabel(ax,'mask nuclei','FontSize',FS); ylabel(ax,'linkable cells','FontSize',FS);
axis(ax,'equal'); ax.FontSize=FS-1; title(ax,'Linkable vs mask nuclei','FontSize',FS); grid(ax,'on');

% (5) REGISTRATION ATTRITION per XY — rounds passed vs attempted. An XY where
%     many rounds registered (attempted) but few passed is the XY28 signature:
%     large stage drift auto-failing rounds so their columns are never written
%     and they silently vanish from that XY's per-XY figure x-axis. Flagged XYs
%     (pass_frac < 0.5 with >=3 attempted) are drawn in red so they stand out.
ax = nexttile(tl);
if isfield(agg,'n_registered') && ~isempty(agg.n_registered)
    xids = agg.xy_ids;
    flagged = agg.pass_frac < 0.5 & agg.n_registered >= 3;
    hold(ax,'on');
    % attempted as a light backdrop bar, passed as the solid foreground bar
    bar(ax, xids, agg.n_registered, 0.8, 'FaceColor',[0.85 0.85 0.85], 'EdgeColor','none', 'DisplayName','attempted');
    bp = bar(ax, xids, agg.n_passed, 0.8, 'FaceColor',[0.35 0.55 0.75], 'EdgeColor','none', 'DisplayName','passed');
    % recolor flagged XYs' passed bars to red via a second overlay bar
    if any(flagged)
        pv = agg.n_passed; pv(~flagged) = 0;
        br = bar(ax, xids, pv, 0.8, 'FaceColor',[0.85 0.2 0.2], 'EdgeColor','none', 'DisplayName','flagged (<50% passed)');
    end
    xlabel(ax,'XY','FontSize',FS); ylabel(ax,'# rounds','FontSize',FS);
    legend(ax,'Location','northoutside','Orientation','horizontal','FontSize',FS-3,'Box','off');
    ax.FontSize=FS-1; grid(ax,'on');
    nflag = sum(flagged);
    title(ax, sprintf('Registration attrition — %d XY(s) flagged', nflag),'FontSize',FS);
else
    axis(ax,'off'); text(ax,0.5,0.5,'attrition data n/a','Units','normalized','HorizontalAlignment','center');
end
end

function local_hist_bars(ax, v, facecolor)
% bar()-based histogram (see local_plot_xy panel 2 note for why histogram() is
% avoided on MATLAB 2022). Data-driven edges; guards the degenerate all-
% identical-value case with a single unit-width bin.
v = double(v(:)); v = v(isfinite(v));
if isempty(v); return; end
lo = min(v); hi = max(v);
if lo == hi
    edges = [lo-0.5, lo+0.5];
else
    nb    = max(5, ceil(sqrt(numel(v))));
    edges = linspace(lo, hi, nb+1);
end
counts  = histcounts(v, edges);
centers = edges(1:end-1) + diff(edges)/2;
bar(ax, centers, counts, 1.0, 'FaceColor', facecolor);
end

function local_save(f, save_dir, stem)
if isempty(save_dir); return; end
if ~exist(save_dir,'dir'); mkdir(save_dir); end
try
    print(f, fullfile(save_dir,[stem '.png']), '-dpng','-r150');
    fprintf('  saved %s.png\n', stem);
catch ME
    warning('cycif_qc_report: could not save %s (%s)', stem, ME.message);
end
end