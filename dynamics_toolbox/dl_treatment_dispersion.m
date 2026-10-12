function R = dl_treatment_dispersion(T, varargin)
% dl_treatment_dispersion  Does treatment define functional state?
%   Partitions CyCIF functional-marker variance into treatment,
%   well-within-treatment and cell components, compares within-group spread
%   to between-centroid separation, and interprets every value in the console.
%
% USAGE
%   R = dl_treatment_dispersion(Tf, 'prenormalized',true);                     % pooled, per-rep z-scored
%   R = dl_treatment_dispersion(T,  'transform','arcsinh', 'dataloc',dataloc); % raw cycif_build_table
%   R = dl_treatment_dispersion(Tf, 'prenormalized',true, 'nan_policy','panel', 'save',false);
%   R = dl_treatment_dispersion(Tf, 'prenormalized',true, 'markers',{'HLA_DR','CD86','CD206'}, 'drop_pcs',[]);
%   R = dl_treatment_dispersion(Tf, 'prenormalized',true, 'verbose','full', 'plots','all');   % everything
%
% INPUT
%   T   table, one row per cell (cycif_build_table output or pooled Tf)
%
% NAME-VALUE
%  Columns
%   'group'         char     grouping column                                ('tx1_label')
%   'unit'          char     unit treatment is assigned to; permutation level ('xy')
%                            combined with the replicate column when present
%   'markers'       cell     marker columns; {} = auto (numeric, non-ID)    ({})
%   'id_cols'       cell     numeric columns never treated as markers
%                            ({'xy','cellid','ptx_is_spike','rep','ExpID','xys','cellidx'})
%   'compartment'   char     '' | 'cyt' | 'nuc' | 'cell'; selects <marker>_<comp>
%                            columns. cycif_build_table strips the suffix, so
%                            leave '' for its output                        ('')
%   'rep_col'       char     replicate column; '' = auto {'rep','ExpID','exp'} ('')
%  Scale
%   'prenormalized' logical  T already per-replicate z-scored -> 'linear'  (false)
%   'transform'     char     'arcsinh' | 'log2' | 'linear'; '' = auto       ('')
%                            auto: 'linear' if prenormalized, else 'arcsinh'
%   'cofactor'      'auto' | scalar | [1xM]  arcsinh cofactor;
%                            'auto' = per-marker 10th pct of positive values ('auto')
%   'drop_pcs'      vector   PCs projected out before multivariate dispersion ([1])
%  Filtering
%   'exclude_spike' logical  drop ptx_is_spike rows                         (true)
%   'nan_policy'    char     'block' | 'panel' | 'complete'                 ('block')
%                            block    = largest co-stained marker x cell block,
%                                       preferring blocks covering all groups
%                            panel    = one analysis per exact stain signature
%                            complete = rows with every selected marker
%   'min_n'         scalar   warn on groups below this n; min panel size     (30)
%  Statistics
%   'nperm'         scalar   label permutations                             (999)
%   'eta_lo'        scalar   treatment fraction below which "does not define" (0.10)
%   'eta_hi'        scalar   treatment fraction above which "largely defines" (0.50)
%   'seed'          scalar   local RandStream seed (global rng untouched)    (0)
%  Output
%   'save'          logical  save SVGs                                       (true)
%   'closefigs'     logical  close figures after save (only if save)         (true)
%   'dataloc'       struct   fig dir from dataloc.fold.fig, else dataloc.exp.fig_dir ([])
%   'savepath'      char     fallback folder if no dataloc                   ('')
%   'max_plot_pts'  scalar   per-group subsample for point plots             (2000)
%   'plots'         'summary' = one step-by-step figure | 'all' = + diagnostics | 'none' ('summary')
%   'verbose'       'summary' = compact table + verdict | 'full' = every value interpreted | 'none' ('summary')
%
% OUTPUT
%   R  struct array (numel > 1 only for 'panel'), fields:
%     .set_tag .row_idx (rows of input T) .markers .n .groups
%     .transform_label .cofactor .drop_pcs .var_expl .pc1_uniform_sign
%     .partition   table  component | SS | frac
%     .eta2        table  marker | eta2 | eta2_unit (well|treatment)
%     .group       table  group | n | n_units | spread | nearest_dist | nearest_group | ratio
%     .ratio_global        mean within-group spread / mean centroid distance
%     .centroid_dist       [k x k] distances between group centroids (drop_pcs removed)
%     .ratio_global_keepPC same, without drop_pcs (sensitivity)
%     .F .p_perm .nperm .n_distinct .perm_level .null_frac (treatment fraction under each shuffle)
%     .verdict     'not-defining' | 'biasing' | 'defining'
%     .figs        figure handles (empty if closed)
%
% NOTES
%   New in DatalocHandler v5 (no v4 counterpart). Statistics are computed on
%   per-marker standardized values, so every marker weighs equally.

%% ---- Parse inputs -------------------------------------------------------
p = inputParser; p.FunctionName = mfilename;
isTxt = @(x) ischar(x) || isstring(x);
addRequired (p, 'T', @istable);
addParameter(p, 'group',         'tx1_label', isTxt);
addParameter(p, 'unit',          'xy',        isTxt);
addParameter(p, 'markers',       {},          @(x) iscell(x) || isstring(x));
addParameter(p, 'id_cols',       {'xy','cellid','ptx_is_spike','rep','ExpID','xys','cellidx'});
addParameter(p, 'compartment',   '',          isTxt);
addParameter(p, 'rep_col',       '',          isTxt);
addParameter(p, 'prenormalized', false,       @(x) islogical(x) || isnumeric(x));
addParameter(p, 'transform',     '',          isTxt);
addParameter(p, 'cofactor',      'auto');
addParameter(p, 'drop_pcs',      1,           @isnumeric);
addParameter(p, 'exclude_spike', true);
addParameter(p, 'nan_policy',    'block',     @(x) any(strcmpi(x, {'block','panel','complete'})));
addParameter(p, 'min_n',         30);
addParameter(p, 'nperm',         999);
addParameter(p, 'eta_lo',        0.10);
addParameter(p, 'eta_hi',        0.50);
addParameter(p, 'seed',          0);
addParameter(p, 'save',          true);
addParameter(p, 'closefigs',     true);
addParameter(p, 'dataloc',       []);
addParameter(p, 'savepath',      '',          isTxt);
addParameter(p, 'max_plot_pts',  2000);
addParameter(p, 'plots',         'summary',   @(x) any(strcmpi(x, {'summary','all','none'})));
addParameter(p, 'verbose',       'summary',   @(x) any(strcmpi(x, {'summary','full','none'})));
parse(p, T, varargin{:});
o = p.Results;
o.group       = char(o.group);       o.unit       = char(o.unit);
o.compartment = char(o.compartment); o.rep_col    = char(o.rep_col);
o.transform   = lower(char(o.transform));
o.nan_policy  = lower(char(o.nan_policy));
o.markers     = cellstr(o.markers);
o.prenormalized = logical(o.prenormalized);
o.plots = lower(char(o.plots));
o.vb = find(strcmpi(o.verbose, {'none','summary','full'})) - 1;   % 0 none, 1 summary, 2 full
rs  = RandStream('mt19937ar', 'Seed', o.seed);   % local stream
tag = '[dl_treatment_dispersion]';
vn  = T.Properties.VariableNames;

assert(ismember(o.group, vn), '%s group column ''%s'' not in table. Columns: %s', ...
       tag, o.group, strjoin(vn, ', '));

%% ---- Row filtering ------------------------------------------------------
keep = true(height(T), 1);
if o.exclude_spike && ismember('ptx_is_spike', vn)
    ns = nnz(T.ptx_is_spike);
    if ns > 0 && o.vb >= 2, fprintf('%s excluding %d spike-treatment cells (ptx_is_spike).\n', tag, ns); end
    keep = keep & ~logical(T.ptx_is_spike);
end
g   = string(T.(o.group));
bad = ismissing(g) | strlength(strtrim(g)) == 0;
if any(bad & keep) && o.vb >= 2, fprintf('%s dropping %d cells with empty %s.\n', tag, nnz(bad & keep), o.group); end
keep = keep & ~bad;

%% ---- Replicate + unit keys ---------------------------------------------
rep_col = o.rep_col;
if isempty(rep_col)
    cand = {'rep','ExpID','exp'};
    hit  = cand(ismember(cand, vn));
    if ~isempty(hit), rep_col = hit{1}; end
end
if ~isempty(rep_col)
    assert(ismember(rep_col, vn), '%s rep_col ''%s'' not in table.', tag, rep_col);
    nrep = numel(unique(string(T.(rep_col)(keep))));
    if nrep > 1 && ~o.prenormalized
        warning('%s %d replicates pooled without prenormalization: batch offsets inflate within-group spread (biases toward "treatment does not define state"). Prefer the per-rep z-scored table with ''prenormalized'',true.', tag, nrep);
    end
end

u = [];
if ~isempty(o.unit) && ismember(o.unit, vn)
    u = string(T.(o.unit));
    % xy numbers repeat across replicates: key units by rep|xy
    if ~isempty(rep_col), u = string(T.(rep_col)) + "|" + u; end
else
    warning('%s unit column ''%s'' absent: permutation falls back to cells (pseudo-replicated; p anti-conservative).', tag, o.unit);
end

%% ---- Markers + scale ----------------------------------------------------
[mk, mk_disp] = i_resolve_markers(T, o, rep_col, tag);

if o.prenormalized
    if ~isempty(o.transform) && ~strcmp(o.transform, 'linear')
        warning('%s prenormalized=true: ''%s'' on z-scores is invalid (negatives). Forcing ''linear''.', tag, o.transform);
    end
    xf = 'linear';
elseif isempty(o.transform)
    xf = 'arcsinh';
else
    xf = o.transform;
end

X       = double(T{keep, mk});
row_idx = find(keep);
gk      = g(keep);
uk      = [];  if ~isempty(u), uk = u(keep); end

% Drop markers with no data in the kept rows
allnan = all(isnan(X), 1);
if any(allnan)
    if o.vb >= 1, fprintf('%s dropping all-NaN markers: %s\n', tag, strjoin(mk_disp(allnan), ', ')); end
    X(:, allnan) = []; mk_disp(allnan) = [];
end
assert(size(X,2) >= 2, '%s fewer than 2 usable markers.', tag);

%% ---- NaN policy: define analysis sets -----------------------------------
S = i_nan_sets(X, mk_disp, gk, o.nan_policy, o.min_n, tag, o.vb);

%% ---- Analyze each set ---------------------------------------------------
fig_dir = i_figdir(o);
if o.save && isempty(fig_dir)
    warning('%s save=true but no dataloc fig dir or savepath: figures not saved.', tag);
end
if o.vb >= 2
    fprintf('%s eta_lo/eta_hi = %.2f/%.2f are interpretive conventions, not derived thresholds.\n', ...
            tag, o.eta_lo, o.eta_hi);
end

Rc = {};
all_groups = unique(gk);
for s = 1:numel(S)
    rows = S(s).rows;  cols = S(s).cols;
    gs   = gk(rows);
    miss = setdiff(all_groups, unique(gs));
    if ~isempty(miss)
        warning('%s [%s] groups absent from this marker set: %s', tag, S(s).tag, strjoin(miss, ', '));
    end
    if numel(unique(gs)) < 2
        warning('%s [%s] <2 groups; skipped.', tag, S(s).tag); continue
    end
    us = [];  if ~isempty(uk), us = uk(rows); end

    [Rk, P] = i_analyze(X(rows, cols), gs, us, mk_disp(cols), xf, o, rs, tag);
    Rk.set_tag = S(s).tag;
    Rk.row_idx = row_idx(rows);

    if o.vb == 2, i_report(Rk, o, tag); elseif o.vb == 1, i_brief(Rk, o, tag); end
    if ~strcmp(o.plots, 'none')
        Rk.figs = i_story(Rk, P, o, rs, fig_dir);
        if strcmp(o.plots, 'all'), Rk.figs = [Rk.figs, i_plots(Rk, P, o, rs, fig_dir)]; end
        if o.closefigs && o.save && ~isempty(fig_dir)
            close(Rk.figs);  Rk.figs = gobjects(0);
        end
    end
    Rc{end+1} = orderfields(Rk); %#ok<AGROW>
end
assert(~isempty(Rc), '%s no analyzable marker set.', tag);
R = [Rc{:}];
end


%% =========================================================================
%  Core analysis
%% =========================================================================
function [Rk, P] = i_analyze(X, gs, us, mk, xf, o, rs, tag)
M = size(X, 2);

% --- Transform (raw stored; transform at analysis) -----------------------
c = [];
switch xf
    case 'arcsinh'
        if ischar(o.cofactor) || isstring(o.cofactor)      % 'auto'
            c = ones(1, M);
            for j = 1:M
                v = X(X(:,j) > 0, j);
                if ~isempty(v), c(j) = prctile(v, 10); end
            end
            clab = 'auto p10';
        else
            c = o.cofactor;  if isscalar(c), c = repmat(c, 1, M); end
            clab = 'user';
        end
        Xt   = asinh(X ./ c);
        xlab = sprintf('arcsinh(x/c), c=%s', clab);
    case 'log2'
        Xt = log2(max(X, 1));  xlab = 'log2(max(x,1))';
    case 'linear'
        Xt = X;  xlab = 'linear';
        if o.prenormalized, xlab = 'linear (per-rep z)'; end
    otherwise
        error('%s unknown transform ''%s''.', tag, xf);
end

% --- Standardize per marker ----------------------------------------------
sd   = std(Xt, 0, 1);
flat = ~(sd > 0) | ~isfinite(sd);
if any(flat)
    warning('%s zero-variance markers dropped: %s', tag, strjoin(mk(flat), ', '));
    Xt(:, flat) = []; X(:, flat) = []; mk(flat) = []; sd(flat) = [];
    if ~isempty(c), c(flat) = []; end
end
assert(numel(mk) >= 2, '%s fewer than 2 markers after variance check.', tag);
Z = (Xt - mean(Xt, 1)) ./ sd;
N = size(Z, 1);

% --- PCA; project out drop_pcs (technical amplitude axis convention) ----
[coef, score, latent] = pca(Z);
ve  = latent / sum(latent);
nPC = size(coef, 2);
dp  = o.drop_pcs(o.drop_pcs >= 1 & o.drop_pcs <= nPC - 1);
if numel(dp) < numel(o.drop_pcs)
    warning('%s drop_pcs trimmed to %s (must leave >=1 PC).', tag, mat2str(dp));
end
pc1_uniform = all(coef(:,1) > 0) || all(coef(:,1) < 0);
Zd     = Z - score(:, dp) * coef(:, dp)';
keepPC = setdiff(1:nPC, dp);

% --- Groups and units ----------------------------------------------------
[Gi, glab] = findgroups(gs);
k  = numel(glab);
ng = accumarray(Gi, 1, [k 1]);

Ui = []; gOfU = []; perm_level = 'cell';
if ~isempty(us)
    Ui    = findgroups(us);
    nGinU = accumarray(Ui, Gi, [], @(x) numel(unique(x)));
    if any(nGinU > 1)
        warning('%s %d units carry >1 group: unit-level partition/permutation disabled.', tag, nnz(nGinU > 1));
        Ui = [];
    else
        gOfU = accumarray(Ui, Gi, [], @(x) x(1));
        perm_level = 'unit';
    end
end

% --- Variance partition (drop_pcs removed) --------------------------------
[SST, SSB, SSU] = i_ss(Zd, Gi, k, Ui, gOfU);
if isempty(Ui)
    comp = ["treatment"; "within"];          SS = [SSB; SST - SSB];
else
    comp = ["treatment"; "well|treatment"; "cell"];  SS = [SSB; SSU; SST - SSB - SSU];
end
partition = table(comp, SS, SS / SST, 'VariableNames', {'component','SS','frac'});
F = (SSB / (k - 1)) / ((SST - SSB) / (N - k));

% --- Permutation null on SSB (labels shuffled at unit level) -------------
if strcmp(perm_level, 'unit')
    nU  = max(Ui);
    nu  = accumarray(Ui, 1, [nU 1]);
    Su  = i_gsum(Zd, Ui, nU);
    upg = accumarray(gOfU, 1, [k 1]);
    ssb = @(gU) sum(sum(i_gsum(Su, gU, k).^2, 2) ./ accumarray(gU, nu, [k 1])) - N * sum(mean(Zd,1).^2);
    obs = ssb(gOfU);
    nul = zeros(o.nperm, 1);
    for b = 1:o.nperm, nul(b) = ssb(gOfU(randperm(rs, nU))); end
    cnt = nnz(nul >= obs - 1e-10 * abs(obs));
    logD = gammaln(nU + 1) - sum(gammaln(upg + 1));
else
    upg = nan(k, 1);
    ssb = @(G) sum(sum(i_gsum(Zd, G, k).^2, 2) ./ ng) - N * sum(mean(Zd,1).^2);
    obs = ssb(Gi);
    nul = zeros(o.nperm, 1);
    for b = 1:o.nperm, nul(b) = ssb(Gi(randperm(rs, N))); end
    cnt = nnz(nul >= obs - 1e-10 * abs(obs));
    logD = gammaln(N + 1) - sum(gammaln(ng + 1));
end
p_perm = (1 + cnt) / (o.nperm + 1);

% --- Per-marker eta^2 (on Z, before PC removal: columns stay interpretable)
mu   = mean(Z, 1);
gmZ  = i_gsum(Z, Gi, k) ./ ng;
SSTj = sum((Z - mu).^2, 1);
eta2 = (ng' * (gmZ - mu).^2) ./ SSTj;
eta2u = nan(1, numel(mk));
if ~isempty(Ui)
    umZ   = i_gsum(Z, Ui, max(Ui)) ./ accumarray(Ui, 1);
    eta2u = (accumarray(Ui, 1)' * (umZ - gmZ(gOfU, :)).^2) ./ SSTj;
end
eta2_tbl = table(string(mk(:)), eta2(:), eta2u(:), 'VariableNames', {'marker','eta2','eta2_unit'});

% --- Group geometry ------------------------------------------------------
[spread, nearest, jn, rg, Dcen] = i_geom(Zd, Gi, k, ng);
[~, ~, ~, rg_keep]        = i_geom(Z,  Gi, k, ng);    % sensitivity: no PC removal
grp = table(glab(:), ng, upg, spread, nearest, glab(jn), spread ./ nearest, ...
      'VariableNames', {'group','n','n_units','spread','nearest_dist','nearest_group','ratio'});

% --- Verdict -------------------------------------------------------------
ft = partition.frac(1);
if     ft < o.eta_lo, verdict = 'not-defining';
elseif ft > o.eta_hi, verdict = 'defining';
else,                 verdict = 'biasing';
end

% --- Pack ----------------------------------------------------------------
Rk = struct();
Rk.markers  = mk;               Rk.n = N;              Rk.groups = cellstr(glab);
Rk.transform_label = xlab;      Rk.cofactor = c;
Rk.drop_pcs = dp;               Rk.var_expl = ve;      Rk.pc1_uniform_sign = pc1_uniform;
Rk.partition = partition;       Rk.eta2 = eta2_tbl;    Rk.group = grp;
Rk.ratio_global = rg;           Rk.ratio_global_keepPC = rg_keep;
Rk.centroid_dist = Dcen;
Rk.F = F;  Rk.p_perm = p_perm;  Rk.nperm = o.nperm;    Rk.n_distinct = exp(logD);
Rk.perm_level = perm_level;     Rk.verdict = verdict;    Rk.null_frac = nul / SST;
Rk.set_tag = '';  Rk.row_idx = [];  Rk.figs = gobjects(0);

% Plot payload (not returned)
P.Xt = Xt;  P.Gi = Gi;  P.score = score(:, keepPC);  P.ve = ve(keepPC);  P.pcid = keepPC;
end


%% =========================================================================
%  Helpers: statistics
%% =========================================================================
function S = i_gsum(Z, idx, n)
% Row sums of Z per index (sparse indicator product)
S = full(sparse(idx(:), (1:numel(idx))', 1, n, numel(idx)) * Z);
end

function [SST, SSB, SSU] = i_ss(Z, Gi, k, Ui, gOfU)
mu  = mean(Z, 1);
SST = sum((Z - mu).^2, 'all');
ng  = accumarray(Gi, 1, [k 1]);
gm  = i_gsum(Z, Gi, k) ./ ng;
SSB = sum(ng .* sum((gm - mu).^2, 2));
SSU = NaN;
if ~isempty(Ui)
    nU  = max(Ui);
    nu  = accumarray(Ui, 1, [nU 1]);
    um  = i_gsum(Z, Ui, nU) ./ nu;
    SSU = sum(nu .* sum((um - gm(gOfU, :)).^2, 2));
end
end

function [spread, nearest, jn, rg, Dcen] = i_geom(Z, Gi, k, ng)
% spread  = mean cell distance to own centroid; nearest = closest other centroid
gm     = i_gsum(Z, Gi, k) ./ ng;
spread = accumarray(Gi, vecnorm(Z - gm(Gi, :), 2, 2), [k 1], @mean);
Dc     = squareform(pdist(gm));
Dcen   = Dc;                                % full centroid distance matrix
Dc(1:k+1:end) = Inf;
[nearest, jn] = min(Dc, [], 2);
rg = (sum(ng .* spread) / sum(ng)) / mean(pdist(gm));
end


%% =========================================================================
%  Helpers: inputs
%% =========================================================================
function [mk, mk_disp] = i_resolve_markers(T, o, rep_col, tag)
vn  = T.Properties.VariableNames;
sfx = '';  if ~isempty(o.compartment), sfx = ['_' o.compartment]; end
if ~isempty(o.markers)
    mk = o.markers(:)';
    if ~isempty(sfx)
        need = ~endsWith(mk, sfx);  mk(need) = strcat(mk(need), sfx);
    end
    miss = setdiff(mk, vn, 'stable');
    assert(isempty(miss), '%s marker columns not found: %s', tag, strjoin(miss, ', '));
else
    isnum = varfun(@(x) isnumeric(x) && ~islogical(x), T, 'OutputFormat', 'uniform');
    excl  = ismember(vn, [o.id_cols(:)', {o.group, o.unit, rep_col}]);
    mk    = vn(isnum & ~excl);
    if ~isempty(sfx)
        hit = endsWith(mk, sfx);
        assert(any(hit), ['%s no columns end in ''%s''. cycif_build_table strips the ', ...
               'compartment suffix (fixed at build time): pass ''compartment'','''' for its output.'], tag, sfx);
        mk = mk(hit);
    end
end
mk_disp = mk;
if ~isempty(sfx), mk_disp = regexprep(mk, [sfx '$'], ''); end
if o.vb >= 2, fprintf('%s %d candidate markers: %s\n', tag, numel(mk), strjoin(mk_disp, ', ')); end
end

function S = i_nan_sets(X, mk, g, policy, min_n, tag, vb)
% Stain-signature diagnostic + analysis-set selection
P = ~isnan(X);
[usig, ~, ic] = unique(P, 'rows');
nG = numel(unique(g));
if vb >= 2, fprintf('%s stain signatures (%d):\n', tag, size(usig, 1)); end
for s = 1:size(usig, 1) * (vb >= 2)
    lst = strjoin(mk(usig(s,:)), ', ');  if isempty(lst), lst = '(none)'; end
    fprintf('   [%d] n=%7d  groups=%d/%d  %s\n', s, nnz(ic == s), ...
            numel(unique(g(ic == s))), nG, lst);
end

S = struct('cols', {}, 'rows', {}, 'tag', {});
switch policy
    case 'complete'
        S(1) = struct('cols', true(1, size(X,2)), 'rows', all(P, 2), 'tag', 'complete');
    case 'block'
        % Candidate blocks = observed signatures; score (covers all groups, n*|S|, |S|)
        best = [-1 -1 -1];  bc = [];  br = [];
        for s = 1:size(usig, 1)
            cols = usig(s, :);  if nnz(cols) < 2, continue; end
            rows = all(P(:, cols), 2);
            sc   = [numel(unique(g(rows))) == nG, nnz(rows) * nnz(cols), nnz(cols)];
            if i_lexgt(sc, best), best = sc; bc = cols; br = rows; end
        end
        assert(~isempty(bc), '%s no signature with >=2 markers.', tag);
        S(1) = struct('cols', bc, 'rows', br, 'tag', 'block');
        if any(~bc) && vb >= 1, fprintf('%s block excludes markers: %s\n', tag, strjoin(mk(~bc), ', ')); end
        if vb >= 2, fprintf('%s block retains %d/%d cells, %d/%d markers.\n', tag, nnz(br), numel(br), nnz(bc), numel(bc)); end
    case 'panel'
        for s = 1:size(usig, 1)
            cols = usig(s, :);  rows = ic == s;
            if nnz(cols) < 2 || nnz(rows) < min_n, continue; end
            S(end+1) = struct('cols', cols, 'rows', rows, 'tag', sprintf('panel%d', s)); %#ok<AGROW>
        end
        assert(~isempty(S), '%s no panel with >=2 markers and >=%d cells.', tag, min_n);
end
end

function tf = i_lexgt(a, b)
% Lexicographic a > b
d  = find(a ~= b, 1);
tf = ~isempty(d) && a(d) > b(d);
end

function d = i_figdir(o)
% Priority: dataloc.fold.fig > dataloc.exp.fig_dir > savepath
d = '';
if ~isempty(o.dataloc)
    try, d = char(o.dataloc.fold.fig); catch, end
    if isempty(d), try, d = char(o.dataloc.exp.fig_dir); catch, end, end
end
if isempty(d), d = char(o.savepath); end
end


%% =========================================================================
%  Console interpretation
%% =========================================================================
function i_report(R, o, tag)
fprintf('\n%s ===== %s =====\n', tag, R.set_tag);
fprintf('  %d cells | %d groups (%s) | %d markers: %s\n', R.n, numel(R.groups), ...
        o.group, numel(R.markers), strjoin(R.markers, ', '));
fprintf('  scale: %s | PCs removed: %s | permutation level: %s\n', ...
        R.transform_label, mat2str(R.drop_pcs), R.perm_level);

% PC1 convention check
if ismember(1, R.drop_pcs)
    if R.pc1_uniform_sign
        fprintf('  PC1 (%.1f%% var) has uniform-sign loadings: consistent with the technical amplitude axis; removal justified.\n', 100*R.var_expl(1));
    else
        fprintf('  CAUTION: PC1 (%.1f%% var) has MIXED-sign loadings: it may carry biology in this marker set; the drop_pcs=[1] convention is not supported here. Compare the keepPC ratio below.\n', 100*R.var_expl(1));
    end
end

% Variance partition
fprintf('\n  Variance partition (standardized functional space):\n');
for i = 1:height(R.partition)
    fprintf('    %-16s %6.1f%%\n', R.partition.component(i), 100*R.partition.frac(i));
end
ft = R.partition.frac(1);
switch R.verdict
    case 'not-defining'
        fprintf('  -> Treatment explains %.1f%% (< %.0f%%): treatment does NOT define functional state. Cells under one condition occupy largely the same functional space as cells under the others; the "treatment defines state" default fails for this marker set.\n', 100*ft, 100*o.eta_lo);
    case 'defining'
        fprintf('  -> Treatment explains %.1f%% (> %.0f%%): treatment largely determines the functional configuration; within-group heterogeneity is the minority component.\n', 100*ft, 100*o.eta_hi);
    otherwise
        fprintf('  -> Treatment explains %.1f%% (between %.0f%% and %.0f%%): treatment BIASES functional state without defining it. Most variance lies within conditions.\n', 100*ft, 100*o.eta_lo, 100*o.eta_hi);
end
if height(R.partition) == 3
    fw = R.partition.frac(2);
    if fw >= 0.5 * ft
        fprintf('  -> Well|treatment = %.1f%% (%.2fx the treatment component, >= 0.5x): well-to-well variation (plating, density, per-well stain/registration) rivals treatment; the treatment effect is weakly separated from well effects at %s wells/group.\n', ...
                100*fw, fw/ft, mat2str(R.group.n_units'));
    else
        fprintf('  -> Well|treatment = %.1f%% (%.2fx treatment): well effects are small relative to treatment.\n', 100*fw, fw/ft);
    end
end
fprintf('  -> The residual cell component includes measurement noise (segmentation, registration drift, stain variability): biological within-group heterogeneity is <= this value. A reproducibility ceiling is needed to split them.\n');

% Permutation
fprintf('\n  Pseudo-F = %.2f; permutation p = %.3g (%d %s-level permutations).\n', R.F, R.p_perm, R.nperm, R.perm_level);
if R.n_distinct < R.nperm
    fprintf('  CAUTION: only %d distinct label arrangements exist: minimum achievable p ~ %.3g; p is coarse. Add wells, not cells, to sharpen it.\n', round(R.n_distinct), 1/R.n_distinct);
end
fprintf('  -> p tests whether treatment shifts the functional centroid at all. With n=%d it is expected to be small; the effect size above, not p, answers whether treatment DEFINES state.\n', R.n);

% Dispersion ratio
fprintf('\n  Dispersion ratio = %.2f (mean cell distance to own treatment centroid / mean distance between centroids).\n', R.ratio_global);
if R.ratio_global > 1
    fprintf('  -> > 1: a typical cell lies farther from its own centroid than centroids lie from each other; treatment groups overlap heavily in functional space.\n');
else
    fprintf('  -> <= 1: centroids are separated relative to within-group spread; treatments occupy distinguishable regions.\n');
end
if ~isempty(R.drop_pcs)
    fprintf('  Sensitivity (PCs retained): ratio = %.2f.', R.ratio_global_keepPC);
    if abs(log(R.ratio_global_keepPC / R.ratio_global)) > log(1.25)
        fprintf(' Differs by >25%%: the removed axis materially changes the conclusion; justify the removal biologically.\n');
    else
        fprintf(' Conclusion robust to PC removal.\n');
    end
end

% Per group
fprintf('\n  %-28s %7s %6s %8s %8s  %-22s %6s\n', 'group', 'n', 'wells', 'spread', 'nearest', 'nearest group', 'ratio');
for i = 1:height(R.group)
    flag = '';  if R.group.n(i) < o.min_n, flag = '  <min_n'; end
    fprintf('  %-28s %7d %6g %8.2f %8.2f  %-22s %6.2f%s\n', R.group.group(i), R.group.n(i), ...
            R.group.n_units(i), R.group.spread(i), R.group.nearest_dist(i), ...
            R.group.nearest_group(i), R.group.ratio(i), flag);
end
[~, imx] = max(R.group.ratio);  [~, imn] = min(R.group.ratio);
fprintf('  -> Most dispersed relative to its neighbor: %s (ratio %.2f); most distinct: %s (ratio %.2f).\n', ...
        R.group.group(imx), R.group.ratio(imx), R.group.group(imn), R.group.ratio(imn));

% Per marker
[~, ord] = sort(R.eta2.eta2, 'descend');
e = R.eta2(ord, :);
fprintf('\n  Per-marker treatment-explained variance (eta^2):\n');
for i = 1:height(e)
    fprintf('    %-14s %.3f', e.marker(i), e.eta2(i));
    if ~isnan(e.eta2_unit(i)), fprintf('   (well|tx %.3f)', e.eta2_unit(i)); end
    fprintf('\n');
end
hi = e.marker(e.eta2 > o.eta_hi);  lo = e.marker(e.eta2 < o.eta_lo);
if ~isempty(hi), fprintf('  -> Strongly treatment-determined functions: %s.\n', strjoin(hi, ', ')); end
if ~isempty(lo), fprintf('  -> Functions essentially unconstrained by treatment: %s. Their variation arises within conditions (candidate targets for dynamics-based prediction).\n', strjoin(lo, ', ')); end
fprintf('\n');
end


function i_brief(R, o, tag)
% Compact console summary: one line per question, cautions only when triggered
P  = R.partition;  ft = P.frac(1);
fw = NaN;  if height(P) == 3, fw = P.frac(2); end
fprintf('\n%s %s | %d cells | %d groups (%s) | %d markers | removed PC%s\n', tag, R.set_tag, R.n, ...
        numel(R.groups), o.group, numel(R.markers), mat2str(R.drop_pcs));
fprintf('  Variance split : treatment %.1f%% | well %.1f%% | cell %.1f%%  -> %s (perm p = %.2g, %s-level)\n', ...
        100*ft, 100*fw, 100*P.frac(end), upper(R.verdict), R.p_perm, R.perm_level);
[~, imn] = min(R.group.ratio);  [~, imx] = max(R.group.ratio);
fprintf('  Separability   : dispersion ratio %.2f (>1 = overlap); most distinct %s (%.2f), least %s (%.2f)\n', ...
        R.ratio_global, R.group.group(imn), R.group.ratio(imn), R.group.group(imx), R.group.ratio(imx));
hi = R.eta2.marker(R.eta2.eta2 >= o.eta_lo);
lo = R.eta2.marker(R.eta2.eta2 <  o.eta_lo);
fprintf('  Per marker     : treatment-set (eta2 >= %.2f) %s | unconstrained %s\n', o.eta_lo, i_lst(hi), i_lst(lo));
% Cautions (only when triggered)
if ismember(1, R.drop_pcs) && ~R.pc1_uniform_sign
    fprintf('  ! PC1 has mixed-sign loadings: removing it may remove biology.\n');
end
if ~isnan(fw) && fw >= ft
    fprintf('  ! Well-to-well variance >= treatment: effect weakly separated from plating/staining.\n');
end
if R.n_distinct < R.nperm
    fprintf('  ! Only %d distinct well shuffles: p is coarse (min %.2g).\n', round(R.n_distinct), 1/R.n_distinct);
end
end

function s = i_lst(c)
c = cellstr(c);
if isempty(c), s = '{}'; else, s = ['{' strjoin(c, ', ') '}']; end
end

%% =========================================================================
%  Step-by-step summary figure
%% =========================================================================
function figs = i_story(R, P, o, rs, fig_dir)
k   = numel(R.groups);  M = numel(R.markers);  col = i_colors(k);
f   = figure('Color', 'w', 'Name', ['Q1 treatment dispersion ' R.set_tag], 'Units', 'pixels');
f.Position(3:4) = [1350 760];
tl  = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q1  Does treatment define functional state?   (%s, %d cells, %d markers)', ...
      o.group, R.n, M), 'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
cP = [0.20 0.45 0.70; 0.90 0.60 0.20; 0.82 0.82 0.82];

% Step 1: variance partition
ax = nexttile(tl);  hold(ax, 'on');
fr = R.partition.frac(:)';  x0 = [0 cumsum(fr(1:end-1))];
for i = 1:numel(fr)
    rectangle(ax, 'Position', [x0(i) 0.75 fr(i) 0.5], 'FaceColor', cP(i,:), 'EdgeColor', 'w');
    if fr(i) > 0.04
        text(ax, x0(i) + fr(i)/2, 1, sprintf('%s\n%.1f%%', R.partition.component(i), 100*fr(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Interpreter', 'none');
    end
end
xline(ax, [o.eta_lo o.eta_hi], ':k');
xlim(ax, [0 1]);  ylim(ax, [0.5 1.5]);  set(ax, 'YTick', []);
xlabel(ax, 'fraction of functional variance');
i_step_title(ax, 'Step 1  Size', 'How is variance split between treatment, wells and cells?');
i_sty(ax);

% Step 2: permutation null
ax = nexttile(tl);  hold(ax, 'on');
histogram(R.null_frac, 30, 'Parent', ax, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
xline(ax, fr(1), '-', 'Color', cP(1,:), 'LineWidth', 2);
text(ax, fr(1), ax.YLim(2)*0.95, sprintf('  observed %.1f%%', 100*fr(1)), 'Color', cP(1,:), 'FontSize', 8);
xlabel(ax, 'treatment fraction under shuffled labels');  ylabel(ax, 'shuffles');
i_step_title(ax, 'Step 2  Reliability', sprintf('Beyond %s-level label shuffles?  p = %.2g', R.perm_level, R.p_perm));
i_sty(ax);

% Step 3: per-marker eta^2
ax = nexttile(tl);  hold(ax, 'on');
E = R.eta2;  [~, ord] = sort(E.eta2);  E = E(ord, :);
Yb = E.eta2;  lg = {'treatment'};
if ~all(isnan(E.eta2_unit)), Yb = [Yb E.eta2_unit];  lg{end+1} = 'well within treatment'; end
b = barh(ax, Yb, 'grouped', 'EdgeColor', 'none', 'BarWidth', 0.85);
for i = 1:numel(b), b(i).FaceColor = cP(i,:); end
xline(ax, o.eta_lo, ':k');
set(ax, 'YTick', 1:M, 'YTickLabel', E.marker);
legend(ax, lg, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7.5);
xlabel(ax, '\eta^2 (fraction of marker variance)');
i_step_title(ax, 'Step 3  Which functions', 'Which markers does treatment actually set?');
i_sty(ax);

% Step 4: separability per group
ax = nexttile(tl);  hold(ax, 'on');
G = R.group;  [~, ord] = sort(G.ratio);  G = G(ord, :);  y = 1:k;
for i = 1:k, plot(ax, [G.nearest_dist(i) G.spread(i)], [y(i) y(i)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 2); end
h1 = plot(ax, G.spread, y, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
h2 = plot(ax, G.nearest_dist, y, 'o', 'MarkerFaceColor', cP(1,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
for i = 1:k
    text(ax, max(G.spread(i), G.nearest_dist(i)), y(i), sprintf('  x%.1f', G.ratio(i)), 'FontSize', 7.5);
end
set(ax, 'YTick', y, 'YTickLabel', G.group);  ylim(ax, [0.5 k+0.5]);
legend(ax, [h1 h2], {'typical cell to own centroid', 'to nearest other condition'}, 'Location', 'southoutside', 'Box', 'off', 'FontSize', 7.5);
xlabel(ax, 'distance (standardized functional units)');
i_step_title(ax, 'Step 4  Separability', 'Is each condition its own region? (x = spread / gap)');
i_sty(ax);

% Functional map
ax = nexttile(tl);  hold(ax, 'on');
if size(P.score, 2) >= 2
    sub = i_subsample(P.Gi, k, min(o.max_plot_pts, 800), rs);
    h = gobjects(k, 1);
    for gi = 1:k
        idx   = sub(P.Gi(sub) == gi);
        h(gi) = scatter(ax, P.score(idx,1), P.score(idx,2), 4, col(gi,:), 'filled', 'MarkerFaceAlpha', 0.25);
    end
    for gi = 1:k
        cxy = mean(P.score(P.Gi == gi, 1:2), 1);
        plot(ax, cxy(1), cxy(2), 'o', 'MarkerFaceColor', col(gi,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 9);
    end
    legend(ax, h, R.groups, 'Interpreter', 'none', 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 7);
    xlabel(ax, sprintf('PC%d (%.0f%%)', P.pcid(1), 100*P.ve(1)));
    ylabel(ax, sprintf('PC%d (%.0f%%)', P.pcid(2), 100*P.ve(2)));
end
i_step_title(ax, 'Map', 'Cells and condition centroids (filled circles)');
i_sty(ax);

% Reading guide
ax = nexttile(tl);
ft = fr(1);
i_textpanel(ax, { ...
  'WHY THIS ORDER', ...
  '1  Size first: how much of the functional variance does', ...
  '   treatment explain, compared with wells and single cells?', ...
  '2  Then reliability: is that share larger than shuffling', ...
  '   whole wells between conditions would give? (size without', ...
  '   reliability is noise; reliability without size is trivial', ...
  '   at this n)', ...
  '3  Then which functions carry the treatment effect.', ...
  '4  Then geometry: are conditions separate regions, or does', ...
  '   a typical cell sit closer to another condition (x > 1)?', ...
  '', ...
  sprintf('VERDICT  %s', upper(R.verdict)), ...
  sprintf('treatment %.1f%% (thresholds %.0f%% / %.0f%%), p = %.2g', 100*ft, 100*o.eta_lo, 100*o.eta_hi, R.p_perm), ...
  sprintf('dispersion ratio %.2f (> 1 = conditions overlap)', R.ratio_global), ...
  '', ...
  'NEXT  markers with low treatment eta2 vary within', ...
  '      conditions: targets for Q3 (do dynamics predict them?)'});

ctag  = sprintf('%dmk', M);  if M <= 4, ctag = strjoin(R.markers, '-'); end
figs = f;  i_save(f, o, fig_dir, 'Q1-summary', ctag, sprintf('%s_%s', o.group, R.set_tag));
end

%% =========================================================================
%  Plots
%% =========================================================================
function figs = i_plots(R, P, o, rs, fig_dir)
k   = numel(R.groups);
M   = numel(R.markers);
col = i_colors(k);
sub = i_subsample(P.Gi, k, o.max_plot_pts, rs);
ctag = sprintf('%dmk', M);  if M <= 4, ctag = strjoin(R.markers, '-'); end
facet = sprintf('%s_%s', o.group, R.set_tag);
figs  = gobjects(0);

% 1 --- Per-marker strip + box by group ------------------------------------
f = figure('Color', 'w', 'Name', ['txdisp strip ' R.set_tag]);
tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact');
for j = 1:M
    ax = nexttile(tl);  hold(ax, 'on');
    for gi = 1:k
        idx = sub(P.Gi(sub) == gi);
        scatter(ax, gi + 0.5*(rand(rs, numel(idx), 1) - 0.5), P.Xt(idx, j), 4, col(gi,:), ...
                'filled', 'MarkerFaceAlpha', 0.2);
    end
    boxchart(ax, P.Gi, P.Xt(:, j), 'BoxFaceColor', [0.2 0.2 0.2], 'BoxFaceAlpha', 0, ...
             'MarkerStyle', 'none', 'LineWidth', 1);
    set(ax, 'XTick', 1:k, 'XTickLabel', R.groups, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
    title(ax, sprintf('%s  (eta2 = %.2f)', R.markers{j}, R.eta2.eta2(j)), 'Interpreter', 'none');
    ylabel(ax, R.transform_label, 'Interpreter', 'none');
end
figs(end+1) = f;  i_save(f, o, fig_dir, 'txdisp-strip', ctag, facet);

% 2 --- Dispersion geometry + partition ----------------------------------
f = figure('Color', 'w', 'Name', ['txdisp dispersion ' R.set_tag]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact');
ax = nexttile(tl);
bar(ax, [R.group.spread R.group.nearest_dist]);
set(ax, 'XTickLabel', R.groups, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
legend(ax, {'within spread', 'nearest centroid'}, 'Location', 'best');
ylabel(ax, 'distance (std. units)');  title(ax, 'Within vs between');
ax = nexttile(tl);
b = bar(ax, R.group.ratio, 'FaceColor', 'flat');  b.CData = col;
yline(ax, 1, '--k', 'overlap threshold');
set(ax, 'XTickLabel', R.groups, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
ylabel(ax, 'spread / nearest centroid');
title(ax, sprintf('Dispersion ratio (global %.2f)', R.ratio_global));
ax = nexttile(tl);
bar(ax, [R.partition.frac'; nan(1, height(R.partition))], 'stacked');   % NaN row forces one stacked bar
set(ax, 'XTick', []);  xlim(ax, [0.4 1.6]);  ylim(ax, [0 1]);  ylabel(ax, 'fraction of variance');
legend(ax, cellstr(R.partition.component), 'Location', 'eastoutside');
title(ax, sprintf('Partition (%s)', R.verdict));
figs(end+1) = f;  i_save(f, o, fig_dir, 'txdisp-dispersion', ctag, facet);

% 3 --- Per-marker eta^2 --------------------------------------------------
f  = figure('Color', 'w', 'Name', ['txdisp eta2 ' R.set_tag]);
ax = axes(f);
[~, ord] = sort(R.eta2.eta2, 'ascend');
Y = R.eta2.eta2(ord);
if ~all(isnan(R.eta2.eta2_unit))
    Y = [Y R.eta2.eta2_unit(ord)];
    barh(ax, Y, 'stacked');  legend(ax, {'treatment', 'well|treatment'}, 'Location', 'southeast');
else
    barh(ax, Y);
end
xline(ax, o.eta_lo, ':k');  xline(ax, o.eta_hi, ':k');
set(ax, 'YTick', 1:M, 'YTickLabel', R.eta2.marker(ord), 'TickLabelInterpreter', 'none');
xlabel(ax, 'fraction of marker variance');  title(ax, 'Treatment-explained variance per marker');
figs(end+1) = f;  i_save(f, o, fig_dir, 'txdisp-eta2', ctag, facet);

% 4 --- Functional PCA scatter (after PC removal) ------------------------
if size(P.score, 2) >= 2
    f  = figure('Color', 'w', 'Name', ['txdisp pca ' R.set_tag]);
    ax = axes(f);  hold(ax, 'on');
    h  = gobjects(k, 1);
    for gi = 1:k
        idx   = sub(P.Gi(sub) == gi);
        h(gi) = scatter(ax, P.score(idx,1), P.score(idx,2), 5, col(gi,:), 'filled', 'MarkerFaceAlpha', 0.25);
        cxy   = mean(P.score(P.Gi == gi, 1:2), 1);
        plot(ax, cxy(1), cxy(2), 'x', 'Color', col(gi,:)*0.6, 'MarkerSize', 14, 'LineWidth', 3);
    end
    legend(ax, h, R.groups, 'Interpreter', 'none', 'Location', 'bestoutside');
    xlabel(ax, sprintf('PC%d (%.1f%%)', P.pcid(1), 100*P.ve(1)));
    ylabel(ax, sprintf('PC%d (%.1f%%)', P.pcid(2), 100*P.ve(2)));
    title(ax, sprintf('Functional PCA, removed PCs %s (x = centroid)', mat2str(R.drop_pcs)));
    figs(end+1) = f;  i_save(f, o, fig_dir, 'txdisp-pca', ctag, facet);
end
end

function sub = i_subsample(Gi, k, nmax, rs)
% Per-group subsample indices for point plots
sub = [];
for gi = 1:k
    idx = find(Gi == gi);
    if numel(idx) > nmax, idx = idx(randperm(rs, numel(idx), nmax)); end
    sub = [sub; idx]; %#ok<AGROW>
end
end

function col = i_colors(k)
% Lab 'bright' palette (as cycif_* plotting), lines() beyond 9
base = [0.900 0.100 0.100; 0.100 0.500 0.900; 0.100 0.850 0.300;
        0.000 0.620 0.200; 0.000 0.390 0.100; 0.900 0.100 0.550;
        1.000 0.700 0.000; 0.960 0.480 0.000; 0.980 0.860 0.100];
if k <= size(base, 1), col = base(1:k, :); else, col = lines(k); end
end

function i_save(f, o, d, ptype, ctag, facet)
% <yyyy-mm-dd>__<plottype>_<channels>_<facet>.svg ; overwrite on re-run
if ~o.save || isempty(d), return; end
if ~exist(d, 'dir'), mkdir(d); end
fn = fullfile(d, sprintf('%s__%s_%s_%s.svg', char(datetime('today', 'Format', 'yyyy-MM-dd')), ...
              ptype, ctag, matlab.lang.makeValidName(facet)));
print(f, fn, '-dsvg', '-painters');
fprintf('  saved %s\n', fn);
end

function i_sty(ax)
% Shared axis style
set(ax, 'FontName', 'Helvetica', 'FontSize', 8.5, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 0.75, ...
        'Color', 'none', 'XColor', [0.25 0.25 0.25], 'YColor', [0.25 0.25 0.25], 'TickLabelInterpreter', 'none');
end

function i_step_title(ax, step, question)
% Two-line panel title: bold step label + plain question
title(ax, {['\bf' strrep(step, '_', '\_')], ['\rm' strrep(question, '_', '\_')]}, 'Interpreter', 'tex', 'FontSize', 9);
end

function i_textpanel(ax, lines)
% Reading guide / verdict panel
axis(ax, 'off');
text(ax, 0, 1, lines, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8.5, ...
     'Interpreter', 'none', 'FontName', 'Helvetica');
end