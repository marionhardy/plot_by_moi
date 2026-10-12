function R = dl_functional_dimensionality(T, varargin)
% dl_functional_dimensionality  Is functional-marker space low-dimensional
%   (functions bundle -> M1/M2-type axes) and multimodal (discrete states)?
%   PCA scree + participation ratio + parallel analysis, loading structure
%   vs an M1/M2 contrast, Hartigan dip tests (pooled AND within treatment),
%   GMM BIC/ICL, and a pulsatile-marker (phase-aliasing) audit.
%
% USAGE
%   R = dl_functional_dimensionality(Tf, 'prenormalized',true);
%   R = dl_functional_dimensionality(T,  'transform','arcsinh', 'dataloc',dataloc);
%   R = dl_functional_dimensionality(Tf, 'prenormalized',true, 'pulsatile_weight',0.5, 'save',false);
%   R = dl_functional_dimensionality(Tf, 'prenormalized',true, 'm2_markers',{'CD206','CD163'}, 'drop_pcs',[]);
%   R = dl_functional_dimensionality(Tf, 'prenormalized',true, 'verbose','full', 'plots','all');
%
% INPUT
%   T   table, one row per cell (cycif_build_table output or pooled Tf)
%
% NAME-VALUE
%  Columns / scale / filtering (identical semantics to dl_treatment_dispersion)
%   'group' ('tx1_label') 'unit' ('xy', parity only) 'markers' ({}) 'id_cols' 'compartment' ('') 'rep_col' ('')
%   'prenormalized' (false) 'transform' ('') 'cofactor' ('auto') 'drop_pcs' ([1])
%   'exclude_spike' (true) 'nan_policy' ('block') 'min_n' (30)
%  Biology
%   'pulsatile'        cell  markers whose snapshot level is time-unstable (opt-in, justify per dataset),
%                            matched case/punctuation-insensitive                                ({})
%   'pulsatile_weight' scalar multiplies standardized pulsatile columns before PCA (1)
%   'm1_markers'       cell  M1-side of the hypothesis contrast
%                            ({'HLA_DR','HLADR','CD86','iNOS','NOS2','IL1b','IL1beta','TNFa','TNF','IL6','IL12'})
%   'm2_markers'       cell  M2-side ({'CD206','MRC1','CD163','ARG1','IL10'})
%  Thresholds (interpretive conventions)
%   'lowdim_frac'   scalar  PR / nKeptPCs at or below which space is "low-dim" (0.5)
%   'var_target'    scalar  cumulative variance for the nPC count         (0.80)
%   'align_thresh'  scalar  |cos| of a PC with the M1-M2 contrast          (0.7)
%   'alpha'         scalar  family-wise alpha (Bonferroni within family)   (0.05)
%  Tests
%   'n_axes'        scalar  kept PCs tested for multimodality             (5)
%   'dip_max_n'     scalar  subsample per dip test                        (5000)
%   'dip_nboot'     scalar  uniform-null bootstrap draws per sample size  (500)
%   'within_group'  logical dip tests within each treatment group         (true)
%   'gmm_kmax'      scalar  max GMM components                            (6)
%   'gmm_dims'      scalar  kept PCs used by the GMM                      (3)
%   'gmm_max_n'     scalar  GMM subsample                                 (20000)
%   'pa_nperm'      scalar  parallel-analysis permutations                (50)
%   'pa_max_n'      scalar  parallel-analysis subsample                   (20000)
%   'seed'          scalar  (0)
%  Output
%   'save' (true) 'closefigs' (true) 'dataloc' ([]) 'savepath' ('') 'max_plot_pts' (5000)
%   'verbose'  'summary' = compact table + verdict | 'full' = every value interpreted | 'none' ('summary')
%   'plots'         'all' | 'none' (use dl_state_summary for compact figures) ('all')
%
% OUTPUT
%   R  struct array (numel > 1 only for 'panel'), fields:
%     .set_tag .row_idx .markers .n .transform_label .pulsatile .pulsatile_weight
%     .eigenvalues .var_expl .drop_pcs .kept_pcs .pc1_uniform_sign
%     .pr .pr_frac .pr_all .n_var_target .n_pa .pa_null95
%     .loadings (M x nPC) .m1m2_cos (per kept PC) .m1m2_best_pc .m1m2_best_cos .m1m2_contrast
%     .uniqueness  table  marker | uniqueness (1 - R^2 from all other markers) | pulsatile
%     .dip_axes    table  axis | dip | p | p_bonf | floor_frac | top_marker | top_is_pulsatile
%     .dip_markers table  marker | dip | p | p_bonf | floor_frac | pulsatile
%     .dip_within  table  group | axis | n | dip | p | p_bonf
%     .gmm         table  K | BIC | ICL ; .k_bic .k_icl
%     .nopuls      struct PR / pr_frac / n_var_target without pulsatile markers
%     .multimodal_pooled .multimodal_within .modes_source
%     .verdict     'low-dim-wells' | 'low-dim-continuum' | 'high-dim-dissociated'
%     .figs
%
% NOTES
%   New in DatalocHandler v5 (no v4 counterpart). Hartigan dip: native port
%   of AS 217 (after Mechler); p-values from a bootstrap uniform null built
%   with the same implementation, so any implementation bias cancels.

%% ---- Parse inputs -------------------------------------------------------
p = inputParser; p.FunctionName = mfilename;
isTxt = @(x) ischar(x) || isstring(x);
addRequired (p, 'T', @istable);
addParameter(p, 'group',         'tx1_label', isTxt);
addParameter(p, 'unit',          'xy',        isTxt);   % accepted for interface parity; excluded from markers
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
addParameter(p, 'pulsatile',     {});   % opt-in; no default flagging
addParameter(p, 'pulsatile_weight', 1);
addParameter(p, 'm1_markers',    {'HLA_DR','HLADR','CD86','iNOS','NOS2','IL1b','IL1beta','TNFa','TNF','IL6','IL12'});
addParameter(p, 'm2_markers',    {'CD206','MRC1','CD163','ARG1','IL10'});
addParameter(p, 'lowdim_frac',   0.5);
addParameter(p, 'var_target',    0.80);
addParameter(p, 'align_thresh',  0.7);
addParameter(p, 'alpha',         0.05);
addParameter(p, 'n_axes',        5);
addParameter(p, 'dip_max_n',     5000);
addParameter(p, 'dip_nboot',     500);
addParameter(p, 'within_group',  true);
addParameter(p, 'gmm_kmax',      6);
addParameter(p, 'gmm_dims',      3);
addParameter(p, 'gmm_max_n',     20000);
addParameter(p, 'pa_nperm',      50);
addParameter(p, 'pa_max_n',      20000);
addParameter(p, 'seed',          0);
addParameter(p, 'save',          true);
addParameter(p, 'closefigs',     true);
addParameter(p, 'dataloc',       []);
addParameter(p, 'savepath',      '',          isTxt);
addParameter(p, 'max_plot_pts',  5000);
addParameter(p, 'plots',         'summary',   @(x) any(strcmpi(x, {'summary','all','none'})));
addParameter(p, 'verbose',       'summary',   @(x) any(strcmpi(x, {'summary','full','none'})));
parse(p, T, varargin{:});
o = p.Results;
o.group       = char(o.group);       o.unit        = char(o.unit);
o.compartment = char(o.compartment);
o.rep_col     = char(o.rep_col);     o.transform   = lower(char(o.transform));
o.nan_policy  = lower(char(o.nan_policy));
o.markers     = cellstr(o.markers);  o.pulsatile   = cellstr(o.pulsatile);
o.m1_markers  = cellstr(o.m1_markers); o.m2_markers = cellstr(o.m2_markers);
o.prenormalized = logical(o.prenormalized);
o.plots = lower(char(o.plots));
o.vb = find(strcmpi(o.verbose, {'none','summary','full'})) - 1;   % 0 none, 1 summary, 2 full
rs  = RandStream('mt19937ar', 'Seed', o.seed);
tag = '[dl_functional_dimensionality]';
vn  = T.Properties.VariableNames;

%% ---- Row filtering ------------------------------------------------------
keep = true(height(T), 1);
if o.exclude_spike && ismember('ptx_is_spike', vn)
    ns = nnz(T.ptx_is_spike);
    if ns > 0 && o.vb >= 2, fprintf('%s excluding %d spike-treatment cells (ptx_is_spike).\n', tag, ns); end
    keep = keep & ~logical(T.ptx_is_spike);
end
if ismember(o.group, vn)
    g   = string(T.(o.group));
    bad = ismissing(g) | strlength(strtrim(g)) == 0;
    keep = keep & ~bad;
else
    warning('%s group column ''%s'' absent: block coverage and within-group dip tests disabled.', tag, o.group);
    g = repmat("all", height(T), 1);
    o.within_group = false;
end

rep_col = o.rep_col;
if isempty(rep_col)
    cand = {'rep','ExpID','exp'};  hit = cand(ismember(cand, vn));
    if ~isempty(hit), rep_col = hit{1}; end
end
if ~isempty(rep_col) && ~o.prenormalized && numel(unique(string(T.(rep_col)(keep)))) > 1
    warning('%s replicates pooled without prenormalization: batch offsets add a spurious shared axis and can create replicate "modes". Prefer ''prenormalized'',true.', tag);
end

%% ---- Markers + scale ----------------------------------------------------
[mk, mk_disp] = i_resolve_markers(T, o, rep_col, tag);
if o.prenormalized
    if ~isempty(o.transform) && ~strcmp(o.transform, 'linear')
        warning('%s prenormalized=true: forcing ''linear''.', tag);
    end
    xf = 'linear';
elseif isempty(o.transform), xf = 'arcsinh';
else,                        xf = o.transform;
end

X       = double(T{keep, mk});
row_idx = find(keep);
gk      = g(keep);
allnan  = all(isnan(X), 1);
if any(allnan)
    if o.vb >= 1, fprintf('%s dropping all-NaN markers: %s\n', tag, strjoin(mk_disp(allnan), ', ')); end
    X(:, allnan) = [];  mk_disp(allnan) = [];
end
assert(size(X,2) >= 3, '%s dimensionality needs >= 3 markers.', tag);

S = i_nan_sets(X, mk_disp, gk, o.nan_policy, o.min_n, tag, o.vb);

%% ---- Analyze each set ---------------------------------------------------
fig_dir = i_figdir(o);
if o.save && isempty(fig_dir)
    warning('%s save=true but no dataloc fig dir or savepath: figures not saved.', tag);
end
if o.vb >= 2
    fprintf('%s thresholds lowdim_frac=%.2f, align_thresh=%.2f are interpretive conventions.\n', ...
            tag, o.lowdim_frac, o.align_thresh);
end

dipcache = containers.Map('KeyType', 'double', 'ValueType', 'any');   % uniform nulls by n
Rc = {};
for s = 1:numel(S)
    if nnz(S(s).cols) < 3
        if o.vb >= 1, fprintf('%s [%s] <3 markers; skipped.\n', tag, S(s).tag); end;  continue
    end
    [Rk, P] = i_analyze(X(S(s).rows, S(s).cols), gk(S(s).rows), mk_disp(S(s).cols), ...
                        xf, o, rs, dipcache, tag);
    Rk.set_tag = S(s).tag;
    Rk.row_idx = row_idx(S(s).rows);
    if o.vb == 2, i_report(Rk, o, tag); elseif o.vb == 1, i_brief(Rk, o, tag); end
    if ~strcmp(o.plots, 'none')
        Rk.figs = i_story(Rk, o, fig_dir);
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
function [Rk, P] = i_analyze(X, gs, mk, xf, o, rs, dipcache, tag)

% --- Scale ---------------------------------------------------------------
[Z, ~, c, xlab, mk] = i_scale(X, mk, xf, o, tag);
N = size(Z, 1);  M = numel(mk);
nm   = i_norm(mk);
puls = ismember(nm, i_norm(o.pulsatile));
if o.pulsatile_weight ~= 1 && any(puls)
    Z(:, puls) = Z(:, puls) * o.pulsatile_weight;
end

% --- PCA -----------------------------------------------------------------
[coef, score, latent] = pca(Z);
ve  = latent / sum(latent);
nPC = size(coef, 2);
dp  = o.drop_pcs(o.drop_pcs >= 1 & o.drop_pcs <= nPC - 2);   % keep >= 2 PCs
if numel(dp) < numel(o.drop_pcs)
    warning('%s drop_pcs trimmed to %s (must leave >=2 PCs).', tag, mat2str(dp));
end
kept = setdiff(1:nPC, dp);
nk   = numel(kept);
pc1_uniform = all(coef(:,1) > 0) || all(coef(:,1) < 0);

% --- Uniqueness: 1 - R^2 of each marker on all others, AFTER removing the
%     drop_pcs axes. Axes are re-estimated without marker j (leave-one-out),
%     otherwise the projection makes every marker exactly predictable.
uniq = nan(M, 1);
for j = 1:M
    oj = setdiff(1:M, j);
    Zr = Z;
    if ~isempty(dp)
        [~, sj] = pca(Z(:, oj));
        Aj = sj(:, dp(dp <= size(sj, 2)));
        Zr = Z - Aj * (Aj \ Z);                     % residualize all markers on technical axes
    end
    bj = lsqminnorm(Zr(:, oj), Zr(:, j));
    uniq(j) = sum((Zr(:, j) - Zr(:, oj) * bj).^2) / sum(Zr(:, j).^2);
end
uniq = min(max(uniq, 0), 1);

lk     = latent(kept);
pr     = sum(lk)^2 / sum(lk.^2);          % participation ratio (effective dims)
pr_all = sum(latent)^2 / sum(latent.^2);
cumk   = cumsum(lk) / sum(lk);
n_vt   = find(cumk >= o.var_target, 1);

% --- Parallel analysis on the residual space (drop_pcs removed) ----------
Zd  = Z - score(:, dp) * coef(:, dp)';
idx = 1:N;  if N > o.pa_max_n, idx = randperm(rs, N, o.pa_max_n); end
Zs  = Zd(idx, :);  ns = numel(idx);
obs = sort(eig(cov(Zs)), 'descend');
E   = zeros(o.pa_nperm, M);
for b = 1:o.pa_nperm
    Zp = Zs;
    for j = 1:M, Zp(:, j) = Zs(randperm(rs, ns), j); end   % break covariance, keep marginals
    E(b, :) = sort(eig(cov(Zp)), 'descend')';
end
null95 = prctile(E, 95, 1);
above  = obs(1:nk)' > null95(1:nk);
n_pa   = find(~above, 1) - 1;  if isempty(n_pa), n_pa = nk; end

% --- M1/M2 contrast alignment --------------------------------------------
v  = zeros(M, 1);
i1 = ismember(nm, i_norm(o.m1_markers));  i2 = ismember(nm, i_norm(o.m2_markers));
m1m2_cos = nan(1, nk);  best_pc = NaN;  best_cos = NaN;
if any(i1) && any(i2)
    v(i1) = 1 / nnz(i1);  v(i2) = -1 / nnz(i2);
    m1m2_cos = (coef(:, kept)' * v)' / norm(v);           % loadings are unit-norm
    [best_cos, ib] = max(abs(m1m2_cos));
    best_pc  = kept(ib);
    best_cos = m1m2_cos(ib);
end

% --- Dip tests -----------------------------------------------------------
na  = min(o.n_axes, nk);
axs = kept(1:na);
[~, topm] = max(abs(coef(:, axs)), [], 1);
d_ax = zeros(na,1); p_ax = zeros(na,1); f_ax = zeros(na,1);
for a = 1:na
    [d_ax(a), p_ax(a), f_ax(a)] = i_dip_test(score(:, axs(a)), o, rs, dipcache, tag);
end
dip_axes = table(compose("PC%d", axs(:)), d_ax, p_ax, min(p_ax*na, 1), f_ax, ...
                 reshape(string(mk(topm)), [], 1), reshape(puls(topm), [], 1), 'VariableNames', ...
                 {'axis','dip','p','p_bonf','floor_frac','top_marker','top_is_pulsatile'});

d_mk = zeros(M,1); p_mk = zeros(M,1); f_mk = zeros(M,1);
for j = 1:M
    [d_mk(j), p_mk(j), f_mk(j)] = i_dip_test(Z(:, j), o, rs, dipcache, tag);
end
dip_markers = table(string(mk(:)), d_mk, p_mk, min(p_mk*M, 1), f_mk, puls(:), ...
                    'VariableNames', {'marker','dip','p','p_bonf','floor_frac','pulsatile'});

[Gi, glab] = findgroups(gs);
dip_within = table(strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
                   'VariableNames', {'group','axis','n','dip','p','p_bonf'});
if o.within_group && numel(glab) > 1
    for gi = 1:numel(glab)
        r = Gi == gi;
        if nnz(r) < o.min_n, continue; end
        for a = 1:na
            [dd, pp] = i_dip_test(score(r, axs(a)), o, rs, dipcache, tag);
            dip_within(end+1, :) = {glab(gi), "PC" + axs(a), nnz(r), dd, pp, NaN}; %#ok<AGROW>
        end
    end
    dip_within.p_bonf = min(dip_within.p * height(dip_within), 1);
end

% --- GMM model selection on leading kept PCs ------------------------------
gd  = kept(1:min(o.gmm_dims, nk));
Yg  = score(:, gd);
if N > o.gmm_max_n, Yg = Yg(randperm(rs, N, o.gmm_max_n), :); end
K   = (1:o.gmm_kmax)';
BIC = nan(size(K));  ICL = nan(size(K));
gstate = rng;  rng(o.seed);                       % fitgmdist uses the global stream
for kk = K'
    try
        gm = fitgmdist(Yg, kk, 'Replicates', 3, 'RegularizationValue', 1e-6, ...
                       'Options', statset('MaxIter', 500));
        pz = posterior(gm, Yg);
        BIC(kk) = gm.BIC;
        ICL(kk) = gm.BIC - 2 * sum(pz .* log(max(pz, realmin)), 'all');   % BIC + 2*entropy
    catch
    end
end
rng(gstate);
[~, k_bic] = min(BIC);  [~, k_icl] = min(ICL);

% --- Sensitivity: without pulsatile markers -------------------------------
nopuls = struct('markers', {mk(~puls)}, 'pr', NaN, 'pr_frac', NaN, 'n_var_target', NaN);
if any(puls) && nnz(~puls) >= 3
    [~, ~, l2] = pca(Z(:, ~puls));
    d2 = dp(dp <= numel(l2) - 2);
    l2 = l2(setdiff(1:numel(l2), d2));
    nopuls.pr           = sum(l2)^2 / sum(l2.^2);
    nopuls.pr_frac      = nopuls.pr / numel(l2);
    nopuls.n_var_target = find(cumsum(l2)/sum(l2) >= o.var_target, 1);
end

% --- Verdict -------------------------------------------------------------
lowdim   = pr / nk <= o.lowdim_frac;
mm_pool  = any(dip_axes.p_bonf < o.alpha);
mm_with  = height(dip_within) > 0 && any(dip_within.p_bonf < o.alpha);
if mm_with,     src = 'within-treatment';
elseif mm_pool, src = 'between-treatment-only';
else,           src = 'none';
end
if ~lowdim,      verdict = 'high-dim-dissociated';
elseif mm_with,  verdict = 'low-dim-wells';
else,            verdict = 'low-dim-continuum';
end

% --- Pack ----------------------------------------------------------------
Rk = struct();
Rk.markers = mk;  Rk.n = N;  Rk.transform_label = xlab;  Rk.cofactor = c;
Rk.pulsatile = mk(puls);  Rk.pulsatile_weight = o.pulsatile_weight;
Rk.eigenvalues = latent;  Rk.var_expl = ve;  Rk.drop_pcs = dp;  Rk.kept_pcs = kept;
Rk.pc1_uniform_sign = pc1_uniform;
Rk.pr = pr;  Rk.pr_frac = pr / nk;  Rk.pr_all = pr_all;  Rk.n_var_target = n_vt;
Rk.n_pa = n_pa;  Rk.pa_null95 = null95;
Rk.loadings = coef;
Rk.m1m2_cos = m1m2_cos;  Rk.m1m2_best_pc = best_pc;  Rk.m1m2_best_cos = best_cos;
Rk.m1m2_markers = {mk(i1), mk(i2)};  Rk.m1m2_contrast = v;
Rk.uniqueness = table(string(mk(:)), uniq, puls(:), 'VariableNames', {'marker','uniqueness','pulsatile'});
Rk.dip_axes = dip_axes;  Rk.dip_markers = dip_markers;  Rk.dip_within = dip_within;
Rk.gmm = table(K, BIC, ICL);  Rk.k_bic = k_bic;  Rk.k_icl = k_icl;
Rk.nopuls = nopuls;
Rk.multimodal_pooled = mm_pool;  Rk.multimodal_within = mm_with;  Rk.modes_source = src;
Rk.verdict = verdict;
Rk.set_tag = '';  Rk.row_idx = [];  Rk.figs = gobjects(0);

P.score = score;  P.Gi = Gi;  P.glab = glab;  P.puls = puls;  P.v = v;
end


%% =========================================================================
%  Dip test (Hartigan & Hartigan 1985; AS 217 port after Mechler)
%% =========================================================================
function [d, pval, floor_frac] = i_dip_test(x, o, rs, cache, tag)
% Dip statistic + bootstrap p against the uniform (least-favorable unimodal) null
x = x(~isnan(x));
floor_frac = mean(x == min(x));               % point mass at floor (clipping artifact)
if numel(x) > o.dip_max_n, x = x(randperm(rs, numel(x), o.dip_max_n)); end
if numel(unique(x)) < numel(x)                 % break ties (AS 217 divides by spacings)
    x = x + 1e-9 * (max(x) - min(x)) * rand(rs, size(x));
end
n = numel(x);
if n < 4, d = NaN; pval = NaN; return; end
d = i_dip(x);
if ~isKey(cache, n)
    if o.vb >= 2, fprintf('%s building uniform dip null (n=%d, %d draws)...\n', tag, n, o.dip_nboot); end
    nd = zeros(o.dip_nboot, 1);
    for b = 1:o.dip_nboot, nd(b) = i_dip(rand(rs, n, 1)); end
    cache(n) = sort(nd); %#ok<NASGU>  (containers.Map is a handle)
end
nd   = cache(n);
pval = (1 + nnz(nd >= d)) / (numel(nd) + 1);
end

function dip = i_dip(x)
% Dip statistic of a sample (0.5 x max distance ECDF <-> best unimodal CDF)
x = sort(x(:));
N = numel(x);
if N < 4 || x(N) == x(1), dip = 0; return; end
% Perfectly unimodal input (AS 217 can cycle): <=1 sign change in spacing trend
xs   = -sign(diff(diff(x)));
posi = find(xs > 0);  negi = find(xs < 0);
if isempty(posi) || isempty(negi) || all(posi < min(negi)), dip = 0; return; end

fn  = N;
mn  = zeros(N,1); mj = zeros(N,1); gcm = zeros(N,1); lcm = zeros(N,1);
low = 1;  high = N;  dip = 1/fn;

% Indices for the greatest convex minorant
mn(1) = 1;
for j = 2:N
    mn(j) = j - 1;
    while true
        mnj = mn(j);  mnmnj = mn(mnj);
        if mnj == 1 || (x(j)-x(mnj))*(mnj-mnmnj) < (x(mnj)-x(mnmnj))*(j-mnj), break; end
        mn(j) = mnmnj;
    end
end
% Indices for the least concave majorant
mj(N) = N;
for k = N-1:-1:1
    mj(k) = k + 1;
    while true
        mjk = mj(k);  mjmjk = mj(mjk);
        if mjk == N || (x(k)-x(mjk))*(mjk-mjmjk) < (x(mjk)-x(mjmjk))*(k-mjk), break; end
        mj(k) = mjmjk;
    end
end

itmax = 10 * N;  it = 0;
while true
    it = it + 1;
    if it > itmax, warning('i_dip: iteration cap reached (n=%d).', N); break; end
    % GCM change points high -> low
    gcm(1) = high;  ic = 2;  gcm(2) = mn(gcm(1));
    while gcm(ic) > low, ic = ic + 1; gcm(ic) = mn(gcm(ic-1)); end
    icx = ic;
    % LCM change points low -> high
    lcm(1) = low;  ic = 2;  lcm(2) = mj(lcm(1));
    while lcm(ic) < high, ic = ic + 1; lcm(ic) = mj(lcm(ic-1)); end
    icv = ic;

    ig = icx;  ih = icv;  ix = icx - 1;  iv = 2;  d = 0;
    if icx == 2 && icv == 2
        d = 1/fn;
    else
        while true
            igcmx = gcm(ix);  lcmiv = lcm(iv);
            if igcmx <= lcmiv
                lcmiv1 = lcm(iv-1);
                a  = lcmiv - lcmiv1;  b = igcmx - lcmiv1 - 1;
                dx = (x(igcmx)-x(lcmiv1))*a / (fn*(x(lcmiv)-x(lcmiv1))) - b/fn;
                ix = ix - 1;
                if dx >= d, d = dx; ig = ix + 1; ih = iv; end
            else
                igcm = gcm(ix);  igcm1 = gcm(ix+1);
                a  = lcmiv - igcm1 + 1;  b = igcm - igcm1;
                dx = a/fn - ((x(lcmiv)-x(igcm1))*b) / (fn*(x(igcm)-x(igcm1)));
                iv = iv + 1;
                if dx >= d, d = dx; ig = ix + 1; ih = iv - 1; end
            end
            if ix < 1,   ix = 1;   end
            if iv > icv, iv = icv; end
            if gcm(ix) == lcm(iv), break; end
        end
    end
    if d < dip, break; end

    % Dip on the convex-minorant side
    dl = 0;
    if ig ~= icx
        for j = ig:icx-1
            temp = 1/fn;  jb = gcm(j+1);  je = gcm(j);
            if je - jb > 1 && x(je) ~= x(jb)
                cst = (je - jb) / (fn*(x(je) - x(jb)));
                for jr = jb:je
                    t = (jr - jb + 1)/fn - (x(jr) - x(jb))*cst;
                    if t > temp, temp = t; end
                end
            end
            if dl < temp, dl = temp; end
        end
    end
    % Dip on the concave-majorant side
    du = 0;
    if ih ~= icv
        for k = ih:icv-1
            temp = 1/fn;  kb = lcm(k);  ke = lcm(k+1);
            if ke - kb > 1 && x(ke) ~= x(kb)
                cst = (ke - kb) / (fn*(x(ke) - x(kb)));
                for kr = kb:ke
                    t = (x(kr) - x(kb))*cst - (kr - kb - 1)/fn;
                    if t > temp, temp = t; end
                end
            end
            if du < temp, du = temp; end
        end
    end
    dip  = max(dip, max(dl, du));
    low  = gcm(ig);  high = lcm(ih);
end
dip = 0.5 * dip;
end


%% =========================================================================
%  Helpers: scaling / inputs (shared semantics with dl_treatment_dispersion)
%% =========================================================================
function [Z, Xt, c, xlab, mk] = i_scale(X, mk, xf, o, tag)
M = size(X, 2);  c = [];
switch xf
    case 'arcsinh'
        if ischar(o.cofactor) || isstring(o.cofactor)
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
        Xt = asinh(X ./ c);  xlab = sprintf('arcsinh(x/c), c=%s', clab);
    case 'log2'
        Xt = log2(max(X, 1));  xlab = 'log2(max(x,1))';
    case 'linear'
        Xt = X;  xlab = 'linear';
        if o.prenormalized, xlab = 'linear (per-rep z)'; end
    otherwise
        error('%s unknown transform ''%s''.', tag, xf);
end
sd   = std(Xt, 0, 1);
flat = ~(sd > 0) | ~isfinite(sd);
if any(flat)
    warning('%s zero-variance markers dropped: %s', tag, strjoin(mk(flat), ', '));
    Xt(:, flat) = [];  mk(flat) = [];  sd(flat) = [];
    if ~isempty(c), c(flat) = []; end
end
Z = (Xt - mean(Xt, 1)) ./ sd;
end

function s = i_norm(c)
% Case/punctuation-insensitive marker key
s = lower(regexprep(cellstr(c), '[^a-zA-Z0-9]', ''));
end

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
    fprintf('   [%d] n=%7d  groups=%d/%d  %s\n', s, nnz(ic == s), numel(unique(g(ic == s))), nG, lst);
end
S = struct('cols', {}, 'rows', {}, 'tag', {});
switch policy
    case 'complete'
        S(1) = struct('cols', true(1, size(X,2)), 'rows', all(P, 2), 'tag', 'complete');
    case 'block'
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
d  = find(a ~= b, 1);
tf = ~isempty(d) && a(d) > b(d);
end

function d = i_figdir(o)
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
nk = numel(R.kept_pcs);
fprintf('\n%s ===== %s =====\n', tag, R.set_tag);
fprintf('  %d cells | %d markers: %s\n', R.n, numel(R.markers), strjoin(R.markers, ', '));
fprintf('  scale: %s | PCs removed: %s | pulsatile flagged: %s (weight %.2f)\n', R.transform_label, ...
        mat2str(R.drop_pcs), i_list(R.pulsatile), R.pulsatile_weight);

if ismember(1, R.drop_pcs)
    if R.pc1_uniform_sign
        fprintf('  PC1 (%.1f%% var) uniform-sign: technical amplitude axis; dimensionality is computed on the remaining %d PCs.\n', 100*R.var_expl(1), nk);
    else
        fprintf('  CAUTION: PC1 (%.1f%% var) is MIXED-sign: may be biological here. With it retained, PR = %.2f.\n', 100*R.var_expl(1), R.pr_all);
    end
end

% Dimensionality
fprintf('\n  Effective dimensionality:\n');
fprintf('    participation ratio PR = %.2f of %d kept PCs (PR/n = %.2f; 1 = all functions bundle on one axis, 1.0*n = fully independent)\n', R.pr, nk, R.pr_frac);
fprintf('    %d PCs reach %.0f%% of kept variance; parallel analysis retains %d PCs above the permutation null.\n', R.n_var_target, 100*o.var_target, R.n_pa);
if R.pr_frac <= o.lowdim_frac
    fprintf('  -> PR/n <= %.2f: functions BUNDLE onto few axes; they co-vary, so axis labels (M1/M2-type) can survive as a compression.\n', o.lowdim_frac);
else
    fprintf('  -> PR/n > %.2f: functions DISSOCIATE across many axes; labels that assume a single activation axis impose structure the data lacks.\n', o.lowdim_frac);
end
fprintf('  -> Parallel analysis at n=%d is liberal (null eigenvalues are tight at large n); read PR as the effect size, PA as the floor.\n', R.n);

% M1/M2 hypothesis
if isnan(R.m1m2_best_cos)
    fprintf('\n  M1/M2 contrast: not testable (needs >=1 present marker on each side; present M1: %s | M2: %s).\n', ...
            i_list(R.m1m2_markers{1}), i_list(R.m1m2_markers{2}));
else
    ib = find(R.kept_pcs == R.m1m2_best_pc);
    vk = R.var_expl(R.kept_pcs) / sum(R.var_expl(R.kept_pcs));
    fprintf('\n  M1/M2 contrast (M1: %s | M2: %s): best alignment PC%d, cos = %+.2f, carrying %.1f%% of kept variance.\n', ...
            i_list(R.m1m2_markers{1}), i_list(R.m1m2_markers{2}), R.m1m2_best_pc, R.m1m2_best_cos, 100*vk(ib));
    if abs(R.m1m2_best_cos) >= o.align_thresh
        fprintf('  -> |cos| >= %.2f: the M1-M2 contrast IS a principal axis; the labels compress a real co-variation pattern.\n', o.align_thresh);
        if ib > 1, fprintf('     It is not the leading kept axis: other functional programs dominate the variance.\n'); end
    else
        fprintf('  -> |cos| < %.2f: no principal axis matches the M1-M2 contrast; polarization is not how these functions co-vary.\n', o.align_thresh);
    end
end

% Uniqueness
[~, ord] = sort(R.uniqueness.uniqueness, 'descend');
fprintf('\n  Marker uniqueness (1 - R^2 from all other markers, after removing PCs %s):\n', mat2str(R.drop_pcs));
for i = ord'
    fl = '';  if R.uniqueness.pulsatile(i), fl = '  [pulsatile]'; end
    fprintf('    %-14s %.2f%s\n', R.uniqueness.marker(i), R.uniqueness.uniqueness(i), fl);
end
fprintf('  -> High uniqueness = the function varies independently of every other measured function (or is noise-dominated).\n');

% Multimodality: pooled axes
fprintf('\n  Hartigan dip, pooled cells (Bonferroni over %d axes):\n', height(R.dip_axes));
for i = 1:height(R.dip_axes)
    fl = '';
    if R.dip_axes.floor_frac(i) > 0.01, fl = sprintf('  [%.0f%% at floor: clipping?]', 100*R.dip_axes.floor_frac(i)); end
    if R.dip_axes.top_is_pulsatile(i), fl = [fl '  [top loading pulsatile]']; end %#ok<AGROW>
    fprintf('    %-5s dip=%.4f  p=%.3g  p_bonf=%.3g  top=%s%s\n', R.dip_axes.axis(i), R.dip_axes.dip(i), ...
            R.dip_axes.p(i), R.dip_axes.p_bonf(i), R.dip_axes.top_marker(i), fl);
end
if R.multimodal_pooled
    fprintf('  -> Pooled multimodality on %s: candidate discrete states. Pooled modes can simply be treatment groups; see within-treatment tests.\n', ...
            i_list(cellstr(R.dip_axes.axis(R.dip_axes.p_bonf < o.alpha))));
else
    fprintf('  -> All tested axes unimodal pooled: a functional CONTINUUM; prefer score regression over classification.\n');
end

% Within treatment
if height(R.dip_within) > 0
    sig = R.dip_within(R.dip_within.p_bonf < o.alpha, :);
    fprintf('  Within-treatment dip (%d tests, Bonferroni): %d significant.\n', height(R.dip_within), height(sig));
    for i = 1:height(sig)
        fprintf('    %-24s %-5s n=%d dip=%.4f p_bonf=%.3g\n', sig.group(i), sig.axis(i), sig.n(i), sig.dip(i), sig.p_bonf(i));
    end
    switch R.modes_source
        case 'within-treatment'
            fprintf('  -> Cells under the SAME stimulus split into modes: evidence for intrinsic functional wells, not stimulus-imposed separation.\n');
        case 'between-treatment-only'
            fprintf('  -> Modes exist only across treatments: the discreteness is imposed by the stimulus, not an intrinsic cell-state bifurcation.\n');
        otherwise
            fprintf('  -> No modes within treatments.\n');
    end
end

% Markers
fprintf('\n  Per-marker dip (Bonferroni over %d markers):\n', height(R.dip_markers));
for i = 1:height(R.dip_markers)
    fl = '';  if R.dip_markers.pulsatile(i), fl = '  [pulsatile]'; end
    if R.dip_markers.floor_frac(i) > 0.01, fl = sprintf('%s  [%.0f%% at floor]', fl, 100*R.dip_markers.floor_frac(i)); end
    fprintf('    %-14s dip=%.4f  p_bonf=%.3g%s\n', R.dip_markers.marker(i), R.dip_markers.dip(i), R.dip_markers.p_bonf(i), fl);
end

% GMM
fprintf('\n  GMM on %d kept PCs: BIC selects K=%d, ICL selects K=%d.\n', min(o.gmm_dims, nk), R.k_bic, R.k_icl);
fprintf('  -> BIC favours K>1 for any skewed or heavy-tailed unimodal cloud; ICL penalizes overlapping components. Trust K>1 only when ICL agrees AND the dip test is significant.\n');

% Pulsatile
if ~isempty(R.pulsatile)
    fprintf('\n  Pulsatile markers %s flagged: a single terminal snapshot samples a random pulse phase. Expected artefacts: inflated uniqueness (extra apparent dimension) and on/off bimodality mistaken for wells.\n', i_list(R.pulsatile));
    if ~isnan(R.nopuls.pr_frac)
        fprintf('    Without them: PR/n = %.2f (vs %.2f with).', R.nopuls.pr_frac, R.pr_frac);
        if R.pr_frac - R.nopuls.pr_frac >= 0.1
            fprintf(' Pulsatile markers add dissociated variance: consistent with phase aliasing OR genuinely independent cytokine programs; a snapshot cannot distinguish them (staggered fixation can).\n');
        elseif R.nopuls.pr_frac - R.pr_frac >= 0.1
            fprintf(' Removing them makes the space MORE dissociated: the pulsatile markers co-vary as a shared module, so phase aliasing is not inflating dimensionality here.\n');
        else
            fprintf(' Dimensionality conclusion is robust to the pulsatile markers.\n');
        end
    end
end

fprintf('\n  VERDICT: %s  (modes: %s)\n\n', R.verdict, R.modes_source);
end

function s = i_list(c)
c = cellstr(c);
if isempty(c), s = '{}'; else, s = ['{' strjoin(c, ', ') '}']; end
end


function i_brief(R, o, tag)
% Compact console summary
nk = numel(R.kept_pcs);
fprintf('\n%s %s | %d cells | %d markers | removed PC%s\n', tag, R.set_tag, R.n, numel(R.markers), mat2str(R.drop_pcs));
fprintf('  Dimensionality : PR %.1f of %d axes (PR/n %.2f; > %.2f = dissociated) | %d PCs to %.0f%% | PA keeps %d\n', ...
        R.pr, nk, R.pr_frac, o.lowdim_frac, R.n_var_target, 100*o.var_target, R.n_pa);
if ~isnan(R.m1m2_best_cos)
    vk = R.var_expl(R.kept_pcs) / sum(R.var_expl(R.kept_pcs));
    fprintf('  M1/M2 axis     : best PC%d, cos %+.2f, %.1f%% of kept variance\n', R.m1m2_best_pc, R.m1m2_best_cos, ...
            100*vk(R.kept_pcs == R.m1m2_best_pc));
end
nW = height(R.dip_within);
fprintf('  Modality       : pooled axes %d/%d multimodal | within-group %d/%d | GMM BIC K=%d, ICL K=%d\n', ...
        nnz(R.dip_axes.p_bonf < o.alpha), height(R.dip_axes), nnz(R.dip_within.p_bonf < o.alpha), nW, R.k_bic, R.k_icl);
fprintf('  -> %s (modes: %s)\n', upper(R.verdict), R.modes_source);
if ismember(1, R.drop_pcs) && ~R.pc1_uniform_sign
    fprintf('  ! PC1 has mixed-sign loadings: removing it may remove biology.\n');
end
end

%% =========================================================================
%  Step-by-step summary figure
%% =========================================================================
function figs = i_story(R, o, fig_dir)
M  = numel(R.markers);  kept = R.kept_pcs;  nk = numel(kept);
f  = figure('Color', 'w', 'Name', ['Q2 functional dimensionality ' R.set_tag], 'Units', 'pixels');
f.Position(3:4) = [1350 760];
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q2  Do discrete functional states exist?   (%d cells, %d markers)', R.n, M), ...
      'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
cB = [0.20 0.45 0.70];  cR = [0.80 0.30 0.30];

% Step 1: dimensionality
ax = nexttile(tl);  hold(ax, 'on');
lam = R.eigenvalues(kept);
bar(ax, 1:nk, lam, 0.7, 'FaceColor', cB, 'EdgeColor', 'none');
plot(ax, 1:nk, R.pa_null95(1:nk), 'k--', 'LineWidth', 1);
yyaxis(ax, 'right');  plot(ax, 1:nk, cumsum(lam)/sum(lam), 'o-', 'Color', [0.5 0.5 0.5], 'MarkerSize', 4);
ylim(ax, [0 1]);  ylabel(ax, 'cumulative variance');  ax.YAxis(2).Color = [0.5 0.5 0.5];
yyaxis(ax, 'left');  ylabel(ax, 'eigenvalue');
set(ax, 'XTick', 1:nk, 'XTickLabel', compose('PC%d', kept));
legend(ax, {'kept axes', 'shuffled-marker null 95%'}, 'Box', 'off', 'FontSize', 7.5);
i_step_title(ax, 'Step 1  How many axes', sprintf('PR = %.1f of %d (1 = one shared axis, %d = all independent)', R.pr, nk, nk));
i_sty(ax);

% Step 2: structure of the axes
ax = nexttile(tl);
ns = min(6, nk);  L = R.loadings(:, kept(1:ns));  xl = compose('PC%d', kept(1:ns));
if any(R.m1m2_contrast), L = [L, R.m1m2_contrast / norm(R.m1m2_contrast)];  xl = [xl, {'M1-M2'}]; end
cl = max(abs(L(:)));
imagesc(ax, L, [-cl cl]);  colormap(ax, i_diverging(256));  cb = colorbar(ax);  cb.Box = 'off';
for i = 1:size(L,1), for j = 1:size(L,2)
    text(ax, j, i, sprintf('%.2f', L(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 6.5, ...
         'Color', [0 0 0] + (abs(L(i,j)) > 0.6*cl));
end, end
set(ax, 'XTick', 1:size(L,2), 'XTickLabel', xl, 'XTickLabelRotation', 45, 'YTick', 1:M, 'YTickLabel', R.markers);
ttl2 = 'no M1/M2 contrast testable';
if ~isnan(R.m1m2_best_cos), ttl2 = sprintf('Best match to M1-M2: PC%d (cos %+.2f)', R.m1m2_best_pc, R.m1m2_best_cos); end
i_step_title(ax, 'Step 2  What the axes are', ttl2);
i_sty(ax);

% Step 3: modality, pooled
ax = nexttile(tl);  hold(ax, 'on');
nA = height(R.dip_axes);
plot(ax, 1:nA, -log10(R.dip_axes.p_bonf), 'o', 'MarkerFaceColor', cB, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
yline(ax, -log10(o.alpha), '--', 'Color', cR);
set(ax, 'XTick', 1:nA, 'XTickLabel', R.dip_axes.axis);  xlim(ax, [0.5 nA+0.5]);
ylim(ax, [0 max(2, max(-log10(R.dip_axes.p_bonf)) + 0.5)]);
ylabel(ax, '-log10 p (dip test, Bonferroni)');
i_step_title(ax, 'Step 3  Modes, all cells', sprintf('%d/%d axes multimodal (above red line)', nnz(R.dip_axes.p_bonf < o.alpha), nA));
i_sty(ax);

% Step 3b: modality within groups
ax = nexttile(tl);
if height(R.dip_within) > 0
    [gu, ~, gi] = unique(R.dip_within.group, 'stable');
    [au, ~, ai] = unique(R.dip_within.axis,  'stable');
    Hm = nan(numel(gu), numel(au));
    Hm(sub2ind(size(Hm), gi, ai)) = -log10(R.dip_within.p_bonf);
    imagesc(ax, Hm, [0 max(2, max(Hm(:)))]);
    colormap(ax, [linspace(1,0.80,64)' linspace(1,0.30,64)' linspace(1,0.30,64)']);  cb = colorbar(ax);  cb.Box = 'off';
    set(ax, 'XTick', 1:numel(au), 'XTickLabel', au, 'YTick', 1:numel(gu), 'YTickLabel', gu, 'FontSize', 7);
    i_step_title(ax, 'Step 3b  Modes within a group', sprintf('%d/%d tests multimodal (-log10 p > %.1f)', ...
                 nnz(R.dip_within.p_bonf < o.alpha), height(R.dip_within), -log10(o.alpha)));
else
    axis(ax, 'off');  i_step_title(ax, 'Step 3b  Modes within a group', 'not run');
end
i_sty(ax);

% Step 4: model check
ax = nexttile(tl);  hold(ax, 'on');
G = R.gmm;
plot(ax, G.K, G.BIC - G.BIC(1), 'o-', 'Color', cB, 'LineWidth', 1.2, 'MarkerFaceColor', cB);
plot(ax, G.K, G.ICL - G.ICL(1), 's-', 'Color', cR, 'LineWidth', 1.2, 'MarkerFaceColor', cR);
yline(ax, 0, ':k');
legend(ax, {sprintf('BIC (best K = %d)', R.k_bic), sprintf('ICL (best K = %d)', R.k_icl)}, 'Box', 'off', 'Location', 'best', 'FontSize', 7.5);
xlabel(ax, 'mixture components K');  ylabel(ax, 'criterion vs K = 1 (lower = better)');
i_step_title(ax, 'Step 4  Cross-check', 'Does a mixture model need > 1 cluster?');
i_sty(ax);

% Reading guide
ax = nexttile(tl);
i_textpanel(ax, { ...
  'WHY THIS ORDER', ...
  '1  Count independent axes first. If functions co-vary on a', ...
  '   few axes, labels like M1/M2 could be a fair compression;', ...
  '   if they need many axes, functions vary independently.', ...
  '2  Read what those axes are, and whether any one of them', ...
  '   matches the M1-minus-M2 contrast.', ...
  '3  Only then ask whether cells form separate groups along', ...
  '   the axes (dip test). Modes in pooled cells can just be', ...
  '   the treatments, so 3b repeats the test inside each group:', ...
  '   modes there = candidate intrinsic states.', ...
  '4  Mixture models always find clusters by BIC; ICL penalizes', ...
  '   overlap. Trust K > 1 only when ICL and the dip test agree.', ...
  '', ...
  sprintf('VERDICT  %s  (modes: %s)', upper(R.verdict), R.modes_source), ...
  sprintf('PR/n %.2f (threshold %.2f)', R.pr_frac, o.lowdim_frac), ...
  '', ...
  'NEXT  continuum -> predict marker levels (regression) in Q3;', ...
  '      discrete states -> predict state membership.'});

ctag = sprintf('%dmk', M);  if M <= 4, ctag = strjoin(R.markers, '-'); end
figs = f;  i_save(f, o, fig_dir, 'Q2-summary', ctag, sprintf('%s_%s', o.group, R.set_tag));
end

%% =========================================================================
%  Plots
%% =========================================================================
function figs = i_plots(R, P, o, rs, fig_dir)
M  = numel(R.markers);  nk = numel(R.kept_pcs);
ctag  = sprintf('%dmk', M);  if M <= 4, ctag = strjoin(R.markers, '-'); end
facet = sprintf('pooled_%s', R.set_tag);
figs  = gobjects(0);
lbl   = R.markers;  lbl(P.puls) = strcat(lbl(P.puls), ' *');

% 1 --- Scree + parallel-analysis null ------------------------------------
f  = figure('Color', 'w', 'Name', ['fdim scree ' R.set_tag]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl);  hold(ax, 'on');
bar(ax, 1:numel(R.eigenvalues), R.eigenvalues, 'FaceColor', [0.75 0.75 0.75]);
bar(ax, R.kept_pcs, R.eigenvalues(R.kept_pcs), 'FaceColor', [0.10 0.50 0.90]);
plot(ax, R.kept_pcs, R.pa_null95(1:nk), 'k--o', 'MarkerSize', 4);
legend(ax, {'removed', 'kept', 'PA null 95%'}, 'Location', 'northeast');
xlabel(ax, 'PC');  ylabel(ax, 'eigenvalue');
title(ax, sprintf('Scree (PR = %.2f / %d, PA = %d)', R.pr, nk, R.n_pa));
ax = nexttile(tl);
vk = R.eigenvalues(R.kept_pcs);  plot(ax, 1:nk, cumsum(vk)/sum(vk), 'o-', 'LineWidth', 1.5);
yline(ax, o.var_target, ':k');  ylim(ax, [0 1]);
xlabel(ax, 'kept PC rank');  ylabel(ax, 'cumulative kept variance');  title(ax, 'Cumulative');
figs(end+1) = f;  i_save(f, o, fig_dir, 'fdim-scree', ctag, facet);

% 2 --- Loading heatmap + M1/M2 contrast column ----------------------------
f  = figure('Color', 'w', 'Name', ['fdim loadings ' R.set_tag]);
ax = axes(f);
Lk = R.loadings(:, R.kept_pcs);
xl = compose('PC%d (%.0f%%)', R.kept_pcs(:), 100*R.var_expl(R.kept_pcs));
if any(P.v)
    Lk = [Lk, P.v / norm(P.v)];  xl = [xl; {'M1-M2'}];
end
imagesc(ax, Lk);  colormap(ax, i_diverging(256));
cl = max(abs(Lk(:)));  caxis(ax, [-cl cl]);  colorbar(ax);
set(ax, 'XTick', 1:size(Lk,2), 'XTickLabel', xl, 'XTickLabelRotation', 45, ...
        'YTick', 1:M, 'YTickLabel', lbl, 'TickLabelInterpreter', 'none');
title(ax, 'Loadings (kept PCs; * = pulsatile)');
figs(end+1) = f;  i_save(f, o, fig_dir, 'fdim-loadings', ctag, facet);

% 3 --- Dip tests ---------------------------------------------------------
f  = figure('Color', 'w', 'Name', ['fdim dip ' R.set_tag]);
nt = 2 + (height(R.dip_within) > 0);
tl = tiledlayout(f, 1, nt, 'TileSpacing', 'compact');
ax = nexttile(tl);
bar(ax, -log10(R.dip_axes.p));  yline(ax, -log10(o.alpha / height(R.dip_axes)), '--r', 'Bonferroni');
set(ax, 'XTickLabel', R.dip_axes.axis);  ylabel(ax, '-log10 p');  title(ax, 'Pooled axes');
ax = nexttile(tl);
bar(ax, -log10(R.dip_markers.p));  yline(ax, -log10(o.alpha / M), '--r', 'Bonferroni');
set(ax, 'XTick', 1:M, 'XTickLabel', lbl, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
title(ax, 'Markers');
if nt == 3
    ax = nexttile(tl);
    [gu, ~, gi] = unique(R.dip_within.group, 'stable');
    [au, ~, ai] = unique(R.dip_within.axis,  'stable');
    Hm = nan(numel(gu), numel(au));
    Hm(sub2ind(size(Hm), gi, ai)) = -log10(R.dip_within.p_bonf);
    imagesc(ax, Hm);  colorbar(ax);
    set(ax, 'XTick', 1:numel(au), 'XTickLabel', au, 'YTick', 1:numel(gu), 'YTickLabel', gu, ...
            'TickLabelInterpreter', 'none');
    title(ax, sprintf('Within treatment, -log10 p_{bonf} (sig > %.1f)', -log10(o.alpha)));
end
figs(end+1) = f;  i_save(f, o, fig_dir, 'fdim-dip', ctag, facet);

% 4 --- GMM model selection ----------------------------------------------
f  = figure('Color', 'w', 'Name', ['fdim gmm ' R.set_tag]);
ax = axes(f);  hold(ax, 'on');
plot(ax, R.gmm.K, R.gmm.BIC, 'o-', 'LineWidth', 1.5);
plot(ax, R.gmm.K, R.gmm.ICL, 's-', 'LineWidth', 1.5);
xline(ax, R.k_bic, ':b');  xline(ax, R.k_icl, ':r');
legend(ax, {'BIC', 'ICL', 'BIC min', 'ICL min'}, 'Location', 'best');
xlabel(ax, 'components K');  ylabel(ax, 'criterion (lower = better)');
title(ax, sprintf('GMM on %d kept PCs', min(o.gmm_dims, nk)));
figs(end+1) = f;  i_save(f, o, fig_dir, 'fdim-gmm', ctag, facet);

% 5 --- Axis densities (1D + 2D) -----------------------------------------
na = min(3, nk);  ak = R.kept_pcs(1:na);
sub = 1:size(P.score, 1);
if numel(sub) > o.max_plot_pts * 4, sub = randperm(rs, numel(sub), o.max_plot_pts * 4); end
f  = figure('Color', 'w', 'Name', ['fdim density ' R.set_tag]);
tl = tiledlayout(f, 2, na, 'TileSpacing', 'compact');
for a = 1:na
    ax = nexttile(tl);
    histogram(P.score(:, ak(a)), 100, 'Parent', ax, 'EdgeColor', 'none', 'Normalization', 'pdf');
    title(ax, sprintf('PC%d  dip p_{bonf}=%.2g', ak(a), R.dip_axes.p_bonf(a)));
end
prs = nchoosek(1:na, 2);  if na < 3, prs = prs(1, :); end
for q = 1:size(prs, 1)
    ax = nexttile(tl);
    xa = P.score(sub, ak(prs(q,1)));  ya = P.score(sub, ak(prs(q,2)));
    [Hc, xe, ye] = histcounts2(xa, ya, 80);
    xc = xe(1:end-1) + diff(xe)/2;  yc = ye(1:end-1) + diff(ye)/2;   % bin centers
    imagesc(ax, xc, yc, log10(Hc' + 1));  axis(ax, 'xy');  colormap(ax, parula);
    xlabel(ax, sprintf('PC%d', ak(prs(q,1))));  ylabel(ax, sprintf('PC%d', ak(prs(q,2))));
    title(ax, 'log10 density');
end
figs(end+1) = f;  i_save(f, o, fig_dir, 'fdim-density', ctag, facet);
end

function cm = i_diverging(n)
% Blue - white - red
h  = floor(n/2);
up = [linspace(0.13,1,h)' linspace(0.40,1,h)' linspace(0.67,1,h)'];
dn = [linspace(1,0.70,n-h)' linspace(1,0.09,n-h)' linspace(1,0.17,n-h)'];
cm = [up; dn];
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