function R = dl_dynamics_predict(T, varargin)
% dl_dynamics_predict  Do live-trace dynamics predict the CyCIF functional
%   endpoint BEYOND treatment? (1) conditional MI I(dyn; y | treatment)
%   against a within-treatment shuffle null (existence), then (2) a
%   per-target PLSR bank with nested-CV ncomp, Q2Y vs null/ceiling, block-
%   labelled VIP (only after Q2Y is established), and a horse race against
%   treatment-only and single-marker predictors (characterization).
%
% USAGE
%   dyn = struct('kind','pulse', 'X',DF);                           % LONG table from dl_make_osc_pulse_df
%   R = dl_dynamics_predict(Tf, 'dynamics',dyn, 'prenormalized',true);
%   dyn = struct('kind','latent', 'X',Zlat);                        % [height(Tf) x L], row-aligned
%   R = dl_dynamics_predict(Tf, 'dynamics',dyn, 'prenormalized',true, 'nperm',50, 'save',false);
%   dyn = struct('kind','mixed', 'X',{{Zlat, DF}}, 'block_labels',{{'latent','pulse'}}, 'pulse_flag','isPulsatile');
%   R = dl_dynamics_predict(Tf, 'dynamics',dyn, 'prenormalized',true, 'ceiling',struct('HLA_DR','HLA_DR_rep'));
%   R = dl_dynamics_predict(Tf, 'dynamics',dyn, 'prenormalized',true, 'markers',{'TNFa','IL1b','IL8','CD206'}, ...
%                           'nonlinear','knn', 'cmi',false, 'nperm',100);          % linear vs kNN, pooled + within-well
%   R = dl_dynamics_predict(A, 'dynamics',dyn, 'group','condition', 'breakdown',{'stimulus','glucose','rep'});  % Q2 per level
%   dyn = struct('kind','mixed', 'X',{{Prk, Phy}}, 'block_labels',{{'RAMPKAR2','HYLIGHT'}});
%   R = dl_dynamics_predict(A, 'dynamics',dyn, 'group','condition', 'feature_sets','blocks');      % per sensor
%   R = dl_dynamics_predict(A, 'dynamics',struct('kind','latent','X',Lat), 'group','condition', ...
%                           'feature_sets','vip', 'report_main',false, 'nperm',20);             % latent groups
%   R = dl_dynamics_predict(A, 'dynamics',dyn, 'feature_sets',struct('early',{{'latent_1','latent_4'}}, 'late',[6 12]));
%   R = dl_dynamics_predict(A, 'dynamics',dyn, 'group','condition', 'nonlinear','knn', ...
%                           'tail',{'TNFa','IL1b','IL8'}, 'tail_bg',Abg);   % rare producers: frequency + identity
%
% INPUT
%   T   table, one row per cell: functional markers + treatment + identity keys
%
% NAME-VALUE
%  Dynamics (required)
%   'dynamics'  struct with fields
%       .kind          'pulse' | 'latent' | 'mixed'
%       .X             numeric [height(T) x p] (row-aligned to T) OR table with key
%                      columns (wide, or LONG with a 'channel' column -> pivoted
%                      to <channel>_<metric>). 'mixed': cell of blocks.
%       .names         (opt) feature names (numeric X) or metric subset (table X);
%                      'mixed': cell of cells
%       .keys          (opt) {Tcol, Xcol; ...} join pairs ({'xy','xys'; 'cellid','cellidx'});
%                      add a replicate pair when T pools replicates
%       .block_labels  (opt, mixed) ({'latent','pulse'})
%       .block_weights (opt) per-block multiplier after block scaling (ones)
%       .dyn_nan       (opt) 'indicator' | 'drop' | cell per block; '' = auto
%                      (auto: 'indicator' for pulse/LONG blocks, 'drop' otherwise)
%       .exclude       (opt) table columns never used as features
%                      ({'isOsc','oscClass','oscLabel','tx_num'})
%       .pulse_flag    (opt) CALLER-SUPPLIED per-row logical, or T column name;
%                      appended as its own block 'flag'
%  Columns / scale (semantics as dl_treatment_dispersion)
%   'group' ('tx1_label')  treatment = conditioning variable Z
%   'unit' ('xy') 'rep_col' ('') 'markers' ({}) 'id_cols' 'compartment' ('')
%   'prenormalized' (false) 'transform' ('') 'cofactor' ('auto') 'exclude_spike' (true)
%   'min_n' (30)
%  Targets
%   'targets'       'markers' | 'pcs'                                   ('markers')
%   'n_target_pcs'  functional PCs used as targets for 'pcs'            (3)
%   'drop_pcs'      functional PCs removed before 'pcs' targets         ([1])
%   'canonical'     single-marker predictors for the horse race; {} = all ({})
%  Statistics
%   'nfolds'        outer CV folds (units dealt within treatment)        (5)
%   'inner_folds'   inner CV folds for ncomp                             (5)
%   'max_ncomp'     PLS component ceiling                               (10)
%   'nperm'         within-treatment linkage shuffles, full nested CV   (100)
%   'alpha'         Holm-adjusted across targets                        (0.05)
%   'ceiling'       [] | numeric per target | struct(target -> value or
%                   repeat-measurement column name -> test-retest r)    ([])
%   'cmi'           run conditional MI                                  (true)
%   'cmi_dims'      dynamics PCs used for MI                            (5)
%   'cmi_k'         KSG neighbours                                      (5)
%   'cmi_max_n'     per-stratum subsample for MI                        (1000)
%   'cmi_nperm'     MI within-treatment shuffles                        (100)
%   'nonlinear'     'none' | 'knn': add a kNN regressor beside PLSR      ('none')
%   'knn_k'         candidate neighbour counts, chosen by inner CV       ([10 25 50 100 200])
%   'nperm_nl'      kNN nulls (y shuffled within treatment); [] = nperm  ([])
%   'seed'          (0)
%  Output
%   'save' (true) 'closefigs' (true) 'dataloc' ([]) 'savepath' ('') 'max_plot_pts' (3000)
%   'plots'     'summary' = step-by-step figure + one bar figure per breakdown variable |
%               'all' = + diagnostics | 'none'                                     ('summary')
%   'verbose'   'summary' = one table + verdicts | 'full' = every value interpreted | 'none' ('summary')
%   'breakdown' cell  columns of T: Q2 is reported within each level of each column
%               (e.g. {'stimulus','glucose','rep'}); {} = {group, replicate}      ({})
%   'n_boot'    well-level bootstrap draws for breakdown CIs                       (200)
%  Feature sets (Step 4b: which parts of the dynamics carry the prediction?)
%   'feature_sets'  {} | 'blocks' | 'vip' | struct(name -> feature names / block labels / column idx)  ({})
%               'blocks' = one set per dynamics block (e.g. one per sensor)
%               'vip'    = dynamics columns clustered by their VIP profile across targets,
%                          computed on DISCOVERY wells; Q2 then measured on the other
%                          (VALIDATION) wells, so set choice does not inflate Q2
%               For each set: ALONE = Q2 from that set only; UNIQUE = Q2(all) - Q2(all but set)
%   'vip_k'         number of VIP groups; [] = best silhouette over 2..6       ([])
%   'set_nperm'     within-treatment shuffles per set for the null             (50)
%   'reference'     struct('label',..,'table',T(target,q2)): dotted reference line per target in
%               the sets figure (e.g. the latent Q2 from another run)            ([])
%   'report_main'   false = skip Steps 1-4 console/figures (sets only)          (true)
%  Rare producers (Step T: markers made by few, very bright cells; targets = 'markers' only)
%   'tail'          {} | 'all' | target names. The regression stays CONTINUOUS; the gate is
%                   used only to EVALUATE it: T1 does treatment set producer FREQUENCY
%                   (well = unit)? T2 do out-of-fold predictions rank producers above
%                   non-producers of the SAME well (stratified within-well AUC)?       ({})
%   'tail_k'        gates: median + k x robust SD (1.4826 MAD) of the stained cells, per
%                   replicate; swept for sensitivity                                ([3 5 8])
%   'tail_kref'     k used for frequency, co-positivity, verdicts and the figure    (5)
%   'tail_bg'       table of secondary-only cells ON THE SAME SCALE as T (same marker
%                   columns + replicate column): false-positive rate of each gate   ([])
%   'tail_nperm'    within-well score shuffles for the AUC null; well-label
%                   shuffles for the frequency test                                 (500)
%
% OUTPUT
%   R struct:
%     .targets .blocks (table label | n_features | weight | nan_policy) .features
%     .feature_block .n_rows
%     .cmi    table target | n | cmi | mi_pooled | null_med | null95 | excess | p | p_holm
%     .q2     table target | n | q2_dyn | q2_null_med | q2_null95 | p_dyn | p_dyn_holm | ncomp |
%                   q2_tx | q2_txdyn | dq2 | p_delta | p_delta_holm | ceiling | q2_frac_ceiling |
%                   best_marker | q2_best_marker | verdict
%     .vip    struct array per target: target | established | vip | block_share (table)
%     .single_q2   [nPredictor x nTarget] single-marker Q2 ; .single_names
%     .breakdown struct array per variable: var | table target | level | n | n_wells | q2 |
%                ci_lo | ci_hi | null95 | above_null   (within-level Q2, PLS, outer CV)
%     .sets   struct: mode | names | members (feature names) | table (target | set | n | q2_all |
%                q2_alone | alone_lo | alone_hi | alone_null95 | alone_above | q2_unique | unique_lo |
%                unique_hi | unique_null95 | unique_above) | vip (discovery VIP, feature x target) |
%                n_disc_wells | n_val_wells
%     .pred   struct array per target: y | yhat (outer-CV PLS) | yhat_nl (outer-CV kNN) | G | W (well id)
%     .resolution  struct  Holm floor per test family; nperm needed to resolve
%     .nl     table  target | q2_lin | q2_lin_null95 | q2_knn | q2_knn_null95 | p_knn_holm | knn_k |
%                    q2w_lin | q2w_lin_null95 | p_w_lin_holm | q2w_knn | q2w_knn_null95 | p_w_knn_holm
%                    (q2w = within-well Q2: observed and predicted centred on each well's mean)
%     .tail   ('tail' only) struct: k | kref | targets | has_bg |
%             table  (target | k | n | pct_tail | bg_fpr | auc_tx | auc_lin | auc_lin_w | auc_lin_w_null95 |
%                     p_lin_w | p_lin_w_holm | auc_knn | auc_knn_w | auc_knn_w_null95 | p_knn_w | p_knn_w_holm)
%                     auc_tx = leave-well-out treatment frequency (pooled); *_w = within-well, stratified
%             freq   (var | target | level | n_wells | pct_mean | pct_sd | pct_lo | pct_hi), kref, well = unit
%             freq_effect (var | target | eta2 | p): share of well-to-well frequency variance set by var
%             copos  struct names | is_tail | enrich (P(a&b)/P(a)P(b)) | n_both, kref, all targets
%             verdict table (target | gate | frequency | identity | verdict)
%             label  [n_rows x nTail] producer call at kref (NaN = not measured), rows = .keep_idx
%     .keep_idx  rows of T analysed (R.pred / R.tail.label row space = T(keep_idx,:) )
%     .figs
%
% NOTES
%   New in DatalocHandler v5 (no v4 counterpart). Identity: T.(xy,cellid) and
%   dl_make_osc_pulse_df (xys,cellidx) both carry dataloc.d{xy}.cellindex, so the
%   default join is label-to-label; pooled replicates REQUIRE a replicate key pair.

%% ---- Parse inputs -------------------------------------------------------
p = inputParser; p.FunctionName = mfilename;
isTxt = @(x) ischar(x) || isstring(x);
addRequired (p, 'T', @istable);
addParameter(p, 'dynamics',      [],          @isstruct);
addParameter(p, 'group',         'tx1_label', isTxt);
addParameter(p, 'unit',          'xy',        isTxt);
addParameter(p, 'rep_col',       '',          isTxt);
addParameter(p, 'markers',       {},          @(x) iscell(x) || isstring(x));
addParameter(p, 'id_cols',       {'xy','cellid','ptx_is_spike','rep','ExpID','xys','cellidx'});
addParameter(p, 'compartment',   '',          isTxt);
addParameter(p, 'prenormalized', false,       @(x) islogical(x) || isnumeric(x));
addParameter(p, 'transform',     '',          isTxt);
addParameter(p, 'cofactor',      'auto');
addParameter(p, 'exclude_spike', true);
addParameter(p, 'min_n',         30);
addParameter(p, 'targets',       'markers',   @(x) any(strcmpi(x, {'markers','pcs'})));
addParameter(p, 'n_target_pcs',  3);
addParameter(p, 'drop_pcs',      1,           @isnumeric);
addParameter(p, 'canonical',     {},          @(x) iscell(x) || isstring(x));
addParameter(p, 'nfolds',        5);
addParameter(p, 'inner_folds',   5);
addParameter(p, 'max_ncomp',     10);
addParameter(p, 'nperm',         100);
addParameter(p, 'alpha',         0.05);
addParameter(p, 'ceiling',       []);
addParameter(p, 'cmi',           true);
addParameter(p, 'cmi_dims',      5);
addParameter(p, 'cmi_k',         5);
addParameter(p, 'cmi_max_n',     1000);
addParameter(p, 'cmi_nperm',     100);
addParameter(p, 'seed',          0);
addParameter(p, 'save',          true);
addParameter(p, 'closefigs',     true);
addParameter(p, 'dataloc',       []);
addParameter(p, 'savepath',      '',          isTxt);
addParameter(p, 'max_plot_pts',  3000);
addParameter(p, 'plots',         'summary',   @(x) any(strcmpi(x, {'summary','all','none'})));
addParameter(p, 'verbose',       'summary',   @(x) any(strcmpi(x, {'summary','full','none'})));
addParameter(p, 'breakdown',     {},          @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(p, 'n_boot',        200);
addParameter(p, 'feature_sets',  {});
addParameter(p, 'vip_k',         []);
addParameter(p, 'set_nperm',     50);
addParameter(p, 'reference',     []);
addParameter(p, 'report_main',   true);
addParameter(p, 'nonlinear',     'none',      @(x) any(strcmpi(x, {'none','knn'})));
addParameter(p, 'knn_k',         [10 25 50 100 200], @isnumeric);
addParameter(p, 'nperm_nl',      []);
addParameter(p, 'tail',          {});
addParameter(p, 'tail_k',        [3 5 8],     @isnumeric);
addParameter(p, 'tail_kref',     5,           @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'tail_bg',       [],          @(x) isempty(x) || istable(x));
addParameter(p, 'tail_nperm',    500);
parse(p, T, varargin{:});
o = p.Results;
o.group = char(o.group);  o.unit = char(o.unit);  o.rep_col = char(o.rep_col);
o.compartment = char(o.compartment);  o.transform = lower(char(o.transform));
o.targets = lower(char(o.targets));   o.markers = cellstr(o.markers);
o.canonical = cellstr(o.canonical);   o.prenormalized = logical(o.prenormalized);
o.nonlinear = lower(char(o.nonlinear));  if isempty(o.nperm_nl), o.nperm_nl = o.nperm; end
o.plots = lower(char(o.plots));  o.breakdown = cellstr(o.breakdown);
o.vb = find(strcmpi(o.verbose, {'none','summary','full'})) - 1;   % 0 none, 1 summary, 2 full
o.tail_k = unique([o.tail_k(:)' o.tail_kref]);                     % kref always evaluated
tailOn = ~isempty(o.tail) && ~(ischar(o.tail) && any(strcmpi(o.tail, {'', 'none'})));
assert(~tailOn || strcmpi(o.targets, 'markers'), '[dl_dynamics_predict] ''tail'' needs targets = ''markers''.');
rs  = RandStream('mt19937ar', 'Seed', o.seed);
tag = '[dl_dynamics_predict]';
vn  = T.Properties.VariableNames;
assert(~isempty(o.dynamics), '%s ''dynamics'' struct is required.', tag);
assert(ismember(o.group, vn), '%s group column ''%s'' not in table.', tag, o.group);

%% ---- Row filtering ------------------------------------------------------
keep = true(height(T), 1);
if o.exclude_spike && ismember('ptx_is_spike', vn)
    ns = nnz(T.ptx_is_spike);
    if ns > 0 && o.vb >= 2, fprintf('%s excluding %d spike-treatment cells.\n', tag, ns); end
    keep = keep & ~logical(T.ptx_is_spike);
end
g    = string(T.(o.group));
keep = keep & ~(ismissing(g) | strlength(strtrim(g)) == 0);

rep_col = o.rep_col;
if isempty(rep_col)
    cand = {'rep','ExpID','exp'};  hit = cand(ismember(cand, vn));
    if ~isempty(hit), rep_col = hit{1}; end
end
u = [];
if ~isempty(o.unit) && ismember(o.unit, vn)
    u = string(T.(o.unit));
    if ~isempty(rep_col), u = string(T.(rep_col)) + "|" + u; end
else
    warning('%s unit column absent: CV folds and nulls fall back to cells (optimistic Q2 via shared well effects).', tag);
end

%% ---- Dynamics blocks (aligned to T rows) --------------------------------
B = i_build_dyn(T, o.dynamics, tag, o.vb);
[D, dnames, dblock, okD, Btab] = i_dyn_matrix(B, keep, tag, o.vb);
keep = keep & okD;
if o.vb >= 2
    fprintf('%s %d cells with treatment + dynamics; %d dynamics features in %d blocks.\n', ...
            tag, nnz(keep), size(D,2), height(Btab));
end
D = D(keep, :);

%% ---- Functional targets -------------------------------------------------
[mk, mk_disp] = i_resolve_markers(T, o, rep_col, tag);
if o.prenormalized
    if ~isempty(o.transform) && ~strcmp(o.transform, 'linear')
        warning('%s prenormalized=true: forcing ''linear''.', tag);
    end
    xf = 'linear';
elseif isempty(o.transform), xf = 'arcsinh';
else,                        xf = o.transform;
end
[Xt, xlab] = i_transform(double(T{keep, mk}), xf, o, tag);
gk = g(keep);  uk = [];  if ~isempty(u), uk = u(keep); end

switch o.targets
    case 'markers'
        Y = Xt;  tnames = mk_disp;
    case 'pcs'
        cc = all(~isnan(Xt), 2);
        if o.vb >= 2, fprintf('%s ''pcs'' targets use %d/%d complete-case cells.\n', tag, nnz(cc), numel(cc)); end
        Zc = (Xt(cc,:) - mean(Xt(cc,:),1)) ./ std(Xt(cc,:),0,1);
        [~, sc] = pca(Zc);
        kp = setdiff(1:size(sc,2), o.drop_pcs);  kp = kp(1:min(o.n_target_pcs, numel(kp)));
        Y  = nan(numel(cc), numel(kp));  Y(cc, :) = sc(:, kp);
        tnames = arrayfun(@(k) sprintf('fPC%d', k), kp, 'UniformOutput', false);
end
nT = numel(tnames);

% Canonical single-marker predictors (markers targets only)
if isempty(o.canonical), can = mk_disp; else, can = o.canonical(:)'; end
[~, ican] = ismember(i_norm(can), i_norm(mk_disp));
can = mk_disp(ican(ican > 0));  ican = ican(ican > 0);

ceil_v = i_resolve_ceiling(o.ceiling, tnames, T, keep, xf, o, Xt, tag);
fig_dir = i_figdir(o);
if o.save && isempty(fig_dir), warning('%s save=true but no fig dir/savepath: not saved.', tag); end

%% ---- Breakdown variables (Step 3: where is it predictable?) ---------------
bdv = o.breakdown(~cellfun(@isempty, o.breakdown));
if isempty(bdv), bdv = {o.group};  if ~isempty(rep_col), bdv{end+1} = rep_col; end; end
okb = ismember(bdv, vn);
if any(~okb), warning('%s breakdown columns not in table (skipped): %s', tag, strjoin(bdv(~okb), ', ')); end
bdv = bdv(okb);  nbd = numel(bdv);
BL  = cell(1, nbd);
for i = 1:nbd, BL{i} = string(T.(bdv{i})(keep));  BL{i}(ismissing(BL{i})) = "<missing>"; end
if o.vb >= 1
    fprintf('%s %d cells | %d targets | dynamics: %s | conditioned on %s | CV on held-out wells\n', tag, ...
            nnz(keep), nT, strjoin(compose('%s (%d)', Btab.label, Btab.n_features), ' + '), o.group);
end

%% ---- Per-target analysis -------------------------------------------------
res = cell(nT, 1);
for t = 1:nT
    rows = ~isnan(Y(:, t));
    if nnz(rows) < o.min_n * 2
        warning('%s %s: only %d cells; skipped.', tag, tnames{t}, nnz(rows));  continue
    end
    t0 = tic;
    y  = Y(rows, t);  y = (y - mean(y)) / std(y);
    Dt = D(rows, :);
    [Gt, ~] = findgroups(gk(rows));
    Ut = [];
    if ~isempty(uk)
        Ut = findgroups(uk(rows));
        if any(accumarray(Ut, Gt, [], @(x) numel(unique(x))) > 1)
            warning('%s %s: units span >1 treatment; cell-level folds.', tag, tnames{t});  Ut = [];
        end
    end

    r = struct('target', tnames{t}, 'n', numel(y), 'rows', rows, 'tcol', t);

    % --- (1) Conditional MI: existence test ------------------------------
    r.cmi = NaN; r.mi_pooled = NaN; r.cmi_null = NaN;
    if o.cmi
        [~, sc] = pca(Dt);
        q  = min(o.cmi_dims, size(sc, 2));
        Dp = sc(:, 1:q) ./ std(sc(:, 1:q), 0, 1);       % unit-variance PCs for max-norm
        [r.cmi, r.cmi_null, r.mi_pooled] = i_cmi(Dp, y, Gt, o, rs);
    end

    % --- (2) PLSR bank: nested CV, fixed folds for all models/permutations ---
    f   = i_folds(Gt, Ut, o.nfolds, rs);
    fin = cell(max(f), 1);
    for k = 1:max(f)
        tr = f ~= k;
        Ui = [];  if ~isempty(Ut), Ui = Ut(tr); end
        fin{k} = i_folds(Gt(tr), Ui, o.inner_folds, rs);
    end
    OH = dummyvar(Gt);
    if size(OH, 2) < 2
        OH = zeros(numel(y), 0);                                   % single treatment: no tx block
    else
        OH = (OH - mean(OH,1)) ./ std(OH,0,1) / sqrt(size(OH,2)); % block-scaled treatment
    end

    [yh_dyn, kch]  = i_nested_pls(Dt, y, f, fin, o.max_ncomp);
    r.q2_dyn       = i_q2(y, yh_dyn, f);
    r.q2_tx        = i_q2(y, i_cv_groupmean(Gt, y, f), f);
    r.q2_txdyn     = i_q2(y, i_nested_pls([Dt OH], y, f, fin, o.max_ncomp), f);
    r.ncomp        = mode(kch(~isnan(kch)));

    % Honest null: shuffle live<->fixed linkage WITHIN treatment, rerun nested CV
    Lt = cell(1, nbd);  Lu = cell(1, nbd);  qbn = cell(1, nbd);
    for i = 1:nbd
        [Lu{i}, ~, Lt{i}] = unique(BL{i}(rows));
        qbn{i} = nan(o.nperm, numel(Lu{i}));
    end
    qn = nan(o.nperm, 1);  dn = nan(o.nperm, 1);  qwn = nan(o.nperm, 1);
    for b = 1:o.nperm
        pp    = i_perm_within(Gt, rs);
        yb    = i_nested_pls(Dt(pp,:), y, f, fin, o.max_ncomp);
        qn(b) = i_q2(y, yb, f);
        if ~isempty(Ut), qwn(b) = i_q2w(y, yb, Ut); end
        for i = 1:nbd, qbn{i}(b, :) = i_q2_levels(y, yb, Lt{i}, numel(Lu{i})); end
        dn(b) = i_q2(y, i_nested_pls([Dt(pp,:) OH], y, f, fin, o.max_ncomp), f) - r.q2_tx;
    end
    r.q2_null = qn;
    r.p_dyn   = (1 + nnz(qn >= r.q2_dyn)) / (o.nperm + 1);
    r.dq2     = r.q2_txdyn - r.q2_tx;
    r.p_delta = (1 + nnz(dn >= r.dq2)) / (o.nperm + 1);

    % --- Step 3 breakdown: within-level Q2, well-bootstrap CI, null 95th ------
    r.bd = cell(1, nbd);
    for i = 1:nbd
        [qv, nv, nwv, lo, hi] = i_q2_levels_boot(y, yh_dyn, Lt{i}, numel(Lu{i}), Ut, o.n_boot, rs);
        n95 = prctile(qbn{i}, 95, 1)';
        r.bd{i} = table(repmat(string(tnames{t}), numel(Lu{i}), 1), Lu{i}(:), nv, nwv, qv, lo, hi, n95, qv > n95, ...
            'VariableNames', {'target','level','n','n_wells','q2','ci_lo','ci_hi','null95','above_null'});
    end

    % --- Within-well Q2 (removes well means from observed AND predicted) -----
    r.q2w_lin = NaN;  r.q2w_lin_null = qwn;  r.p_w_lin = NaN;
    if ~isempty(Ut)
        r.q2w_lin = i_q2w(y, yh_dyn, Ut);
        r.p_w_lin = (1 + nnz(qwn >= r.q2w_lin)) / (o.nperm + 1);
    end

    % --- Nonlinear comparator: kNN regression, k by inner CV ------------------
    %  Same outer/inner folds as PLSR. Null shuffles y within treatment with the
    %  dynamics fixed: equivalent to the linkage shuffle, and lets neighbour lists
    %  be computed once (each permutation is then pure indexing).
    [r.q2_nl, r.p_nl, r.q2w_nl, r.p_w_nl, r.knn_k] = deal(NaN);
    [r.q2_nl_null, r.q2w_nl_null] = deal(nan(o.nperm_nl, 1));
    yh_nl = nan(size(y));
    if strcmp(o.nonlinear, 'knn')
        KP = i_knn_prep(Dt, f, fin, max(o.knn_k));
        [yh_nl, kk] = i_knn_predict(KP, y, f, fin, o.knn_k);
        r.q2_nl = i_q2(y, yh_nl, f);  r.knn_k = mode(kk(~isnan(kk)));
        if ~isempty(Ut), r.q2w_nl = i_q2w(y, yh_nl, Ut); end
        for b = 1:o.nperm_nl
            yp  = y(i_perm_within(Gt, rs));
            yhb = i_knn_predict(KP, yp, f, fin, o.knn_k);
            r.q2_nl_null(b) = i_q2(yp, yhb, f);
            if ~isempty(Ut), r.q2w_nl_null(b) = i_q2w(yp, yhb, Ut); end
        end
        r.p_nl = (1 + nnz(r.q2_nl_null >= r.q2_nl)) / (o.nperm_nl + 1);
        if ~isempty(Ut), r.p_w_nl = (1 + nnz(r.q2w_nl_null >= r.q2w_nl)) / (o.nperm_nl + 1); end
        clear KP
    end

    % --- Horse race: single canonical markers (co-stained rows only) --------
    r.single = nan(numel(can), 1);
    if strcmp(o.targets, 'markers')
        allrows = find(rows);
        for c = 1:numel(can)
            if strcmp(can{c}, tnames{t}), continue; end
            x  = Xt(allrows, ican(c));  ok = ~isnan(x);
            if nnz(ok) < o.min_n * 2, continue; end
            xs = (x(ok) - mean(x(ok))) / std(x(ok));
            Uc = [];  if ~isempty(Ut), Uc = Ut(ok); end
            fc = i_folds(Gt(ok), Uc, o.nfolds, rs);
            r.single(c) = i_q2(y(ok), i_cv_ols(xs, y(ok), fc), fc);
        end
    end

    % --- VIP on the full data at the CV-chosen ncomp -----------------------
    r.vip = i_vip(Dt, y, r.ncomp);
    r.y = y;  r.yhat = yh_dyn;  r.yhat_nl = yh_nl;  r.G = Gt;  r.W = Ut;
    res{t} = r;
    if o.vb >= 2
        fprintf('%s %-10s n=%6d  Q2dyn=%+.3f  Q2tx=%+.3f  Q2tx+dyn=%+.3f  CMI=%.3f', tag, ...
                tnames{t}, r.n, r.q2_dyn, r.q2_tx, r.q2_txdyn, r.cmi);
        if strcmp(o.nonlinear, 'knn'), fprintf('  Q2knn=%+.3f  within-well lin/knn=%+.3f/%+.3f', r.q2_nl, r.q2w_lin, r.q2w_nl); end
        fprintf('  (%.0fs)\n', toc(t0));
    elseif o.vb == 1
        fprintf('  %d/%d %s done (%.0f s)\n', t, nT, tnames{t}, toc(t0));
    end
end
done = ~cellfun(@isempty, res);
assert(any(done), '%s no target analyzable.', tag);
res  = [res{done}];  tnames = tnames(done);  ceil_v = ceil_v(done);

%% ---- Step 4b: feature sets (alone / unique) ----------------------------------
SR = [];
fsOn = isstruct(o.feature_sets) || ((ischar(o.feature_sets) || isstring(o.feature_sets)) && ~any(strcmpi(o.feature_sets, {'', 'none'})));
if fsOn
    SR = i_feature_sets(res, Y, D, dnames, dblock, Btab, gk, uk, o, rs, tag);
end

%% ---- Assemble results ---------------------------------------------------
R = i_assemble(res, tnames, ceil_v, can, dnames, dblock, Btab, o);
R.transform_label = xlab;
R.n_rows = nnz(keep);
R.breakdown = struct('var', bdv, 'table', []);
for i = 1:nbd
    tb = arrayfun(@(rr) rr.bd{i}, res, 'UniformOutput', false);
    R.breakdown(i).table = vertcat(tb{:});
end
R.sets = SR;
R.keep_idx = find(keep);

%% ---- Step T: rare producers (frequency set by treatment? identity by dynamics?) --
R.tail = [];
if tailOn
    repk = repmat("all", nnz(keep), 1);
    if ~isempty(rep_col), repk = string(T.(rep_col)(keep)); end
    R.tail = i_tail(res, Y, mk, repk, rep_col, uk, bdv, BL, xf, o, rs, tag);
end

if o.report_main
    if o.vb == 2, i_report(R, o, tag); elseif o.vb == 1, i_brief(R, o, tag); end
end
if ~isempty(SR) && o.vb >= 1, i_sets_brief(SR, tag); end
if ~isempty(R.tail) && o.vb >= 1, i_tail_brief(R.tail, o, tag); end
R.figs = gobjects(0);
if ~strcmp(o.plots, 'none')
    % A plotting error must not discard hours of CV: warn and return R without the figure.
    try
        if o.report_main
            R.figs = i_story(R, o, fig_dir);
            R.figs = [R.figs, i_breakdown_figs(R, o, fig_dir)];
        end
        if ~isempty(SR), R.figs = [R.figs, i_sets_fig(SR, o, fig_dir)]; end
        if ~isempty(R.tail), R.figs = [R.figs, i_tail_fig(R.tail, o, fig_dir)]; end
        if strcmp(o.plots, 'all'), R.figs = [R.figs, i_plots(R, o, rs, fig_dir)]; end
    catch ME
        warning('%s plotting failed (%s, line %d of %s); results are returned without the remaining figures.', ...
                tag, ME.message, ME.stack(1).line, ME.stack(1).name);
    end
    R.figs = R.figs(isgraphics(R.figs));
    if o.closefigs && o.save && ~isempty(fig_dir), close(R.figs);  R.figs = gobjects(0); end
end
end


%% =========================================================================
%  Dynamics assembly
%% =========================================================================
function B = i_build_dyn(T, dyn, tag, vb)
assert(isfield(dyn, 'kind') && isfield(dyn, 'X'), '%s dynamics needs fields kind and X.', tag);
kind = lower(char(dyn.kind));
switch kind
    case 'pulse',  Xs = {dyn.X};  lab = {'pulse'};
    case 'latent', Xs = {dyn.X};  lab = {'latent'};
    case 'mixed'
        assert(iscell(dyn.X) && numel(dyn.X) >= 2, '%s mixed: X must be a cell of >=2 blocks.', tag);
        Xs  = dyn.X;
        lab = i_getf(dyn, 'block_labels', {'latent','pulse'});
        for b = numel(lab)+1:numel(Xs), lab{b} = sprintf('block%d', b); end
    otherwise
        error('%s dynamics.kind must be pulse | latent | mixed.', tag);
end
names = i_getf(dyn, 'names', {});
keys  = i_getf(dyn, 'keys', {'xy','xys'; 'cellid','cellidx'});
wts   = i_getf(dyn, 'block_weights', ones(1, numel(Xs)));
nanp  = i_getf(dyn, 'dyn_nan', '');
excl  = i_getf(dyn, 'exclude', {'isOsc','oscClass','oscLabel','tx_num'});

B = struct('label', {}, 'X', {}, 'names', {}, 'weight', {}, 'nan', {});
for b = 1:numel(Xs)
    nb = {};
    if strcmp(kind, 'mixed')
        if numel(names) >= b && iscell(names{b}), nb = names{b}; end
    else
        nb = names;
    end
    [Xb, nb, islong] = i_align_block(T, Xs{b}, cellstr(nb), keys, excl, lab{b}, tag, vb);
    np = nanp;  if iscell(np), np = np{b}; end
    if isempty(np)
        if strcmp(lab{b}, 'pulse') || islong, np = 'indicator'; else, np = 'drop'; end
    end
    B(end+1) = struct('label', lab{b}, 'X', Xb, 'names', {nb}, 'weight', wts(b), 'nan', np); %#ok<AGROW>
end

% Caller-supplied pulsatility flag (provenance owned by the caller)
if isfield(dyn, 'pulse_flag') && ~isempty(dyn.pulse_flag)
    pf = dyn.pulse_flag;
    if ischar(pf) || isstring(pf)
        assert(ismember(char(pf), T.Properties.VariableNames), '%s pulse_flag column ''%s'' not in T.', tag, char(pf));
        pf = T.(char(pf));
    end
    assert(numel(pf) == height(T), '%s pulse_flag needs one entry per row of T.', tag);
    B(end+1) = struct('label', 'flag', 'X', double(pf(:)), 'names', {{'pulse_flag'}}, ...
                      'weight', i_getf(dyn, 'flag_weight', 1), 'nan', 'drop');
end
end

function [Xb, nb, islong] = i_align_block(T, Xin, nb, keys, excl, lab, tag, vb)
% Return block features row-aligned to T (NaN where a cell has no dynamics)
nT = height(T);  islong = false;
if isnumeric(Xin)
    assert(size(Xin,1) == nT, ['%s %s block: numeric X has %d rows, T has %d. Pass a table ', ...
           'with key columns to join instead.'], tag, lab, size(Xin,1), nT);
    Xb = double(Xin);
    if isempty(nb), nb = arrayfun(@(i) sprintf('%s_%d', lab, i), 1:size(Xb,2), 'UniformOutput', false); end
    assert(numel(nb) == size(Xb,2), '%s %s block: names/columns mismatch.', tag, lab);
    return
end
assert(istable(Xin), '%s %s block: X must be numeric or a table.', tag, lab);
xv = Xin.Properties.VariableNames;  tv = T.Properties.VariableNames;

% Resolve key pairs present in both tables
tk = {};  xk = {};
for i = 1:size(keys, 1)
    tc = keys{i,1};  xc = keys{i,2};
    if ~ismember(xc, xv) && ismember(tc, xv), xc = tc; end
    if ismember(tc, tv) && ismember(xc, xv), tk{end+1} = tc; xk{end+1} = xc; end %#ok<AGROW>
end
assert(~isempty(tk), '%s %s block: no key pair present in both tables.', tag, lab);
kT = i_key(T, tk);
assert(numel(unique(kT)) == nT, ['%s %s block: keys {%s} are not unique in T (pooled replicates?). ', ...
       'Add a replicate pair to dynamics.keys, e.g. {''rep'',''exp''}.'], tag, lab, strjoin(tk, ','));
kX = i_key(Xin, xk);

isnum = varfun(@(x) isnumeric(x) && ~islogical(x), Xin, 'OutputFormat', 'uniform');
feat  = xv(isnum & ~ismember(xv, [xk, excl(:)']));
if ~isempty(nb), feat = intersect(feat, nb, 'stable'); end
assert(~isempty(feat), '%s %s block: no numeric feature columns.', tag, lab);

islong = ismember('channel', xv);
if islong
    ch = string(Xin.channel);
    [ku, ~, ki] = unique(kX);  [cu, ~, ci] = unique(ch);
    cnt = accumarray([ki ci], 1, [numel(ku) numel(cu)]);
    assert(all(cnt(:) <= 1), '%s %s block: duplicate (cell, channel) rows in LONG table.', tag, lab);
    Xw = [];  names = {};
    for i = 1:numel(feat)
        Wm = nan(numel(ku), numel(cu));
        Wm(sub2ind(size(Wm), ki, ci)) = double(Xin.(feat{i}));
        Xw = [Xw Wm]; %#ok<AGROW>
        names = [names, strcat(cellstr(cu)', '_', feat{i})]; %#ok<AGROW>
    end
    kX = ku;
else
    assert(numel(unique(kX)) == numel(kX), '%s %s block: duplicate keys in dynamics table.', tag, lab);
    Xw = double(Xin{:, feat});  names = feat;
end
[tf, loc] = ismember(kT, kX);
Xb = nan(nT, size(Xw, 2));
Xb(tf, :) = Xw(loc(tf), :);
nb = names;
if mean(tf) < 0.5
    warning('%s %s block: only %.0f%% of T rows matched; check key formats (e.g. replicate labels).', tag, lab, 100*mean(tf));
end
if vb >= 1
    fprintf('%s %s block: %d features, %d/%d T rows matched (%.0f%%) via {%s}<->{%s}.\n', tag, lab, ...
            size(Xb,2), nnz(tf), nT, 100*mean(tf), strjoin(tk, ','), strjoin(xk, ','));
end
end

function [D, dnames, dblock, ok, Btab] = i_dyn_matrix(B, keep, tag, vb)
% Row mask, NaN policy, per-column z-score, block scaling (1/sqrt(p_b) x weight)
ok = keep;
for b = 1:numel(B)
    ok = ok & any(~isnan(B(b).X), 2);                         % cell matched to this block
    if strcmp(B(b).nan, 'drop'), ok = ok & all(~isnan(B(b).X), 2); end
end
if vb >= 2, fprintf('%s %d/%d kept cells have all dynamics blocks.\n', tag, nnz(ok), nnz(keep)); end
D = zeros(nnz(ok), 0);  dnames = {};  dblock = [];
lab = cell(numel(B),1);  nf = zeros(numel(B),1);  wt = zeros(numel(B),1);  np = cell(numel(B),1);
for b = 1:numel(B)
    Xb = B(b).X(ok, :);  nb = B(b).names;
    if strcmp(B(b).nan, 'indicator')
        miss = isnan(Xb);  hasm = any(miss, 1);
        if any(hasm)
            if vb >= 2, fprintf('%s %s: undefined values filled (median) + indicator for %s\n', tag, B(b).label, strjoin(nb(hasm), ', ')); end
            med = median(Xb, 1, 'omitnan');
            for j = find(hasm), Xb(miss(:,j), j) = med(j); end
            I = double(miss(:, hasm));  In = strcat(nb(hasm), '_undef');
            [~, iu] = unique(I', 'rows', 'stable');              % dedupe identical patterns
            Xb = [Xb I(:, iu)];  nb = [nb In(iu)]; %#ok<AGROW>
        end
    end
    sd   = std(Xb, 0, 1);  flat = ~(sd > 0);
    if any(flat)
        if vb >= 1, fprintf('%s %s: dropping %d zero-variance features.\n', tag, B(b).label, nnz(flat)); end
        Xb(:, flat) = [];  nb(flat) = [];  sd(flat) = [];
    end
    pb = size(Xb, 2);
    Zb = (Xb - mean(Xb, 1)) ./ sd * (B(b).weight / sqrt(pb));
    D  = [D Zb];  dnames = [dnames nb];  dblock = [dblock repmat(b, 1, pb)]; %#ok<AGROW>
    lab{b} = B(b).label;  nf(b) = pb;  wt(b) = B(b).weight;  np{b} = B(b).nan;
end
Dfull = nan(numel(ok), size(D,2));  Dfull(ok, :) = D;  D = Dfull;
Btab = table(string(lab), nf, wt, string(np), 'VariableNames', {'label','n_features','weight','nan_policy'});
end


%% =========================================================================
%  Statistics
%% =========================================================================
function [cmi, nullv, mi_pool] = i_cmi(Dp, y, G, o, rs)
% I(dyn; y | Z) = sum_z p(z) I_z (KSG alg. 1, max-norm); null: shuffle y within z
Dp = Dp + 1e-10 * randn(rs, size(Dp));  y = y + 1e-10 * randn(rs, size(y));   % break ties
Dx = {};  Dy = {};  w = [];
for gg = unique(G)'
    idx = find(G == gg);  ng = numel(idx);
    if ng < max(o.min_n, o.cmi_k + 2), continue; end
    if ng > o.cmi_max_n, idx = idx(randperm(rs, ng, o.cmi_max_n)); end
    Dx{end+1} = i_cheb(Dp(idx, :));  Dy{end+1} = i_cheb(y(idx));  w(end+1) = ng; %#ok<AGROW>
end
cmi = NaN;  nullv = nan(o.cmi_nperm, 1);  mi_pool = NaN;
if isempty(w), return; end
w   = w / sum(w);
cmi = 0;
for s = 1:numel(w), cmi = cmi + w(s) * i_ksg(Dx{s}, Dy{s}, o.cmi_k); end
for b = 1:o.cmi_nperm
    v = 0;
    for s = 1:numel(w)
        pp = randperm(rs, size(Dy{s}, 1));
        v  = v + w(s) * i_ksg(Dx{s}, Dy{s}(pp, pp), o.cmi_k);
    end
    nullv(b) = v;
end
idx = 1:numel(y);  if numel(idx) > o.cmi_max_n, idx = randperm(rs, numel(y), o.cmi_max_n); end
mi_pool = i_ksg(i_cheb(Dp(idx, :)), i_cheb(y(idx)), o.cmi_k);
end

function D = i_cheb(A)
D = squareform(pdist(A, 'chebychev'));
D(1:size(D,1)+1:end) = Inf;
end

function I = i_ksg(Dx, Dy, k)
% Kraskov-Stoegbauer-Grassberger estimator (algorithm 1), nats
n  = size(Dx, 1);
e  = mink(max(Dx, Dy), k, 2);  e = e(:, k);
nx = sum(Dx < e, 2);  ny = sum(Dy < e, 2);
I  = psi(k) + psi(n) - mean(psi(nx + 1) + psi(ny + 1));
end

function f = i_folds(G, U, K, rs)
% Fold ids; whole units dealt round-robin within each treatment (else cells)
f = zeros(numel(G), 1);
if isempty(U)
    for gg = unique(G)'
        idx = find(G == gg);  idx = idx(randperm(rs, numel(idx)));
        f(idx) = mod((0:numel(idx)-1)' + randi(rs, K), K) + 1;
    end
else
    [~, ~, ui] = unique(U);
    gU = accumarray(ui, G, [], @(x) x(1));
    fu = zeros(numel(gU), 1);
    for gg = unique(gU)'
        iu = find(gU == gg);  iu = iu(randperm(rs, numel(iu)));
        fu(iu) = mod((0:numel(iu)-1)' + randi(rs, K), K) + 1;
    end
    f = fu(ui);
end
end

function [yhat, kch] = i_nested_pls(X, y, f, fin, kmax)
% Outer CV predictions; ncomp chosen per outer fold by inner-CV PRESS
K = max(f);  yhat = nan(size(y));  kch = nan(K, 1);
for k = 1:K
    te = f == k;  tr = ~te;
    if ~any(te) || nnz(tr) < 10, continue; end
    Xtr = X(tr, :);  ytr = y(tr);  fi = fin{k};
    ninner = arrayfun(@(j) nnz(fi ~= j), 1:max(fi));
    km = max(1, min([kmax, size(X, 2), min(ninner) - 2]));
    press = zeros(1, km);
    for j = 1:max(fi)
        ite = fi == j;
        if ~any(ite), continue; end
        Yh = i_pls_path(Xtr(~ite, :), ytr(~ite), Xtr(ite, :), km);
        press = press + sum((ytr(ite) - Yh).^2, 1);
    end
    [~, kb] = min(press);  kch(k) = kb;
    Yh = i_pls_path(Xtr, ytr, X(te, :), kb);
    yhat(te) = Yh(:, kb);
end
end

function Yh = i_pls_path(Xtr, ytr, Xte, km)
% Predictions for ncomp = 1..km from one SIMPLS fit (BETA_k = W(:,1:k)*YL(:,1:k)')
[~, YL, ~, ~, ~, ~, ~, st] = plsregress(Xtr, ytr, km);
Xc = Xte - mean(Xtr, 1);
Bk = zeros(size(Xtr, 2), 1);
Yh = zeros(size(Xte, 1), km);
for k = 1:km
    Bk = Bk + st.W(:, k) * YL(1, k);
    Yh(:, k) = mean(ytr) + Xc * Bk;
end
end

function q2 = i_q2(y, yhat, f)
% Q2 = 1 - PRESS / TSS, TSS about each fold's TRAINING mean
num = 0;  den = 0;
for k = unique(f(~isnan(yhat)))'
    te = f == k & ~isnan(yhat);  tr = f ~= k;
    num = num + sum((y(te) - yhat(te)).^2);
    den = den + sum((y(te) - mean(y(tr))).^2);
end
q2 = 1 - num / den;
end

function yhat = i_cv_groupmean(G, y, f)
% Treatment-only predictor: training-fold group means
yhat = nan(size(y));
for k = unique(f)'
    te = f == k;  tr = ~te;
    gm = accumarray(G(tr), y(tr), [max(G) 1], @mean, NaN);
    yh = gm(G(te));  yh(isnan(yh)) = mean(y(tr));
    yhat(te) = yh;
end
end

function yhat = i_cv_ols(x, y, f)
yhat = nan(size(y));
for k = unique(f)'
    te = f == k;  tr = ~te;
    b  = [ones(nnz(tr),1) x(tr)] \ y(tr);
    yhat(te) = [ones(nnz(te),1) x(te)] * b;
end
end

function pidx = i_perm_within(G, rs)
pidx = (1:numel(G))';
for gg = unique(G)'
    idx = find(G == gg);
    pidx(idx) = idx(randperm(rs, numel(idx)));
end
end

function vip = i_vip(X, y, k)
% Variable importance in projection (Wold); mean(vip.^2) = 1
[XL, YL, XS, ~, ~, ~, ~, st] = plsregress(X, y, k);
W0  = st.W ./ sqrt(sum(st.W.^2, 1));
ss  = sum(XS.^2, 1) .* sum(YL.^2, 1);
vip = sqrt(size(XL, 1) * sum(ss .* W0.^2, 2) ./ sum(ss, 2));
end

function q = i_q2w(y, yhat, W)
% Within-well Q2: observed and predicted each centred on their well mean
ok = ~isnan(yhat);  y = y(ok);  yh = yhat(ok);  W = W(ok);
[~, ~, wi] = unique(W);  n = accumarray(wi, 1);
my = accumarray(wi, y) ./ n;  mh = accumarray(wi, yh) ./ n;
yc = y - my(wi);  hc = yh - mh(wi);
q  = 1 - sum((yc - hc).^2) / sum(yc.^2);
end

function KP = i_knn_prep(X, f, fin, kmax)
% Neighbour lists (global row ids, sorted by distance) for every outer and inner fold
K = max(f);  KP.out = cell(K, 1);  KP.in = cell(K, 1);
for k = 1:K
    te = find(f == k);  tr = find(f ~= k);
    if isempty(te), continue; end
    KP.out{k} = int32(reshape(tr(knnsearch(X(tr,:), X(te,:), 'K', min(kmax, numel(tr)))), numel(te), []));
    fi = fin{k};  KP.in{k} = cell(max(fi), 1);
    for j = 1:max(fi)
        ite = find(fi == j);  itr = find(fi ~= j);
        if isempty(ite), continue; end
        nb = knnsearch(X(tr(itr),:), X(tr(ite),:), 'K', min(kmax, numel(itr)));
        KP.in{k}{j} = int32(reshape(tr(itr(nb)), numel(ite), []));
    end
end
end

function [yhat, kch] = i_knn_predict(KP, y, f, fin, ks)
% kNN mean prediction; k chosen per outer fold by inner-CV PRESS
K = max(f);  yhat = nan(size(y));  kch = nan(K, 1);
for k = 1:K
    if isempty(KP.out{k}), continue; end
    tr = find(f ~= k);  fi = fin{k};  press = zeros(1, numel(ks));
    for j = 1:numel(KP.in{k})
        I = KP.in{k}{j};  if isempty(I), continue; end
        C   = cumsum(reshape(y(I), size(I)), 2);
        kj  = min(ks, size(I, 2));
        press = press + sum((y(tr(fi == j)) - C(:, kj) ./ kj).^2, 1);
    end
    [~, ib] = min(press);  kb = min(ks(ib), size(KP.out{k}, 2));  kch(k) = kb;
    I = KP.out{k}(:, 1:kb);
    yhat(f == k) = mean(reshape(y(I), size(I)), 2);
end
end

function q = i_q2_levels(y, yh, L, nL)
% Q2 within each level: observed and predicted centred on the level mean
q = nan(1, nL);  ok = ~isnan(yh);
for l = 1:nL
    i = ok & L == l;  if nnz(i) < 3, continue; end
    yc = y(i) - mean(y(i));  hc = yh(i) - mean(yh(i));
    q(l) = 1 - sum((yc - hc).^2) / sum(yc.^2);
end
end

function [q, n, nw, lo, hi] = i_q2_levels_boot(y, yh, L, nL, W, B, rs)
% Within-level Q2 with a 95% CI from resampling wells (cells if no well id)
[q, n, nw, lo, hi] = deal(nan(nL, 1));  ok = ~isnan(yh);
for l = 1:nL
    i = find(ok & L == l);  n(l) = numel(i);
    if n(l) < 3, continue; end
    yc = y(i) - mean(y(i));  hc = yh(i) - mean(yh(i));
    e  = (yc - hc).^2;  tt = yc.^2;
    q(l) = 1 - sum(e) / sum(tt);
    if isempty(W), u = (1:numel(i))'; else, [~, ~, u] = unique(W(i)); nw(l) = max(u); end
    E = accumarray(u, e);  Tt = accumarray(u, tt);  m = numel(E);
    if m < 2, continue; end
    Sb = randi(rs, m, m, B);
    qb = 1 - sum(E(Sb), 1) ./ sum(Tt(Sb), 1);
    lo(l) = prctile(qb, 2.5);  hi(l) = prctile(qb, 97.5);
end
end

%% =========================================================================
%  Step 4b: feature sets
%% =========================================================================
function SR = i_feature_sets(res, Y, D, dnames, dblock, Btab, gk, uk, o, rs, tag)
% ALONE and UNIQUE Q2 for named groups of dynamics columns (same CV design as Step 2)
p  = size(D, 2);  nT = numel(res);
SR = struct('mode', '', 'names', {{}}, 'members', {{}}, 'table', table(), 'vip', [], ...
            'vip_targets', {{}}, 'vip_features', {dnames}, 'n_disc_wells', NaN, 'n_val_wells', NaN, ...
            'reference', o.reference);
evalMask = true(size(D, 1), 1);
fs = o.feature_sets;
if isstruct(fs)
    SR.mode = 'user';
    snames  = fieldnames(fs)';
    sets    = cellfun(@(n) i_resolve_set(fs.(n), dnames, Btab, dblock, tag), snames, 'UniformOutput', false);
else
    switch lower(char(fs))
        case 'blocks'
            SR.mode = 'blocks';
            assert(height(Btab) >= 2, '%s feature_sets ''blocks'' needs >= 2 dynamics blocks.', tag);
            sets   = arrayfun(@(b) find(dblock == b), 1:height(Btab), 'UniformOutput', false);
            snames = cellstr(Btab.label)';
        case 'vip'
            SR.mode = 'vip';
            [disc, evalMask, SR.n_disc_wells, SR.n_val_wells] = i_split_wells(gk, uk, rs);
            % Targets that cleared their Step-2 null drive the grouping (all if none did)
            use = arrayfun(@(r) r.q2_dyn > prctile(r.q2_null, 95), res);
            if ~any(use), use = true(1, nT); end
            iu = find(use);  V = nan(p, numel(iu));
            for j = 1:numel(iu)
                r  = res(iu(j));  rows = r.rows & disc;
                yj = Y(rows, r.tcol);  yj = (yj - mean(yj)) / std(yj);
                V(:, j) = i_vip(D(rows, :), yj, max(1, r.ncomp));
            end
            [sets, snames] = i_vip_groups(V, dnames, o.vip_k);
            SR.vip = V;  SR.vip_targets = {res(iu).target};
        otherwise
            error('%s feature_sets must be ''blocks'', ''vip'' or a struct.', tag);
    end
end
SR.names = snames;  SR.members = cellfun(@(c) dnames(c), sets, 'UniformOutput', false);
nS = numel(sets);

rowsOut = cell(nT, 1);
for t = 1:nT
    t0 = tic;
    r    = res(t);  rows = r.rows & evalMask;
    y    = Y(rows, r.tcol);  y = (y - mean(y)) / std(y);
    Dt   = D(rows, :);
    Gt   = findgroups(gk(rows));
    Ut   = [];  if ~isempty(uk), Ut = findgroups(uk(rows)); end
    f    = i_folds(Gt, Ut, o.nfolds, rs);
    fin  = cell(max(f), 1);
    for k = 1:max(f)
        Ui = [];  if ~isempty(Ut), Ui = Ut(f ~= k); end
        fin{k} = i_folds(Gt(f ~= k), Ui, o.inner_folds, rs);
    end
    mtr  = nan(size(y));                                   % training-fold mean per cell
    for k = 1:max(f), mtr(f == k) = mean(y(f ~= k)); end
    yh_all = i_nested_pls(Dt, y, f, fin, o.max_ncomp);
    q_all  = i_q2(y, yh_all, f);
    W  = Ut;  if isempty(W), W = (1:numel(y))'; end
    tb = cell(nS, 1);
    for s = 1:nS
        c  = sets{s};  cw = setdiff(1:size(Dt, 2), c);
        yh_a = i_nested_pls(Dt(:, c), y, f, fin, o.max_ncomp);
        if isempty(cw), yh_w = mtr; else, yh_w = i_nested_pls(Dt(:, cw), y, f, fin, o.max_ncomp); end
        q_a = i_q2(y, yh_a, f);  q_w = i_q2(y, yh_w, f);  q_u = q_all - q_w;
        [alo, ahi, ulo, uhi] = i_set_boot(y, mtr, yh_all, yh_a, yh_w, W, o.n_boot, rs);
        % Nulls: shuffle the SET's rows within treatment (alone: set only; unique: set inside full model)
        na = nan(o.set_nperm, 1);  nu = nan(o.set_nperm, 1);
        for b = 1:o.set_nperm
            pp = i_perm_within(Gt, rs);
            na(b) = i_q2(y, i_nested_pls(Dt(pp, c), y, f, fin, o.max_ncomp), f);
            Dp = Dt;  Dp(:, c) = Dt(pp, c);
            nu(b) = i_q2(y, i_nested_pls(Dp, y, f, fin, o.max_ncomp), f) - q_w;
        end
        a95 = prctile(na, 95);  u95 = prctile(nu, 95);
        tb{s} = table(string(r.target), string(snames{s}), numel(y), q_all, q_a, alo, ahi, a95, q_a > a95, ...
                      q_u, ulo, uhi, u95, q_u > u95, 'VariableNames', {'target','set','n','q2_all','q2_alone', ...
                      'alone_lo','alone_hi','alone_null95','alone_above','q2_unique','unique_lo','unique_hi', ...
                      'unique_null95','unique_above'});
    end
    rowsOut{t} = vertcat(tb{:});
    if o.vb >= 1, fprintf('  sets %d/%d %s done (%.0f s)\n', t, nT, r.target, toc(t0)); end
end
SR.table = vertcat(rowsOut{:});
end

function c = i_resolve_set(spec, dnames, Btab, dblock, tag)
% Set spec -> column indices: numeric idx, feature names, or block labels
if isnumeric(spec)
    c = unique(spec(:)');
    assert(all(c >= 1 & c <= numel(dnames) & c == round(c)), '%s feature set index out of range.', tag);
    return
end
spec = cellstr(spec);  c = [];
for i = 1:numel(spec)
    hit = find(strcmp(dnames, spec{i}));
    if isempty(hit)
        b = find(strcmpi(cellstr(Btab.label), spec{i}));
        assert(~isempty(b), '%s feature set entry ''%s'' is neither a feature nor a block.', tag, spec{i});
        hit = find(dblock == b);
    end
    c = [c hit(:)']; %#ok<AGROW>
end
c = unique(c);
end

function [disc, val, nd, nv] = i_split_wells(gk, uk, rs)
% Half of the wells of each treatment -> discovery, rest -> validation
if isempty(uk), uk = string((1:numel(gk))'); end
[uu, ~, ui] = unique(uk);  gU = strings(numel(uu), 1);
for i = 1:numel(uu), gU(i) = gk(find(ui == i, 1)); end
isd = false(numel(uu), 1);
for g = unique(gU)'
    iw = find(gU == g);  iw = iw(randperm(rs, numel(iw)));
    isd(iw(1:floor(numel(iw)/2))) = true;
end
disc = isd(ui);  val = ~disc;  nd = nnz(isd);  nv = nnz(~isd);
end

function [sets, names] = i_vip_groups(V, dnames, k)
% Cluster features by VIP profile across targets (Ward); groups ordered by mean VIP
p = size(V, 1);
if p < 4, sets = num2cell(1:p);  names = dnames;  return; end
Zl = linkage(V, 'ward');
if isempty(k)
    ks  = 2:min(6, p - 1);
    sil = arrayfun(@(kk) mean(silhouette(V, cluster(Zl, 'maxclust', kk))), ks);
    [~, ib] = max(sil);  k = ks(ib);
end
cl = cluster(Zl, 'maxclust', k);
mv = accumarray(cl, mean(V, 2), [], @mean);  [~, og] = sort(mv, 'descend');
short = regexprep(dnames, '^latent_', 'L');
sets = cell(1, k);  names = cell(1, k);
for g = 1:k
    m = find(cl == og(g));  [~, oi] = sort(mean(V(m, :), 2), 'descend');  m = m(oi)';
    sets{g}  = m;
    more = '';  if numel(m) > 3, more = sprintf(' +%d', numel(m) - 3); end
    names{g} = sprintf('G%d: %s%s', g, strjoin(short(m(1:min(3, end))), ','), more);
end
end

function [alo, ahi, ulo, uhi] = i_set_boot(y, mtr, yh_all, yh_a, yh_w, W, B, rs)
% Well-bootstrap 95% CIs for alone Q2 and unique Q2 (paired: same wells for all models)
ok = ~isnan(yh_all) & ~isnan(yh_a) & ~isnan(yh_w);
[~, ~, u] = unique(W(ok));
Ef = accumarray(u, (y(ok) - yh_all(ok)).^2);  Ea = accumarray(u, (y(ok) - yh_a(ok)).^2);
Ew = accumarray(u, (y(ok) - yh_w(ok)).^2);    Tt = accumarray(u, (y(ok) - mtr(ok)).^2);
m  = numel(Tt);
if m < 2, [alo, ahi, ulo, uhi] = deal(NaN); return; end
Sb = randi(rs, m, m, B);  st = sum(Tt(Sb), 1);
qa = 1 - sum(Ea(Sb), 1) ./ st;
qu = (sum(Ew(Sb), 1) - sum(Ef(Sb), 1)) ./ st;
alo = prctile(qa, 2.5);  ahi = prctile(qa, 97.5);  ulo = prctile(qu, 2.5);  uhi = prctile(qu, 97.5);
end

function i_sets_brief(SR, tag)
% One line per set
Tb = SR.table;  nT = numel(unique(Tb.target));
note = '';
if strcmp(SR.mode, 'vip'), note = sprintf(' (groups from %d discovery wells, Q2 on %d validation wells)', SR.n_disc_wells, SR.n_val_wells); end
fprintf('%s Step 4b feature sets, mode %s%s\n', tag, SR.mode, note);
fprintf('  %-26s %10s %8s %11s %8s\n', 'set', 'alone Q2', 'above', 'unique Q2', 'above');
for s = 1:numel(SR.names)
    i = Tb.set == string(SR.names{s});
    fprintf('  %-26s %+10.3f %5d/%-2d %+11.3f %5d/%-2d\n', SR.names{s}, mean(Tb.q2_alone(i)), nnz(Tb.alone_above(i)), nT, ...
            mean(Tb.q2_unique(i)), nnz(Tb.unique_above(i)), nT);
end
fprintf('  (mean over targets; above = beats the 95th pct of within-treatment shuffles of that set)\n');
end

function figs = i_sets_fig(SR, o, fig_dir)
% Step 4b figure: overview heatmaps + one bar panel per target
Tb = SR.table;  tg = unique(Tb.target, 'stable');  nT = numel(tg);  nS = numel(SR.names);
cS = lines(max(nS, 2));  if nS <= 6, cS = [0.20 0.45 0.70; 0.80 0.30 0.30; 0.30 0.65 0.35; 0.55 0.35 0.70; 0.90 0.60 0.20; 0.40 0.40 0.40]; cS = cS(1:nS, :); end
hasV = strcmp(SR.mode, 'vip') && ~isempty(SR.vip);
nOv  = 1 + hasV;  ncol = 4;  nrow = ceil((nT + nOv) / ncol);
f  = figure('Color', 'w', 'Name', 'Q3 feature sets', 'Units', 'pixels');
f.Position(3:4) = [1500, 320 * nrow + 90];
tl = tiledlayout(f, nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
sub = 'ALONE = Q^2 from this set only (solid);  UNIQUE = Q^2 lost when the set is removed (light);  bars 95% CI over wells;  black tick = shuffled null 95%';
if strcmp(SR.mode, 'vip')
    sub = {sub, sprintf('Groups defined from VIP on %d discovery wells; all Q^2 measured on the other %d validation wells', SR.n_disc_wells, SR.n_val_wells)};
end
title(tl, [{sprintf('Q3 Step 4b  Which parts of the dynamics carry the prediction?  (%s)', SR.mode)}, cellstr(sub)], ...
      'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'tex');

% Overview: alone Q2, sets x targets
ax = nexttile(tl);
H = nan(nS, nT);  A = false(nS, nT);
for s = 1:nS
    for t = 1:nT
        i = Tb.set == string(SR.names{s}) & Tb.target == tg(t);
        if any(i), H(s, t) = Tb.q2_alone(i);  A(s, t) = Tb.alone_above(i); end
    end
end
cl = max(0.02, max(abs(H(:)), [], 'omitnan'));
imagesc(ax, H, [-cl cl]);  colormap(ax, i_div(256));  cb = colorbar(ax);  cb.Box = 'off';  cb.Label.String = 'alone Q^2';
[rr, cc] = find(A);
if ~isempty(rr), text(ax, cc, rr, '*', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold'); end
set(ax, 'XTick', 1:nT, 'XTickLabel', tg, 'XTickLabelRotation', 45, 'YTick', 1:nS, 'YTickLabel', SR.names, 'FontSize', 7.5);
i_step_title(ax, 'Overview', 'Alone Q2 per set and marker (* = above null)');
i_sty(ax);

% VIP map used to define the groups
if hasV
    ax = nexttile(tl);
    ord = [];  for s = 1:nS, ord = [ord, find(ismember(SR.vip_features, SR.members{s}))]; end %#ok<AGROW>
    Vd  = SR.vip(ord, :);
    imagesc(ax, Vd', [0 max(2, max(Vd(:)))]);  colormap(ax, i_seqmap(256));  cb = colorbar(ax);  cb.Box = 'off';  cb.Label.String = 'VIP (discovery)';
    edges = cumsum(cellfun(@numel, SR.members));
    for e = edges(1:end-1), xline(ax, e + 0.5, 'k-', 'LineWidth', 1.5); end
    set(ax, 'XTick', 1:numel(ord), 'XTickLabel', regexprep(SR.vip_features(ord), '^latent_', 'L'), 'XTickLabelRotation', 90, ...
            'YTick', 1:numel(SR.vip_targets), 'YTickLabel', SR.vip_targets, 'FontSize', 7);
    i_step_title(ax, 'Grouping', 'Features ordered by group (lines = group borders)');
    i_sty(ax);
end

% One panel per target
for t = 1:nT
    ax = nexttile(tl);  hold(ax, 'on');
    Tt = Tb(Tb.target == tg(t), :);  xx = 1:height(Tt);  w = 0.36;
    for s = 1:height(Tt)
        ca = cS(s, :);  if ~Tt.alone_above(s), ca = [0.80 0.80 0.80]; end
        cu = 0.45 * cS(s, :) + 0.55;  if ~Tt.unique_above(s), cu = [0.93 0.93 0.93]; end
        bar(ax, xx(s) - w/2, Tt.q2_alone(s),  w, 'FaceColor', ca, 'EdgeColor', 'none');
        bar(ax, xx(s) + w/2, Tt.q2_unique(s), w, 'FaceColor', cu, 'EdgeColor', cS(s, :) * 0.8);
    end
    errorbar(ax, xx - w/2, Tt.q2_alone,  Tt.q2_alone  - Tt.alone_lo,  Tt.alone_hi  - Tt.q2_alone,  'k', 'LineStyle', 'none', 'CapSize', 3);
    errorbar(ax, xx + w/2, Tt.q2_unique, Tt.q2_unique - Tt.unique_lo, Tt.unique_hi - Tt.q2_unique, 'k', 'LineStyle', 'none', 'CapSize', 3);
    plot(ax, xx - w/2, Tt.alone_null95,  '_', 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 1.3);
    plot(ax, xx + w/2, Tt.unique_null95, '_', 'Color', 'k', 'MarkerSize', 8, 'LineWidth', 1.3);
    yline(ax, Tt.q2_all(1), '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    if ~isempty(SR.reference) && isfield(SR.reference, 'table')
        Rf = SR.reference.table;  ir = string(Rf{:, 1}) == tg(t);
        if any(ir), yline(ax, Rf{ir, 2}, ':', 'Color', [0.8 0.3 0.3], 'LineWidth', 1.3); end
    end
    yline(ax, 0, '-', 'Color', [0.7 0.7 0.7]);
    set(ax, 'XTick', xx, 'XTickLabel', Tt.set, 'XTickLabelRotation', 30, 'FontSize', 7.5);
    xlim(ax, [0.4 height(Tt) + 0.6]);  ylabel(ax, 'Q^2');
    ttl = sprintf('%s   (all features %.3f', tg(t), Tt.q2_all(1));
    if ~isempty(SR.reference) && isfield(SR.reference, 'label'), ttl = [ttl sprintf(', dotted = %s', SR.reference.label)]; end
    title(ax, [ttl ')'], 'Interpreter', 'none', 'FontSize', 8.5);
    i_sty(ax);
end
figs = f;  i_save(f, o, fig_dir, 'Q3-sets', sprintf('%dtg', nT), SR.mode);
end

function pa = i_holm(p)
pa = nan(size(p));  v = ~isnan(p);  pv = p(v);  m = numel(pv);
if m == 0, return; end
[ps, ord] = sort(pv(:));
adj = min(1, cummax((m - (1:m)' + 1) .* ps));
tmp = zeros(m, 1);  tmp(ord) = adj;  pa(v) = tmp;
end


%% =========================================================================
%  Step T: rare producers (sparse bright tail)
%% =========================================================================
function TR = i_tail(res, Y, mk, repk, rep_col, uk, bdv, BL, xf, o, rs, tag)
% Producer calls from a per-replicate gate on the stained distribution; evaluate
% (T1) frequency vs treatment with the well as unit and (T2) within-well ranking of
% producers by the out-of-fold predictions already computed in Step 2.
nm  = {res.target};  nR = numel(res);  nKeep = numel(repk);
if (ischar(o.tail) || isstring(o.tail)) && strcmpi(o.tail, 'all')
    sel = 1:nR;
else
    want = cellstr(o.tail);  [tf, loc] = ismember(i_norm(want), i_norm(nm));
    if any(~tf), warning('%s tail targets not analysed (skipped): %s', tag, strjoin(want(~tf), ', ')); end
    sel = reshape(loc(tf), 1, []);
end
assert(~isempty(sel), '%s no tail target among the analysed targets.', tag);
ks = o.tail_k;  kr = o.tail_kref;  nK = numel(ks);
hasBg = ~isempty(o.tail_bg);
if hasBg
    Tbg = o.tail_bg;
    rb  = repmat("all", height(Tbg), 1);
    if ~isempty(rep_col) && ismember(rep_col, Tbg.Properties.VariableNames), rb = string(Tbg.(rep_col)); end
end

% --- producer calls at kref for EVERY analysed target (co-positivity, artefact check)
L = nan(nKeep, nR);
for t = 1:nR
    r = res(t);  x = Y(r.rows, r.tcol);
    l = i_gate(x, repk(r.rows), kr);
    col = nan(nKeep, 1);  col(r.rows) = l;  L(:, t) = col;
end

% --- T2: identity (within-well AUC) and gate sanity, per tail target and k
rows = {};
for t = sel
    r  = res(t);  x = Y(r.rows, r.tcol);  rr = repk(r.rows);  W = r.W;  G = r.G;
    yl = r.yhat;  yk = r.yhat_nl;  hasK = any(~isnan(yk));
    Pm = zeros(numel(x), 0, 'uint32');
    if ~isempty(W)
        Pm = zeros(numel(x), o.tail_nperm, 'uint32');
        for b = 1:o.tail_nperm, Pm(:, b) = i_perm_within(W, rs); end
    end
    if hasBg
        assert(ismember(mk{r.tcol}, Tbg.Properties.VariableNames), '%s tail_bg lacks column %s.', tag, mk{r.tcol});
        xb = i_transform(double(Tbg.(mk{r.tcol})), xf, o, tag);
    end
    for kk = 1:nK
        [l, g, ur] = i_gate(x, rr, ks(kk));
        fpr = NaN;
        if hasBg
            c = 0;  n = 0;
            for j = 1:numel(ur)
                ib = rb == ur(j) & ~isnan(xb);  c = c + nnz(xb(ib) > g(j));  n = n + nnz(ib);
            end
            fpr = 100 * c / max(n, 1);
        end
        a_tx = i_auc(i_lowo_mean(G, W, double(l)), l);
        [a_l, a_lw, n_lw, p_lw] = i_auc_set(yl, l, W, Pm, ~isempty(W));
        [a_k, a_kw, n_kw, p_kw] = deal(NaN);
        if hasK, [a_k, a_kw, n_kw, p_kw] = i_auc_set(yk, l, W, Pm, ~isempty(W)); end
        rows(end+1, :) = {string(r.target), ks(kk), numel(x), 100 * mean(l), fpr, a_tx, ...
                          a_l, a_lw, n_lw, p_lw, NaN, a_k, a_kw, n_kw, p_kw, NaN}; %#ok<AGROW>
    end
    clear Pm
end
Tb = cell2table(rows, 'VariableNames', {'target','k','n','pct_tail','bg_fpr','auc_tx','auc_lin','auc_lin_w', ...
                'auc_lin_w_null95','p_lin_w','p_lin_w_holm','auc_knn','auc_knn_w','auc_knn_w_null95','p_knn_w','p_knn_w_holm'});
for kk = 1:nK                                                   % Holm across tail targets, per gate
    i = Tb.k == ks(kk);
    Tb.p_lin_w_holm(i) = i_holm(Tb.p_lin_w(i));  Tb.p_knn_w_holm(i) = i_holm(Tb.p_knn_w(i));
end

% --- T1: frequency per level (well = unit), kref
F = {};  FE = {};
if isempty(uk)
    warning('%s no well/unit column: producer frequency per well not computed.', tag);
else
    for t = sel
        r = res(t);  l = L(r.rows, t);  [~, ~, wi] = unique(uk(r.rows));
        fw = 100 * accumarray(wi, l) ./ accumarray(wi, 1);  nw = numel(fw);
        for i = 1:numel(bdv)
            lvAll = BL{i}(r.rows);  lw = strings(nw, 1);
            first = accumarray(wi, (1:numel(wi))', [], @min);  lw(:) = lvAll(first);
            [lu, ~, li] = unique(lw);
            m  = accumarray(li, fw, [], @mean);  s = accumarray(li, fw, [], @std);  n = accumarray(li, 1);
            se = s ./ sqrt(max(n, 1));
            F{end+1} = table(repmat(string(bdv{i}), numel(lu), 1), repmat(string(r.target), numel(lu), 1), lu, n, m, s, ...
                             m - 1.96 * se, m + 1.96 * se, 'VariableNames', {'var','target','level','n_wells','pct_mean', ...
                             'pct_sd','pct_lo','pct_hi'}); %#ok<AGROW>
            e2 = i_eta2(fw, li);  nul = nan(o.tail_nperm, 1);
            for b = 1:o.tail_nperm, nul(b) = i_eta2(fw, li(randperm(rs, nw))); end
            pe = (1 + nnz(nul >= e2)) / (o.tail_nperm + 1);  if isnan(e2), pe = NaN; end   % no variation (e.g. 0 producers)
            FE{end+1} = table(string(bdv{i}), string(r.target), e2, pe, ...
                              'VariableNames', {'var','target','eta2','p'}); %#ok<AGROW>
        end
    end
end
freq = table();  fe = table();
if ~isempty(F), freq = vertcat(F{:});  fe = vertcat(FE{:}); end
Tr = Tb(Tb.k == kr, :);  few = Tr.pct_tail * nKeep / 100 < 20;      % < ~20 producer cells at kref
if any(few)
    warning(['%s gate k = %g calls < 20 cells for %s (%% producers %s). T1/T2 are undefined there. ', ...
             'On a log-like scale (asinh) bright tails are compressed: lower ''tail_k''/''tail_kref'' ', ...
             '(R.tail.table lists %% called per k).'], tag, kr, strjoin(Tr.target(few), ', '), ...
             mat2str(round(Tr.pct_tail(few)', 3)));
end

% --- co-positivity at kref (all targets; enrichment over independence)
E = nan(nR);  N = zeros(nR);
for a = 1:nR
    for b = 1:nR
        ok = ~isnan(L(:, a)) & ~isnan(L(:, b));
        pa = mean(L(ok, a));  pb = mean(L(ok, b));  both = L(ok, a) == 1 & L(ok, b) == 1;
        N(a, b) = nnz(both);  E(a, b) = mean(both) / (pa * pb);
    end
end
isTail = false(1, nR);  isTail(sel) = true;

TR = struct('k', ks, 'kref', kr, 'targets', {nm(sel)}, 'has_bg', hasBg, 'table', Tb, 'freq', freq, ...
            'freq_effect', fe, 'copos', struct('names', {nm}, 'is_tail', isTail, 'enrich', E, 'n_both', N), ...
            'label', L(:, sel), 'group', o.group, 'alpha', o.alpha);
TR.verdict = i_tail_verdict(TR, bdv, o);
end

function [l, g, ur] = i_gate(x, rr, k)
% Producer call: x above median + k x robust SD of its replicate (MAD -> SD; SD if MAD = 0)
ur = unique(rr);  g = nan(numel(ur), 1);  l = false(size(x));
for j = 1:numel(ur)
    i = rr == ur(j) & ~isnan(x);
    if ~any(i), continue; end
    s = 1.4826 * mad(x(i), 1);  if ~(s > 0), s = std(x(i)); end
    g(j) = median(x(i)) + k * s;  l(i) = x(i) > g(j);
end
end

function a = i_auc(s, l)
% Pooled AUC (Mann-Whitney), NaN scores ignored
ok = ~isnan(s);  s = s(ok);  l = logical(l(ok));
n1 = nnz(l);  n0 = nnz(~l);
if n1 == 0 || n0 == 0, a = NaN; return; end
rk = tiedrank(s);
a  = (sum(rk(l)) - n1 * (n1 + 1) / 2) / (n1 * n0);
end

function a = i_auc_strat(s, l, W)
% Within-well AUC: P(a producer outranks a non-producer OF THE SAME WELL), pooled over
% wells by pair count. One sortrows per call (ranks within well; continuous scores).
ok = ~isnan(s);  s = s(ok);  l = double(l(ok));  [~, ~, wi] = unique(W(ok));
n  = numel(s);  [~, ord] = sortrows([wi s]);
pos = zeros(n, 1);  pos(ord) = (1:n)';
cnt = accumarray(wi, 1);  st = [0; cumsum(cnt(1:end-1))];
rk  = pos - st(wi);                                        % rank within well
n1  = accumarray(wi, l);  n0 = cnt - n1;
U   = accumarray(wi, rk .* l) - n1 .* (n1 + 1) / 2;
good = n1 > 0 & n0 > 0;
a = sum(U(good)) / sum(n1(good) .* n0(good));
if ~any(good), a = NaN; end
end

function [a, aw, n95, p] = i_auc_set(s, l, W, Pm, hasW)
% Pooled AUC, within-well AUC, and its within-well shuffle null
a = i_auc(s, l);  [aw, n95, p] = deal(NaN);
if ~hasW || all(isnan(s)), return; end
aw  = i_auc_strat(s, l, W);
nul = nan(size(Pm, 2), 1);
for b = 1:size(Pm, 2), nul(b) = i_auc_strat(s(Pm(:, b)), l, W); end
n95 = prctile(nul, 95);  p = (1 + nnz(nul >= aw)) / (numel(nul) + 1);
end

function m = i_lowo_mean(G, W, l)
% Treatment-only producer score: frequency in the same treatment, leaving the cell's
% own well out (in-sample group mean when no wells)
if isempty(W)
    gs = accumarray(G, l);  gn = accumarray(G, 1);  m = gs(G) ./ gn(G);  return
end
[~, ~, wi] = unique(W);
gs = accumarray(G, l);  gn = accumarray(G, 1);  ws = accumarray(wi, l);  wn = accumarray(wi, 1);
m  = (gs(G) - ws(wi)) ./ max(gn(G) - wn(wi), 1);
m(gn(G) - wn(wi) < 1) = NaN;
end

function e2 = i_eta2(v, g)
% Between-level share of variance (one-way)
mu = mean(v);  mg = accumarray(g, v, [], @mean);  ng = accumarray(g, 1);
e2 = sum(ng .* (mg - mu).^2) / sum((v - mu).^2);
end

function V = i_tail_verdict(TR, bdv, o)
% One line per tail target at kref: gate sanity | frequency | identity
Tb = TR.table(TR.table.k == TR.kref, :);  nt = height(Tb);
[gate, fr, idn, vd] = deal(strings(nt, 1));
v1 = '';  if ~isempty(bdv), v1 = bdv{1}; end
for i = 1:nt
    if TR.has_bg && Tb.bg_fpr(i) > 0.2 * Tb.pct_tail(i), gate(i) = "noisy (bg FPR > 20% of calls)";
    elseif TR.has_bg, gate(i) = "clean";
    else, gate(i) = "no bg control"; end
    fr(i) = "n/a";
    if ~isempty(TR.freq_effect)
        e = TR.freq_effect(TR.freq_effect.var == string(v1) & TR.freq_effect.target == Tb.target(i), :);
        if ~isempty(e)
            if isnan(e.eta2), fr(i) = sprintf("undefined (no producers at k = %g)", TR.kref);
            elseif e.p < o.alpha, fr(i) = sprintf("set by %s (eta2 %.2f)", v1, e.eta2);
            else, fr(i) = sprintf("not set by %s (eta2 %.2f)", v1, e.eta2); end
        end
    end
    lin = Tb.p_lin_w_holm(i) < o.alpha && Tb.auc_lin_w(i) > 0.5;
    knn = Tb.p_knn_w_holm(i) < o.alpha && Tb.auc_knn_w(i) > 0.5;
    if Tb.pct_tail(i) == 0
        idn(i) = "identity undefined (no producers)";
    elseif lin || knn
        idn(i) = sprintf("dynamics pick producers (in-well AUC %.2f)", max([Tb.auc_lin_w(i) Tb.auc_knn_w(i)], [], 'omitnan'));
    else
        idn(i) = "producers not identified by dynamics";
    end
    vd(i) = fr(i) + "; " + idn(i);
    if startsWith(gate(i), "noisy"), vd(i) = "UNRELIABLE GATE: " + vd(i); end
end
V = table(Tb.target, gate, fr, idn, vd, 'VariableNames', {'target','gate','frequency','identity','verdict'});
end

function i_tail_brief(TR, o, tag)
% Console: one row per tail target at kref + frequency + co-positivity lines
Tb = TR.table(TR.table.k == TR.kref, :);
fprintf('\n%s Step T rare producers: gate = median + %g x robust SD of stained cells per replicate (sweep %s in R.tail.table)\n', ...
        tag, TR.kref, mat2str(TR.k));
fprintf('  T1 = does %s set HOW MANY cells produce (wells as units)?  T2 = do dynamics pick WHICH cells produce (same well)?\n', o.group);
fprintf('  %-8s %7s %7s %7s %7s %9s %8s %7s %9s\n', 'target', '% prod', 'bg FPR', 'AUC tx', 'AUC dyn', 'in-well', 'null95', 'p_holm', 'kNN in-w');
for i = 1:height(Tb)
    fprintf('  %-8s %7.2f %7.2f %7.3f %7.3f %9.3f %8.3f %7.3g %9.3f\n', Tb.target(i), Tb.pct_tail(i), Tb.bg_fpr(i), ...
            Tb.auc_tx(i), Tb.auc_lin(i), Tb.auc_lin_w(i), Tb.auc_lin_w_null95(i), Tb.p_lin_w_holm(i), Tb.auc_knn_w(i));
end
if ~isempty(TR.freq)
    v1 = TR.freq.var(1);
    for i = 1:height(Tb)
        Fi = TR.freq(TR.freq.var == v1 & TR.freq.target == Tb.target(i), :);
        if isempty(Fi), continue; end
        [~, hi] = max(Fi.pct_mean);  [~, lo] = min(Fi.pct_mean);
        fprintf('  T1 %-8s by %s: highest %s %.1f%%, lowest %s %.1f%% (x%.1f)\n', Tb.target(i), v1, Fi.level(hi), ...
                Fi.pct_mean(hi), Fi.level(lo), Fi.pct_mean(lo), Fi.pct_mean(hi) / max(Fi.pct_mean(lo), 0.01));
    end
end
C = TR.copos;  it = find(C.is_tail);  nt = find(~C.is_tail);
pr = {};
for a = 1:numel(it), for b = a+1:numel(it), pr{end+1} = sprintf('%s-%s x%.1f', C.names{it(a)}, C.names{it(b)}, C.enrich(it(a), it(b))); end, end %#ok<AGROW>
if ~isempty(pr), fprintf('  co-producers (enrichment over chance): %s\n', strjoin(pr, ', ')); end
if ~isempty(nt)
    [mx, im] = max(C.enrich(it, nt), [], 2, 'omitnan');
    sus = mx > 2;
    if any(sus)
        s1 = string(C.names(it(sus)));  s2 = string(C.names(nt(im(sus))));  mv = mx(sus);
        fprintf('  ! producer calls also enriched in the top tail of non-producer markers (pan-bright cells = artefact suspect): %s\n', ...
                strjoin(reshape(compose('%s~%s x%.1f', s1(:), s2(:), mv(:)), 1, []), ', '));
    end
end
for i = 1:height(TR.verdict), fprintf('  -> %-8s %s\n', TR.verdict.target(i), TR.verdict.verdict(i)); end
end

function figs = i_tail_fig(TR, o, fig_dir)
% Step T figure: gate sanity | T1 frequency | T2 identity | k sensitivity | co-positivity | guide
Tb = TR.table;  tg = string(TR.targets(:))';  nt = numel(tg);  ks = TR.k;  nK = numel(ks);
cR = [0.80 0.30 0.30];  cB = [0.20 0.45 0.70];  cG = [0.55 0.55 0.55];
f  = figure('Color', 'w', 'Name', 'Q3 rare producers', 'Units', 'pixels');
f.Position(3:4) = [1400 780];
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q3 Step T  Rare producers: does %s set HOW MANY cells produce, do dynamics pick WHICH ones?', o.group), ...
      'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

% Gate sanity: % producers per k vs secondary-only false-positive rate
ax = nexttile(tl);  hold(ax, 'on');
cK = [0.75 0.75 0.75; 0.35 0.35 0.35; 0.10 0.10 0.10];  if nK ~= 3, cK = gray(nK + 2); cK = cK(2:end-1, :); end
w = 0.8 / nK;  h = gobjects(0);  lg = {};
for kk = 1:nK
    Tk = Tb(Tb.k == ks(kk), :);  [~, ix] = ismember(tg, Tk.target);
    xx = (1:nt) - 0.4 + w * (kk - 0.5);
    pv = Tk.pct_tail(ix);                                         % log axis: 0 % cannot be drawn
    h(end+1) = bar(ax, xx, max(pv, 0.01), w, 'FaceColor', cK(kk, :), 'EdgeColor', 'none', 'BaseValue', 0.01); %#ok<AGROW>
    z = pv == 0;  if any(z), text(ax, xx(z), 0.013 * ones(1, nnz(z)), '0', 'HorizontalAlignment', 'center', 'FontSize', 7); end
    lg{end+1} = sprintf('k = %g', ks(kk)); %#ok<AGROW>
    if TR.has_bg, plot(ax, xx, Tk.bg_fpr(ix), '_', 'Color', cR, 'MarkerSize', 9, 'LineWidth', 1.6); end
end
if TR.has_bg, h(end+1) = plot(ax, nan, nan, '_', 'Color', cR, 'MarkerSize', 9, 'LineWidth', 1.6); lg{end+1} = 'secondary-only above gate'; end
set(ax, 'XTick', 1:nt, 'XTickLabel', tg, 'YScale', 'log');  ylabel(ax, '% cells above gate');  xlim(ax, [0.5 nt + 0.5]);
yl = [Tb.pct_tail; Tb.bg_fpr];  yl = yl(yl > 0);  if isempty(yl), yl = 1; end
ylim(ax, [0.01, 10^ceil(log10(max(yl) * 1.5))]);
legend(ax, h, lg, 'Box', 'off', 'FontSize', 7, 'Location', 'best');
i_step_title(ax, 'Gate  Are producer calls above background?', 'Bars = stained cells called; red = secondary-only cells passing the same gate');
i_sty(ax);

% T1 frequency by the first breakdown variable (wells as units)
ax = nexttile(tl);  hold(ax, 'on');
if ~isempty(TR.freq)
    v1 = TR.freq.var(1);  Fv = TR.freq(TR.freq.var == v1, :);  lv = unique(Fv.level, 'stable');  nl = numel(lv);
    cT = [cR; cB; 0.30 0.65 0.35; 0.55 0.35 0.70; 0.90 0.60 0.20];  cT = cT(1 + mod(0:nt-1, 5), :);
    wb = 0.8 / nt;
    for t = 1:nt
        Ft = Fv(Fv.target == tg(t), :);  [~, ix] = ismember(lv, Ft.level);  ok = ix > 0;
        xx = (1:nl) - 0.4 + wb * (t - 0.5);
        m = nan(nl, 1);  lo = m;  hi = m;  m(ok) = Ft.pct_mean(ix(ok));  lo(ok) = Ft.pct_lo(ix(ok));  hi(ok) = Ft.pct_hi(ix(ok));
        bar(ax, xx, m, wb, 'FaceColor', cT(t, :), 'EdgeColor', 'none');
        errorbar(ax, xx, m, m - max(lo, 0), hi - m, 'k', 'LineStyle', 'none', 'CapSize', 2);
    end
    set(ax, 'XTick', 1:nl, 'XTickLabel', lv, 'XTickLabelRotation', 35);  xlim(ax, [0.5 nl + 0.5]);
    ylabel(ax, sprintf('%% producers per well (k = %g)', TR.kref));
    legend(ax, cellstr(tg), 'Box', 'off', 'FontSize', 7, 'Location', 'best');
    E = TR.freq_effect(TR.freq_effect.var == v1, :);
    sub = strjoin(reshape(compose('%s eta2 %.2f (p %.2g)', E.target, E.eta2, E.p), 1, []), ' | ');
    i_step_title(ax, sprintf('T1  How many produce, by %s?', v1), sub);
else
    axis(ax, 'off');  i_step_title(ax, 'T1  How many produce?', 'No well column: not computed');
end
i_sty(ax);

% T2 identity at kref
ax = nexttile(tl);  hold(ax, 'on');
Tk = Tb(Tb.k == TR.kref, :);  [~, ix] = ismember(tg, Tk.target);  Tk = Tk(ix, :);  x = 1:nt;
h = gobjects(0);  lg = {};
h(end+1) = plot(ax, x, Tk.auc_tx,  'o', 'MarkerFaceColor', cG, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);   lg{end+1} = 'treatment frequency (pooled)';
h(end+1) = plot(ax, x, Tk.auc_lin, 'o', 'MarkerFaceColor', cR, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);   lg{end+1} = 'dynamics PLS (pooled)';
h(end+1) = plot(ax, x, Tk.auc_lin_w, 'd', 'Color', cR, 'MarkerSize', 9, 'LineWidth', 1.4);                     lg{end+1} = 'dynamics PLS, same well';
if any(~isnan(Tk.auc_knn_w))
    h(end+1) = plot(ax, x, Tk.auc_knn_w, 'd', 'Color', cB, 'MarkerSize', 9, 'LineWidth', 1.4);                 lg{end+1} = 'dynamics kNN, same well';
end
h(end+1) = plot(ax, x, Tk.auc_lin_w_null95, '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);           lg{end+1} = 'same-well shuffle 95%';
sig = Tk.p_lin_w_holm < o.alpha | Tk.p_knn_w_holm < o.alpha;
if any(sig)
    xs = x(sig);  ys = max([Tk.auc_lin_w(sig) Tk.auc_knn_w(sig)], [], 2, 'omitnan');
    text(ax, xs(:), ys(:), '  *', 'FontSize', 11, 'FontWeight', 'bold');
end
yline(ax, 0.5, '-', 'Color', [0.7 0.7 0.7]);
set(ax, 'XTick', x, 'XTickLabel', tg);  xlim(ax, [0.5 nt + 0.5]);  ylabel(ax, 'AUC (producer vs non-producer)');
legend(ax, h, lg, 'Box', 'off', 'FontSize', 7, 'Location', 'best');
if all(isnan(Tk.auc_lin_w)) && all(isnan(Tk.auc_knn_w))
    text(ax, mean(xlim(ax)), 1, sprintf('no producers at k = %g:\nlower ''tail_kref''', TR.kref), ...
         'HorizontalAlignment', 'center', 'Color', cR, 'FontSize', 9);
end
i_step_title(ax, 'T2  Do dynamics pick which cells produce?', 'Same-well AUC removes treatment and well: * = Holm p < alpha');
i_sty(ax);

% k sensitivity
ax = nexttile(tl);  hold(ax, 'on');
for t = 1:nt
    Tt = Tb(Tb.target == tg(t), :);  [~, o2] = sort(Tt.k);  Tt = Tt(o2, :);
    plot(ax, Tt.k, Tt.auc_lin_w, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
end
yline(ax, 0.5, '-', 'Color', [0.7 0.7 0.7]);  xline(ax, TR.kref, ':', 'Color', [0.4 0.4 0.4]);
set(ax, 'XTick', ks);  xlabel(ax, 'gate k (robust SDs above the replicate median)');  ylabel(ax, 'same-well AUC (PLS)');
legend(ax, cellstr(tg), 'Box', 'off', 'FontSize', 7, 'Location', 'best');
i_step_title(ax, 'Sensitivity  Does the answer depend on the gate?', 'Flat lines = conclusion robust to where "producer" starts');
i_sty(ax);

% Co-positivity
ax = nexttile(tl);
C = TR.copos;  Lg = log2(C.enrich);  Lg(1:size(Lg, 1) + 1:end) = NaN;
fin = isfinite(Lg);  cl = 1;  if any(fin(:)), cl = max(1, max(abs(Lg(fin)))); end
Lg(Lg == -Inf) = -cl;  Lg(Lg == Inf) = cl;                      % never / always co-called -> colour limit
% clims cannot precede name-value args in imagesc: set CLim separately
imagesc(ax, Lg, 'AlphaData', ~isnan(Lg));  set(ax, 'CLim', [-cl cl]);  colormap(ax, i_div(256));
cb = colorbar(ax);  cb.Box = 'off';
cb.Label.String = 'log_2 co-call enrichment';
lab = C.names;  lab(C.is_tail) = strcat(lab(C.is_tail), ' *');
set(ax, 'XTick', 1:numel(lab), 'XTickLabel', lab, 'XTickLabelRotation', 45, 'YTick', 1:numel(lab), 'YTickLabel', lab, 'FontSize', 7.5);
i_step_title(ax, 'Co-production  Same cells?', 'Top tails of every marker; * = tail target. Enrichment with non-producer markers = artefact suspect');
i_sty(ax);

% Reading guide
ax = nexttile(tl);
V = TR.verdict;
i_textpanel(ax, [{ ...
  'WHY THIS VIEW', ...
  'Few cells make these markers, brightly; the rest are negative.', ...
  'Q2 on intensity mostly measures the diffuse bulk, which is', ...
  'where the antibody-free technical coupling lives. Here the', ...
  'model is unchanged (continuous); only its EVALUATION asks', ...
  'about the bright tail:', ...
  'Gate  calls must clear the secondary-only false-positive rate', ...
  'T1    treatment may set the PROBABILITY of producing ...', ...
  'T2    ... while dynamics decide WHICH cell does (same well).', ...
  'k     robust conclusions do not move with the gate.', ...
  ''}, cellstr(compose('%s: %s', V.target, V.verdict))']);

figs = f;  i_save(f, o, fig_dir, 'Q3-tail', sprintf('%dtg', nt), o.group);
end


%% =========================================================================
%  Assembly + helpers
%% =========================================================================
function R = i_assemble(res, tnames, ceil_v, can, dnames, dblock, Btab, o)
nT = numel(res);
cmi = [res.cmi]';  mip = [res.mi_pooled]';
nmed = arrayfun(@(r) median(r.cmi_null, 'omitnan'), res)';
n95  = arrayfun(@(r) prctile(r.cmi_null, 95), res)';
pc   = arrayfun(@(r) (1 + nnz(r.cmi_null >= r.cmi)) / (numel(r.cmi_null) + 1), res)';
pc(isnan(cmi)) = NaN;
R.cmi = table(string(tnames(:)), [res.n]', cmi, mip, nmed, n95, cmi - nmed, pc, i_holm(pc), ...
    'VariableNames', {'target','n','cmi','mi_pooled','null_med','null95','excess','p','p_holm'});

q2d  = [res.q2_dyn]';
qmed = arrayfun(@(r) median(r.q2_null), res)';
q95  = arrayfun(@(r) prctile(r.q2_null, 95), res)';
pd   = [res.p_dyn]';  pdl = [res.p_delta]';
S    = [res.single];                                         % [nCan x nT]
if isempty(S), qbest = nan(1, nT); ib = ones(1, nT); else, [qbest, ib] = max(S, [], 1, 'omitnan'); end
bestm = strings(nT, 1);
for t = 1:nT, if ~isnan(qbest(t)), bestm(t) = can{ib(t)}; end, end
pdh = i_holm(pd);  pch = R.cmi.p_holm;
% Permutation resolution: with K targets, Holm p cannot fall below K/(nperm+1)
qres = nT / (o.nperm + 1) < o.alpha;   cres = nT / (o.cmi_nperm + 1) < o.alpha;
R.resolution = struct('q2_floor', nT / (o.nperm + 1), 'cmi_floor', nT / (o.cmi_nperm + 1), ...
                      'q2_resolved', qres, 'cmi_resolved', cres, 'nperm_needed', ceil(nT / o.alpha));
verdict = strings(nT, 1);
for t = 1:nT
    qok = pdh(t) < o.alpha && q2d(t) > 0;
    cok = ~isnan(pch(t)) && pch(t) < o.alpha;
    qun = ~qres && pd(t) <= 1 / (o.nperm + 1) + 1e-12 && q2d(t) > 0;       % beat every permutation
    cun = ~cres && ~isnan(pc(t)) && pc(t) <= 1 / (o.cmi_nperm + 1) + 1e-12;
    if     qok && cok,   verdict(t) = "predictive";
    elseif (qok || qun) && (cok || cun || isnan(pch(t))) && (qun || cun), verdict(t) = "unresolved (raise nperm)";
    elseif ~qok && cok,  verdict(t) = "nonlinear-only";
    elseif qok && isnan(pch(t)), verdict(t) = "predictive (MI not run)";
    elseif qok,          verdict(t) = "linear-only (MI underpowered)";
    else,                verdict(t) = "not-predictive";
    end
end
R.q2 = table(string(tnames(:)), [res.n]', q2d, qmed, q95, pd, pdh, [res.ncomp]', ...
    [res.q2_tx]', [res.q2_txdyn]', [res.dq2]', pdl, i_holm(pdl), ceil_v(:), q2d ./ ceil_v(:), ...
    bestm, qbest(:), verdict, 'VariableNames', {'target','n','q2_dyn','q2_null_med','q2_null95', ...
    'p_dyn','p_dyn_holm','ncomp','q2_tx','q2_txdyn','dq2','p_delta','p_delta_holm','ceiling', ...
    'q2_frac_ceiling','best_marker','q2_best_marker','verdict'});

% VIP with block shares (interpreted only when Q2 established)
nb = height(Btab);
R.vip = struct('target', {}, 'established', {}, 'vip', {}, 'block_share', {});
for t = 1:nT
    v2 = res(t).vip.^2;
    sh = accumarray(dblock(:), v2, [nb 1]) / sum(v2);
    cs = Btab.n_features / sum(Btab.n_features);
    bs = table(Btab.label, sh, cs, sh ./ cs, 'VariableNames', {'block','vip_share','column_share','enrichment'});
    R.vip(t) = struct('target', tnames{t}, 'established', pdh(t) < o.alpha && q2d(t) > 0, ...
                      'vip', res(t).vip, 'block_share', bs);
end
R.pred = struct('target', tnames(:)', 'y', {res.y}, 'yhat', {res.yhat}, 'yhat_nl', {res.yhat_nl}, ...
                'G', {res.G}, 'W', {res.W});
p95 = @(c) cellfun(@(v) prctile(v, 95), c);
R.nl = table(string(tnames(:)), q2d, q95, [res.q2_nl]', p95({res.q2_nl_null})', i_holm([res.p_nl]'), [res.knn_k]', ...
    [res.q2w_lin]', p95({res.q2w_lin_null})', i_holm([res.p_w_lin]'), ...
    [res.q2w_nl]', p95({res.q2w_nl_null})', i_holm([res.p_w_nl]'), ...
    'VariableNames', {'target','q2_lin','q2_lin_null95','q2_knn','q2_knn_null95','p_knn_holm','knn_k', ...
    'q2w_lin','q2w_lin_null95','p_w_lin_holm','q2w_knn','q2w_knn_null95','p_w_knn_holm'});
R.single_q2 = S;  R.single_names = can;
R.targets = tnames;  R.features = dnames;  R.feature_block = dblock;  R.blocks = Btab;
R.cmi_null = {res.cmi_null};  R.q2_null = {res.q2_null};
end

function c = i_resolve_ceiling(cv, tnames, T, keep, xf, o, Xt, tag) %#ok<INUSD>
% Ceiling on explainable variance per target: numeric, or test-retest r from a repeat column
nT = numel(tnames);  c = nan(nT, 1);
if isempty(cv), return; end
if isnumeric(cv)
    assert(numel(cv) == nT, '%s numeric ceiling needs one value per target (%d).', tag, nT);
    c = cv(:);  return
end
assert(isstruct(cv), '%s ceiling must be [], numeric, or struct.', tag);
for t = 1:nT
    fn = matlab.lang.makeValidName(tnames{t});
    if ~isfield(cv, fn), continue; end
    v = cv.(fn);
    if isnumeric(v), c(t) = v; continue; end
    rc = char(v);
    assert(ismember(rc, T.Properties.VariableNames), '%s ceiling repeat column ''%s'' not in T.', tag, rc);
    y1 = Xt(:, t);  y2 = i_transform(double(T.(rc)(keep)), xf, o, tag);
    ok = ~isnan(y1) & ~isnan(y2);
    c(t) = corr(y1(ok), y2(ok));        % test-retest reliability = max explainable fraction
    fprintf('%s ceiling %s: test-retest r = %.3f vs %s (n=%d).\n', tag, tnames{t}, c(t), rc, nnz(ok));
end
end

function [Xt, xlab] = i_transform(X, xf, o, tag)
% NaN-preserving transform (standardization happens per target)
M = size(X, 2);
switch xf
    case 'arcsinh'
        if ischar(o.cofactor) || isstring(o.cofactor)
            c = ones(1, M);
            for j = 1:M
                v = X(X(:,j) > 0, j);
                if ~isempty(v), c(j) = prctile(v, 10); end
            end
            xlab = 'arcsinh(x/c), c=auto p10';
        else
            c = o.cofactor;  if isscalar(c), c = repmat(c, 1, M); end
            xlab = 'arcsinh(x/c), c=user';
        end
        Xt = asinh(X ./ c);
    case 'log2',   Xt = log2(max(X, 1));  xlab = 'log2(max(x,1))';
    case 'linear', Xt = X;  xlab = 'linear';  if o.prenormalized, xlab = 'linear (per-rep z)'; end
    otherwise,     error('%s unknown transform ''%s''.', tag, xf);
end
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
        assert(any(hit), '%s no columns end in ''%s'' (cycif_build_table strips the suffix).', tag, sfx);
        mk = mk(hit);
    end
end
mk_disp = mk;
if ~isempty(sfx), mk_disp = regexprep(mk, [sfx '$'], ''); end
if o.vb >= 2, fprintf('%s %d functional markers: %s\n', tag, numel(mk), strjoin(mk_disp, ', ')); end
end

function k = i_key(Tb, cols)
k = string(Tb.(cols{1}));
for i = 2:numel(cols), k = k + "|" + string(Tb.(cols{i})); end
k = k(:);
end

function v = i_getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
if iscell(d) && (ischar(v) || isstring(v)), v = cellstr(v); end
end

function s = i_norm(c)
s = lower(regexprep(cellstr(c), '[^a-zA-Z0-9]', ''));
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
fprintf('\n%s ===== RESULTS (%d cells, scale: %s) =====\n', tag, R.n_rows, R.transform_label);
fprintf('  Dynamics blocks (block-scaled, each column z-scored then x weight/sqrt(p_b)):\n');
for b = 1:height(R.blocks)
    fprintf('    %-8s %3d features  weight %.2f  NaN policy %s\n', R.blocks.label(b), ...
            R.blocks.n_features(b), R.blocks.weight(b), R.blocks.nan_policy(b));
end
fprintf('  Nulls shuffle the live<->fixed linkage WITHIN treatment: they keep P(dyn|tx) and P(y|tx), so exceeding them means information beyond treatment.\n');
if ~R.resolution.q2_resolved || ~R.resolution.cmi_resolved
    fprintf('  CAUTION: with %d targets, the smallest Holm-adjusted p is %.3g (Q2, nperm=%d) / %.3g (CMI, cmi_nperm=%d), at or above alpha=%.2g.\n', ...
            height(R.q2), R.resolution.q2_floor, o.nperm, R.resolution.cmi_floor, o.cmi_nperm, o.alpha);
    fprintf('           Significance is unreachable at this resolution: use nperm and cmi_nperm >= %d. Targets that beat every permutation are labelled UNRESOLVED, not negative.\n', R.resolution.nperm_needed);
end

for t = 1:height(R.q2)
    q = R.q2(t, :);  c = R.cmi(t, :);
    fprintf('\n  --- %s (n=%d) ---\n', q.target, q.n);
    % (1) existence
    if ~isnan(c.cmi)
        fprintf('  I(dyn; %s | tx) = %.3f nats (null median %.3f, 95th %.3f; p_holm = %.3g). Pooled I(dyn; %s) = %.3f.\n', ...
                q.target, c.cmi, c.null_med, c.null95, c.p_holm, q.target, c.mi_pooled);
        if c.p_holm < o.alpha
            fprintf('  -> CMI exceeds the null: dynamics carry information about %s beyond what treatment dictates; the sensors reflect (or drive) a cell-level decision.\n', q.target);
        elseif ~R.resolution.cmi_resolved && c.p <= 1 / (o.cmi_nperm + 1) + 1e-12
            fprintf('  -> CMI exceeds all %d permutations (excess %.3f nats) but cannot reach Holm significance at this resolution: UNRESOLVED.\n', o.cmi_nperm, c.excess);
        else
            fprintf('  -> CMI ~ null: conditional on treatment, dynamics are EPIPHENOMENAL for %s.\n', q.target);
        end
        if c.mi_pooled > 0
            rr = max(c.excess, 0) / c.mi_pooled;
            if rr <= 1
                fprintf('     %.0f%% of the pooled dependence survives conditioning on treatment (the rest is treatment-mediated).\n', 100 * rr);
            else
                fprintf('     Conditional dependence exceeds the pooled estimate (x%.1f): pooling across treatments partly masks within-treatment structure (KSG estimates at different n are only roughly comparable).\n', rr);
            end
        end
    end
    % (2) characterization
    fprintf('  Q2Y: dyn %+.3f (within-tx null median %+.3f, 95th %+.3f; p_holm = %.3g; ncomp = %d)\n', ...
            q.q2_dyn, q.q2_null_med, q.q2_null95, q.p_dyn_holm, q.ncomp);
    fprintf('       tx-only %+.3f | tx+dyn %+.3f | dQ2 = %+.3f (p_holm = %.3g)', q.q2_tx, q.q2_txdyn, q.dq2, q.p_delta_holm);
    if strlength(q.best_marker) > 0
        fprintf(' | best single marker %s %+.3f', q.best_marker, q.q2_best_marker);
    end
    fprintf('\n');
    if ~isnan(q.ceiling)
        fprintf('  Ceiling (test-retest) = %.3f: dynamics reach %.0f%% of the explainable variance.', q.ceiling, 100*q.q2_frac_ceiling);
        if q.q2_frac_ceiling >= 0.8, fprintf(' -> near-optimal prediction of this function.\n');
        else, fprintf(' -> substantial explainable variance remains unpredicted.\n'); end
    else
        fprintf('  No ceiling supplied: Q2 is benchmarked against 0 only; measurement noise caps it at an unknown value < 1.\n');
    end
    nl = R.nl(t, :);
    if ~isnan(nl.q2w_lin)
        fprintf('  Within-well Q2 (well means removed): PLS %+.3f (null 95th %+.3f, p_holm %.3g)', nl.q2w_lin, nl.q2w_lin_null95, nl.p_w_lin_holm);
        if ~isnan(nl.q2w_knn), fprintf(' | kNN %+.3f (null 95th %+.3f, p_holm %.3g)', nl.q2w_knn, nl.q2w_knn_null95, nl.p_w_knn_holm); end
        fprintf('\n');
        if nl.q2w_lin > nl.q2w_lin_null95 || (~isnan(nl.q2w_knn) && nl.q2w_knn > nl.q2w_knn_null95)
            fprintf('  -> Dynamics rank cells WITHIN a well: which cells express more %s is not set by the well environment alone.\n', q.target);
        else
            fprintf('  -> No within-well signal: the pooled Q2 for %s reflects well-to-well differences the dynamics share with the endpoint.\n', q.target);
        end
    end
    if ~isnan(nl.q2_knn)
        fprintf('  kNN (k=%d) pooled Q2 %+.3f vs PLS %+.3f (null 95th %+.3f, p_holm %.3g).', nl.knn_k, nl.q2_knn, nl.q2_lin, nl.q2_knn_null95, nl.p_knn_holm);
        gap = nl.q2_knn - nl.q2_lin;
        if gap > 0.02
            fprintf(' -> kNN gains %+.3f: part of the dependence is NONLINEAR in latent space (e.g. a subpopulation, not a gradient).\n', gap);
        elseif gap < -0.02
            fprintf(' -> PLS ahead by %.3f: the dependence is a smooth linear gradient; kNN pays a variance cost.\n', -gap);
        else
            fprintf(' -> comparable: no evidence of nonlinear structure beyond the linear projection.\n');
        end
    end
    switch char(q.verdict)
        case 'predictive'
            fprintf('  -> PREDICTIVE: dependence exists (CMI) and is linearly recoverable (Q2 beyond the within-tx null).\n');
        case 'nonlinear-only'
            fprintf('  -> NONLINEAR: CMI > 0 but PLSR Q2 is not established; the linear projection misses the dependence (try the latent block or a nonlinear predictor).\n');
        case 'unresolved (raise nperm)'
            fprintf('  -> UNRESOLVED: dynamics beat every within-treatment permutation, but %d permutations cannot certify it after Holm. Rerun with nperm >= %d.\n', o.nperm, R.resolution.nperm_needed);
        case {'linear-only (MI underpowered)', 'predictive (MI not run)'}
            fprintf('  -> LINEAR signal without CMI support: the kNN MI on %d PCs/%d cells per stratum may be underpowered; treat as provisional.\n', o.cmi_dims, o.cmi_max_n);
        otherwise
            fprintf('  -> NOT PREDICTIVE: dynamics add nothing detectable to treatment for %s.\n', q.target);
    end
    if q.dq2 > 0 && q.p_delta_holm < o.alpha
        fprintf('  -> Adding dynamics to treatment improves prediction (dQ2 %+.3f): state is not reducible to the stimulus for %s.\n', q.dq2, q.target);
    end
    if strlength(q.best_marker) > 0 && ~isnan(q.q2_best_marker)
        if q.q2_dyn > q.q2_best_marker
            fprintf('  -> Dynamics out-predict every single canonical marker (best %s): the live view is the better state-definer for %s.\n', q.best_marker, q.target);
        else
            fprintf('  -> %s alone predicts %s better than dynamics: the endpoint co-variation is tighter than the live link (note: marker models use co-stained cells only).\n', q.best_marker, q.target);
        end
    end
    % VIP (gated)
    v = R.vip(t);
    if v.established
        bs = v.block_share;
        fprintf('  VIP block shares (share of sum VIP^2 vs share of columns):\n');
        for b = 1:height(bs)
            fprintf('    %-8s %.2f vs %.2f (x%.2f)\n', bs.block(b), bs.vip_share(b), bs.column_share(b), bs.enrichment(b));
        end
        [~, ord] = sort(v.vip, 'descend');  ord = ord(1:min(5, numel(ord)));
        fprintf('    top features: %s\n', strjoin(compose('%s (%.2f)', string(R.features(ord))', v.vip(ord)), ', '));
        fl = find(bs.block == "flag");
        if ~isempty(fl)
            if bs.enrichment(fl) > 1
                fprintf('  -> pulse_flag carries more weight than an average column (x%.2f): pulsatility the encoder may smooth away IS function-relevant for %s.\n', bs.enrichment(fl), q.target);
            else
                fprintf('  -> pulse_flag carries <= average weight (x%.2f): the encoder''s treatment of pulses is not costing predictive information for %s.\n', bs.enrichment(fl), q.target);
            end
        end
    else
        fprintf('  VIP not interpreted: Q2Y not established for %s (Q2-before-VIP).\n', q.target);
    end
end

% Across targets: degeneracy signature
pred = R.q2.target(R.q2.verdict == "predictive" | startsWith(R.q2.verdict, "linear") | R.q2.verdict == "predictive (MI not run)");
nonl = R.q2.target(R.q2.verdict == "nonlinear-only");
none = R.q2.target(R.q2.verdict == "not-predictive");
unr  = R.q2.target(R.q2.verdict == "unresolved (raise nperm)");
fprintf('\n  SUMMARY  predictable: %s | nonlinear-only: %s | unresolved: %s | not predictable: %s\n', ...
        i_list(pred), i_list(nonl), i_list(unr), i_list(none));
if ~isempty(unr)
    fprintf('  -> %d targets beat every permutation but are unresolved at nperm=%d: no conclusion until rerun.\n', numel(unr), o.nperm);
end
if ~isempty(pred) && ~isempty(none)
    fprintf('  -> Dynamics predict some functions and not others: "state" is not one monolithic variable (degeneracy signature).\n\n');
elseif isempty(pred) && isempty(nonl) && isempty(unr)
    fprintf('  -> No function is predictable from dynamics beyond treatment.\n\n');
else
    fprintf('\n');
end
end

function i_brief(R, o, tag)
% Compact console summary: one row per target, then Step 3 and verdicts
Q = R.q2;  C = R.cmi;  N = R.nl;
hasW = any(~isnan(N.q2w_lin));  hasK = any(~isnan(N.q2_knn));  hasC = any(~isnan(C.cmi));
fprintf('\n%s Step 1 = any dependence beyond %s (CMI) | Step 2 = linear prediction on held-out wells (PLS Q2)\n', tag, o.group);
fprintf('  %-10s', 'target');
if hasC, fprintf(' %7s %7s |', 'CMI', 'p_holm'); end
fprintf(' %7s %7s %7s %7s', 'Q2 dyn', 'null95', 'Q2 tx', 'dQ2');
if hasW, fprintf(' %8s', 'in-well'); end
if hasK, fprintf(' %7s', 'Q2 kNN'); end
fprintf('  verdict\n');
for t = 1:height(Q)
    fprintf('  %-10s', Q.target(t));
    if hasC, fprintf(' %7.3f %7.3g |', C.cmi(t), C.p_holm(t)); end
    fprintf(' %+7.3f %+7.3f %+7.3f %+7.3f', Q.q2_dyn(t), Q.q2_null95(t), Q.q2_tx(t), Q.dq2(t));
    if hasW, fprintf(' %+8.3f', N.q2w_lin(t)); end
    if hasK, fprintf(' %+7.3f', N.q2_knn(t)); end
    fprintf('  %s\n', Q.verdict(t));
end
% Step 3: where (summarised over targets that passed Step 2)
est = Q.target([R.vip.established]);
for i = 1:numel(R.breakdown)
    Tb = R.breakdown(i).table;
    if isempty(est) || isempty(Tb), continue; end
    Tb = Tb(ismember(Tb.target, est), :);
    [lv, ~, li] = unique(Tb.level);
    mq = accumarray(li, Tb.q2, [], @(x) mean(x, 'omitnan'));
    na = accumarray(li, Tb.above_null);
    [~, ib] = max(mq);  [~, iw] = min(mq);
    fprintf('  Step 3 by %-9s: highest %s (mean Q2 %.3f), lowest %s (%.3f); %d/%d levels above null for every predictable target\n', ...
            R.breakdown(i).var, lv(ib), mq(ib), lv(iw), mq(iw), nnz(na == numel(est)), numel(lv));
end
pred = Q.target(startsWith(Q.verdict, "predictive") | startsWith(Q.verdict, "linear"));
fprintf('  -> predictable beyond %s: %s | nonlinear-only: %s | not: %s\n', o.group, i_list(pred), ...
        i_list(Q.target(Q.verdict == "nonlinear-only")), i_list(Q.target(Q.verdict == "not-predictive")));
unr = Q.target(Q.verdict == "unresolved (raise nperm)");
if ~isempty(unr), fprintf('  ! unresolved at nperm = %d: %s (use nperm >= %d)\n', o.nperm, i_list(unr), R.resolution.nperm_needed); end
end

%% =========================================================================
%  Step-by-step summary figure + breakdown bar figures
%% =========================================================================
function figs = i_story(R, o, fig_dir)
Q = R.q2;  C = R.cmi;  N = R.nl;  nT = height(Q);
[~, ord] = sort(Q.q2_dyn, 'descend');  x = 1:nT;  tl_ = cellstr(Q.target(ord));
cR = [0.80 0.30 0.30];  cB = [0.20 0.45 0.70];  cG = [0.55 0.55 0.55];  cP = [0.50 0.30 0.70];
f  = figure('Color', 'w', 'Name', 'Q3 dynamics predict function', 'Units', 'pixels');
f.Position(3:4) = [1400 780];
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q3  Do live dynamics predict function beyond %s?   (%d cells, %d targets)', o.group, R.n_rows, nT), ...
      'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

% Step 1: CMI
ax = nexttile(tl);  hold(ax, 'on');
if any(~isnan(C.cmi))
    hn = plot(ax, x, C.null95(ord), '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);
    hc = plot(ax, x, C.cmi(ord), 'o', 'MarkerFaceColor', cR, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    sig = C.p_holm(ord) < o.alpha;
    xs = x(sig);  if any(sig), text(ax, xs(:), C.cmi(ord(sig)), '  *', 'FontSize', 11, 'FontWeight', 'bold'); end
    legend(ax, [hc hn], {'observed', 'shuffled within condition, 95%'}, 'Box', 'off', 'FontSize', 7.5, 'Location', 'best');
    ylabel(ax, 'I(dynamics; marker | condition), nats');
else
    text(ax, 0.5, 0.5, 'CMI not run', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
set(ax, 'XTick', x, 'XTickLabel', tl_, 'XTickLabelRotation', 45);  xlim(ax, [0.5 nT+0.5]);
i_step_title(ax, 'Step 1  Any dependence?', 'Model-free test, conditioned on treatment (* = Holm p < alpha)');
i_sty(ax);

% Step 2: Q2 horse race
ax = nexttile(tl);  hold(ax, 'on');
h = gobjects(0);  lg = {};
h(end+1) = plot(ax, x, Q.q2_tx(ord),    'o', 'MarkerFaceColor', cG, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);  lg{end+1} = 'treatment only';
h(end+1) = plot(ax, x, Q.q2_dyn(ord),   'o', 'MarkerFaceColor', cR, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);  lg{end+1} = 'dynamics (PLS)';
h(end+1) = plot(ax, x, Q.q2_txdyn(ord), 'o', 'MarkerFaceColor', cP, 'MarkerEdgeColor', 'none', 'MarkerSize', 6);  lg{end+1} = 'treatment + dynamics';
if any(~isnan(N.q2_knn))
    h(end+1) = plot(ax, x, N.q2_knn(ord), 'o', 'Color', cB, 'MarkerSize', 8, 'LineWidth', 1.3);  lg{end+1} = 'dynamics (kNN)';
end
if any(~isnan(N.q2w_lin))
    h(end+1) = plot(ax, x, N.q2w_lin(ord), 'd', 'Color', cR, 'MarkerSize', 7, 'LineWidth', 1.1);  lg{end+1} = 'dynamics, within well';
end
h(end+1) = plot(ax, x, Q.q2_null95(ord), '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);  lg{end+1} = 'dynamics null 95%';
yline(ax, 0, '-', 'Color', [0.7 0.7 0.7]);
set(ax, 'XTick', x, 'XTickLabel', tl_, 'XTickLabelRotation', 45);  xlim(ax, [0.5 nT+0.5]);
legend(ax, h, lg, 'Box', 'off', 'FontSize', 7, 'Location', 'northeast');
ylabel(ax, 'Q^2 (variance predicted in held-out wells)');
i_step_title(ax, 'Step 2  How much is predictable?', 'Held-out wells; dynamics vs treatment vs both');
i_sty(ax);

% Reading guide
ax = nexttile(tl);
pred = Q.target(startsWith(Q.verdict, "predictive") | startsWith(Q.verdict, "linear"));
i_textpanel(ax, { ...
  'WHY THIS ORDER', ...
  '1  Existence first, with no model assumed: do dynamics share', ...
  '   any information with the marker among cells that got the', ...
  '   same treatment? (CMI vs within-treatment shuffles)', ...
  '2  Then size, with a linear model tested on wells it never saw:', ...
  '   how much marker variance do dynamics predict, and is it', ...
  '   more than treatment alone? 1 yes + 2 no = nonlinear link.', ...
  '3  Then where: is the prediction carried by some conditions /', ...
  '   replicates and absent in others? (bars, next figure)', ...
  '4  Only for targets that pass 2: which dynamic features carry', ...
  '   the prediction (VIP). VIP of a failed model means nothing.', ...
  '', ...
  sprintf('PREDICTABLE  %s', i_list(pred)), ...
  sprintf('NOT          %s', i_list(Q.target(Q.verdict == "not-predictive"))), ...
  sprintf('NONLINEAR    %s', i_list(Q.target(Q.verdict == "nonlinear-only")))});

% Step 3: breakdown heatmaps (first two variables)
for i = 1:2
    ax = nexttile(tl);
    if i > numel(R.breakdown) || isempty(R.breakdown(i).table)
        axis(ax, 'off');  continue
    end
    Tb = R.breakdown(i).table;
    [lv, ~, li] = unique(Tb.level);
    [~, ti] = ismember(Tb.target, Q.target(ord));
    H = nan(numel(lv), nT);  A = false(numel(lv), nT);
    H(sub2ind(size(H), li, ti)) = Tb.q2;  A(sub2ind(size(A), li, ti)) = Tb.above_null;
    cl = max(0.02, max(abs(H(:)), [], 'omitnan'));
    imagesc(ax, H, [-cl cl]);  colormap(ax, i_div(256));  cb = colorbar(ax);  cb.Box = 'off';  cb.Label.String = 'Q^2 within level';
    [rr, cc] = find(A);
    if ~isempty(rr), text(ax, cc, rr, '*', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold'); end
    set(ax, 'XTick', 1:nT, 'XTickLabel', tl_, 'XTickLabelRotation', 45, 'YTick', 1:numel(lv), 'YTickLabel', lv, 'FontSize', 7.5);
    i_step_title(ax, sprintf('Step 3  Where? (by %s)', R.breakdown(i).var), 'Q2 inside each level (* = above shuffled null)');
    i_sty(ax);
end

% Step 4: VIP for established targets
ax = nexttile(tl);
est = find([R.vip.established]);
if ~isempty(est)
    V = cell2mat(arrayfun(@(k) R.vip(k).vip(:)', est(:), 'UniformOutput', false));
    fn = R.features;
    if numel(fn) > 30, [~, top] = maxk(mean(V, 1), 30);  top = sort(top);  V = V(:, top);  fn = fn(top); end
    imagesc(ax, V, [0 max(2, max(V(:)))]);  colormap(ax, i_seqmap(256));  cb = colorbar(ax);  cb.Box = 'off';  cb.Label.String = 'VIP';
    set(ax, 'XTick', 1:numel(fn), 'XTickLabel', fn, 'XTickLabelRotation', 90, 'YTick', 1:numel(est), ...
            'YTickLabel', {R.vip(est).target}, 'FontSize', 7);
    i_step_title(ax, 'Step 4  Which features?', 'VIP > 1 = above-average contribution (passed Step 2 only)');
else
    axis(ax, 'off');
    i_step_title(ax, 'Step 4  Which features?', 'No target passed Step 2: VIP not interpreted');
end
i_sty(ax);

figs = f;  i_save(f, o, fig_dir, 'Q3-summary', sprintf('%dtg', nT), o.group);
end

function figs = i_breakdown_figs(R, o, fig_dir)
% One figure per breakdown variable: bars of within-level Q2 per target
figs = gobjects(0);
cB = [0.20 0.45 0.70];  cG = [0.78 0.78 0.78];
for i = 1:numel(R.breakdown)
    Tb = R.breakdown(i).table;  if isempty(Tb), continue; end
    tg = unique(Tb.target, 'stable');  nT = numel(tg);
    f  = figure('Color', 'w', 'Name', ['Q3 breakdown ' R.breakdown(i).var], 'Units', 'pixels');
    nl = numel(unique(Tb.level));
    f.Position(3:4) = [min(1600, 260 * min(nT, 4) + 25 * nl * min(nT, 4)), 300 * ceil(nT / 4) + 80];
    tl = tiledlayout(f, ceil(nT / 4), min(nT, 4), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, {sprintf('Q3 Step 3  Where are markers predictable from dynamics?  Breakdown by %s', R.breakdown(i).var), ...
               'Q^2 inside each level (held-out wells); bars 95% CI over wells; black tick = shuffled null 95%; blue = above null'}, ...
          'FontWeight', 'bold', 'FontSize', 10);
    for t = 1:nT
        ax = nexttile(tl);  hold(ax, 'on');
        Tt = Tb(Tb.target == tg(t), :);  xx = 1:height(Tt);
        b  = bar(ax, xx, Tt.q2, 0.7, 'FaceColor', 'flat', 'EdgeColor', 'none');
        b.CData = repmat(cG, height(Tt), 1);  b.CData(Tt.above_null, :) = repmat(cB, nnz(Tt.above_null), 1);
        errorbar(ax, xx, Tt.q2, Tt.q2 - Tt.ci_lo, Tt.ci_hi - Tt.q2, 'k', 'LineStyle', 'none', 'CapSize', 3);
        plot(ax, xx, Tt.null95, '_', 'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5);
        yline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
        set(ax, 'XTick', xx, 'XTickLabel', Tt.level, 'XTickLabelRotation', 45, 'FontSize', 7.5);
        xlim(ax, [0.4 height(Tt) + 0.6]);
        ylabel(ax, 'Q^2 within level');
        title(ax, char(tg(t)), 'Interpreter', 'none', 'FontSize', 9);
        i_sty(ax);
    end
    figs(end+1) = f; %#ok<AGROW>
    i_save(f, o, fig_dir, 'Q3-breakdown', sprintf('%dtg', nT), R.breakdown(i).var);
end
end

function cm = i_div(n)
% Blue - white - red
h  = floor(n/2);
cm = [linspace(0.13,1,h)' linspace(0.40,1,h)' linspace(0.67,1,h)'; ...
      linspace(1,0.70,n-h)' linspace(1,0.09,n-h)' linspace(1,0.17,n-h)'];
end

function cm = i_seqmap(n)
% White -> blue
cm = [linspace(1, 0.12, n)' linspace(1, 0.33, n)' linspace(1, 0.60, n)'];
end

function s = i_list(c)
c = cellstr(c);
if isempty(c), s = '{}'; else, s = ['{' strjoin(c, ', ') '}']; end
end


%% =========================================================================
%  Plots
%% =========================================================================
function figs = i_plots(R, o, rs, fig_dir)
nT = height(R.q2);  figs = gobjects(0);
blk = cellstr(R.blocks.label);
ctag = sprintf('%dtg', nT);  facet = sprintf('%s_%s', strjoin(blk, '-'), o.group);
bcol = lines(numel(blk));

% 1 --- CMI vs within-treatment null --------------------------------------
if o.cmi && any(~isnan(R.cmi.cmi))
    f  = figure('Color', 'w', 'Name', 'dynpred cmi');
    tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact');
    for t = 1:nT
        ax = nexttile(tl);
        i_hist(ax, R.cmi_null{t}, 30, [0.6 0.6 0.6]);
        xline(ax, R.cmi.cmi(t), 'r', 'LineWidth', 2);
        title(ax, sprintf('%s  p_holm=%.2g', R.cmi.target(t), R.cmi.p_holm(t)), 'Interpreter', 'none');
        xlabel(ax, 'I(dyn; y | tx), nats');
    end
    figs(end+1) = f;  i_save(f, o, fig_dir, 'dynpred-cmi', ctag, facet);
end

% 2 --- Q2 vs null and ceiling -------------------------------------------
f  = figure('Color', 'w', 'Name', 'dynpred q2');
ax = axes(f);  hold(ax, 'on');
bar(ax, R.q2.q2_dyn, 'FaceColor', [0.10 0.50 0.90]);
plot(ax, 1:nT, R.q2.q2_null95, 'kd', 'MarkerFaceColor', 'k');
if any(~isnan(R.q2.ceiling)), plot(ax, 1:nT, R.q2.ceiling, 'rv', 'MarkerFaceColor', 'r'); end
yline(ax, 0, '-k');
set(ax, 'XTick', 1:nT, 'XTickLabel', R.q2.target, 'TickLabelInterpreter', 'none', 'XTickLabelRotation', 45);
lg = {'Q2 dynamics', 'within-tx null 95%'};  if any(~isnan(R.q2.ceiling)), lg{end+1} = 'ceiling'; end
legend(ax, lg, 'Location', 'best');  ylabel(ax, 'Q^2Y (outer CV)');
title(ax, 'Predictability of each function from dynamics');
figs(end+1) = f;  i_save(f, o, fig_dir, 'dynpred-q2', ctag, facet);

% 3 --- VIP (established targets only), coloured by block -----------------
est = find([R.vip.established]);
if ~isempty(est)
    f  = figure('Color', 'w', 'Name', 'dynpred vip');
    tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact');
    for t = est
        ax = nexttile(tl);  hold(ax, 'on');
        for b = 1:numel(blk)
            ib = find(R.feature_block == b);
            if isempty(ib), continue; end
            bar(ax, ib, R.vip(t).vip(ib), 'FaceColor', bcol(b,:), 'EdgeColor', 'none', 'BarWidth', 1);
        end
        yline(ax, 1, '--k');
        nf = numel(R.features);
        if nf <= 40
            set(ax, 'XTick', 1:nf, 'XTickLabel', R.features, 'TickLabelInterpreter', 'none', 'XTickLabelRotation', 90);
        end
        legend(ax, blk, 'Location', 'best', 'Interpreter', 'none');
        title(ax, sprintf('VIP: %s (Q2=%.2f)', R.q2.target(t), R.q2.q2_dyn(t)), 'Interpreter', 'none');
    end
    figs(end+1) = f;  i_save(f, o, fig_dir, 'dynpred-vip', ctag, facet);

    % 4 --- Predicted vs actual (outer CV) ---------------------------------
    f  = figure('Color', 'w', 'Name', 'dynpred predicted');
    tl = tiledlayout(f, 'flow', 'TileSpacing', 'compact');
    for t = est
        ax = nexttile(tl);  hold(ax, 'on');
        y = R.pred(t).y;  yh = R.pred(t).yhat;  ok = find(~isnan(yh));
        if numel(ok) > o.max_plot_pts, ok = ok(randperm(rs, numel(ok), o.max_plot_pts)); end
        scatter(ax, y(ok), yh(ok), 5, R.pred(t).G(ok), 'filled', 'MarkerFaceAlpha', 0.3);
        lim = [min([y(ok); yh(ok)]) max([y(ok); yh(ok)])];
        plot(ax, lim, lim, '--k');
        xlabel(ax, 'observed (z)');  ylabel(ax, 'predicted, outer CV (z)');
        title(ax, sprintf('%s  Q2=%.2f', R.q2.target(t), R.q2.q2_dyn(t)), 'Interpreter', 'none');
    end
    figs(end+1) = f;  i_save(f, o, fig_dir, 'dynpred-predobs', ctag, facet);
end

% 5 --- Horse race ---------------------------------------------------------
f  = figure('Color', 'w', 'Name', 'dynpred horserace');
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl);
bar(ax, [R.q2.q2_tx R.q2.q2_dyn R.q2.q2_txdyn R.q2.q2_best_marker]);
yline(ax, 0, '-k');
set(ax, 'XTickLabel', R.q2.target, 'TickLabelInterpreter', 'none', 'XTickLabelRotation', 45);
legend(ax, {'treatment', 'dynamics', 'tx + dyn', 'best single marker'}, 'Location', 'best');
ylabel(ax, 'Q^2Y');  title(ax, 'Which view predicts each function?');
ax = nexttile(tl);
imagesc(ax, R.single_q2, 'AlphaData', ~isnan(R.single_q2));  colorbar(ax);
set(ax, 'XTick', 1:nT, 'XTickLabel', R.q2.target, 'YTick', 1:numel(R.single_names), ...
        'YTickLabel', R.single_names, 'TickLabelInterpreter', 'none', 'XTickLabelRotation', 45);
xlabel(ax, 'target');  ylabel(ax, 'single-marker predictor');  title(ax, 'Single-marker Q^2Y');
figs(end+1) = f;  i_save(f, o, fig_dir, 'dynpred-horserace', ctag, facet);
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

function i_hist(ax, x, nb, col, normpdf)
% Histogram via histcounts + bar (robust to a shadowed histogram.m on the path)
x = x(isfinite(x));  if isempty(x), return; end
[c, e] = histcounts(x, nb);
if nargin > 4 && normpdf, c = c / (sum(c) * mean(diff(e))); end
bar(ax, e(1:end-1) + diff(e)/2, c, 1, 'FaceColor', col, 'EdgeColor', 'none');
end