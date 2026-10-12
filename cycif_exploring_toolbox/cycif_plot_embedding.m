function [fig, emb, emb_cell_idx, embed_features] = cycif_plot_embedding(T, varargin)
% cycif_plot_embedding  Dimensionality reduction of CyCIF markers.
%                       Supports UMAP, tSNE, and PCA. Embedding is computed
%                       once on the full (subsampled) dataset; facet panels
%                       share the same coordinate space.
%
% INPUT
%   T  — table from cycif_build_table
%
% NAME-VALUE OPTIONS
%   'markers'       {1xM} cell  marker names; [] = auto-detect
%   'group_by'      string      T column to color by          (default 'tx1_label')
%   'facet_by'      string      T column to split panels      (default '' = single)
%   'method'        string      'tsne'|'pca'|'umap'           (default 'tsne')
%   'colormap'      string      'bright'|'parula'|any MATLAB map (default 'bright')
%   'use_log2'      logical     log2(max(x,1)) before embed   (default true)
%   'zscore_embed'  logical     per-marker z-score before embed (default true)
%                               set FALSE for pre-scaled input (e.g. a pooled
%                               table already per-replicate z-scored by
%                               cycif_pool_replicates) to avoid double-scaling
%   'exclude_spike' logical     drop ptx_is_spike wells       (default false)
%                               NOTE: ptx_is_spike flags Ptx-at-timepoint
%                               wells. For designs where Ptx is a core
%                               condition (e.g. glucose starve), these are NOT
%                               disposable controls -- keep them (default).
%   'max_pts'       scalar      subsample cap (Inf = none)    (default Inf)
%   'n_neighbors'   scalar      UMAP n_neighbors              (default 15)
%   'min_dist'      scalar      UMAP min_dist                 (default 0.1)
%   'alpha'         scalar      scatter alpha                 (default 0.4)
%   'mksz'          scalar      scatter marker size           (default 8)
%   'max_nan_frac'  scalar      drop markers > frac NaN       (default 1.0 = keep all)
%   'impute_nan'    string      NaN fill strategy for embedding:
%                               'lowerquantile' — fill with per-marker global 25th pct
%                               'median'        — per-condition median (original)
%                               'none'          — no imputation (rows with NaN dropped)
%                               (default 'none')
%   'perplexity'    scalar      tSNE perplexity               (default 30)
%   'seed'          scalar      RNG seed for reproducibility  (default 0)
%   'fig_pos'       [x y w h]   figure position               (default auto)
%   'save_path'     string      full path to save SVG; '' = no save
%   'umap_handler'  UmapHandler optional handler instance to reuse across calls
%                               (avoids re-init); [] constructs one on first
%                               UMAP use. Only consulted when method='umap'.
%
% OUTPUT
%   fig          — figure handle
%   emb          — [nCells x 2] embedding coordinates
%                  rows correspond to T after spike filtering and subsampling
%   emb_cell_idx — table [nCells x 2] with xy and cellid columns
%                  use this to join emb rows back to dataloc or original T
%   embed_features — [nCells x nFeat] the prep feature matrix the embedding was
%                  built on (arcsinh/PCA/PC-dropped). Row-aligned to emb and
%                  emb_cell_idx. Pass this to cycif_cluster('features',...) so
%                  clustering and embedding share identical cells + features.
%
% USAGE
%   cycif_plot_embedding(T);
%   [fig, emb_tsne, idx] = cycif_plot_embedding(T, 'method','tsne', 'seed',0);
%   cycif_plot_embedding(T, 'group_by','tx1_label', 'facet_by','ptx_pool_label');
%   % UMAP (reuse one handler across many calls):
%   uh = UmapHandler;
%   cycif_plot_embedding(T, 'method','umap', 'umap_handler',uh, 'seed',0);

p = inputParser();
addParameter(p, 'markers',       []);
addParameter(p, 'group_by',      'tx1_label');
addParameter(p, 'facet_by',      '');
addParameter(p, 'method',        'tsne');
addParameter(p, 'colormap',      'bright');
addParameter(p, 'use_log2',      true);
addParameter(p, 'zscore_embed',  true);   % per-marker z-score before embedding;
                                          % set false when input is pre-scaled
                                          % (e.g. pooled per-replicate z-score)
addParameter(p, 'exclude_spike', false);
addParameter(p, 'drop_ctrl',     false);
addParameter(p, 'max_pts',       Inf);
addParameter(p, 'n_neighbors',   15);
addParameter(p, 'min_dist',      0.1);
addParameter(p, 'alpha',         0.4);
addParameter(p, 'mksz',          8);
addParameter(p, 'max_nan_frac',  1.0);
addParameter(p, 'min_marker_frac', 0.5);   % per-cell completeness floor (fraction of markers measured)
addParameter(p, 'keep_cells',    []);      % optional external cell keep-list (e.g. from QC.survived_good), overrides min_marker_frac
addParameter(p, 'impute_nan',    'none');  % 'none'|'lowerquantile'|'median'
% --- feature-prep options (delegated to cycif_prep_features) ---
addParameter(p, 'transform',     '');      % ''|'arcsinh'|'log2'|'linear'; '' => map from use_log2
addParameter(p, 'cofactor',      0.4);     % arcsinh cofactor (data-driven default)
addParameter(p, 'n_pcs',         7);       % PCA components for embed features; 0 = off
addParameter(p, 'drop_pcs',      [1]);     % PCs to drop (amplitude axis); [] = none
addParameter(p, 'normalize_cells','none'); % per-cell norm (invalid on signed data; leave none)
addParameter(p, 'perplexity',    30);
addParameter(p, 'seed',          0);
addParameter(p, 'fig_pos',       []);
addParameter(p, 'save_path',     '');
addParameter(p, 'umap_handler',  []);   % optional UmapHandler instance to reuse;
                                        % [] = construct one on first UMAP call
parse(p, varargin{:});

markers        = p.Results.markers;
group_by       = p.Results.group_by;
facet_by       = p.Results.facet_by;
method_req     = lower(p.Results.method);
palette        = p.Results.colormap;
use_log2       = p.Results.use_log2;
zscore_embed   = logical(p.Results.zscore_embed);
exclude_spike  = p.Results.exclude_spike;
drop_ctrl      = p.Results.drop_ctrl;
max_pts        = p.Results.max_pts;
n_neighbors    = p.Results.n_neighbors;
min_dist       = p.Results.min_dist;
umap_handler   = p.Results.umap_handler;
alpha          = p.Results.alpha;
mksz           = p.Results.mksz;
max_nan_frac   = p.Results.max_nan_frac;
min_marker_frac = p.Results.min_marker_frac;
keep_cells     = p.Results.keep_cells;
impute_nan     = lower(p.Results.impute_nan);
% --- feature-prep params (delegated to cycif_prep_features) ---
transform      = lower(string(p.Results.transform));
cofactor       = p.Results.cofactor;
n_pcs          = p.Results.n_pcs;
drop_pcs       = p.Results.drop_pcs;
normalize_cells= lower(string(p.Results.normalize_cells));
% Back-compat: if 'transform' not given, derive it from legacy use_log2
% (true->log2, false->linear). Explicit 'transform' always wins.
if transform == ""
    if use_log2, transform = "log2"; else, transform = "linear"; end
end
perplexity     = p.Results.perplexity;
seed           = p.Results.seed;
fig_pos        = p.Results.fig_pos;

if isempty(markers), markers = i_auto_markers(T); end

% --- Drop markers exceeding NaN threshold --------------------------------
nan_frac = cellfun(@(m) mean(isnan(T.(m))), markers);
keep     = nan_frac <= max_nan_frac;
if any(~keep)
    fprintf('[cycif_plot_embedding] Dropping %d markers with >%.0f%% NaN:\n', ...
            sum(~keep), 100*max_nan_frac);
    fprintf('  %s\n', markers{~keep});
end
markers = markers(keep);

% --- Always drop degenerate columns (independent of max_nan_frac) --------
% All-NaN columns (e.g. unacquired control stains) and zero-variance columns
% carry no information AND break downstream steps: an all-NaN column makes the
% impute_nan='none' row-drop (all(~isnan,2)) reject EVERY cell, collapsing the
% feature matrix and the PCA. max_nan_frac's default (1.0) does not catch
% 100%-NaN columns, so filter them explicitly here.
allnan = cellfun(@(m) all(isnan(T.(m))), markers);
zerov  = cellfun(@(m) local_nzvar(T.(m)) == 0, markers);
degen  = allnan | zerov;
if any(degen)
    fprintf('[cycif_plot_embedding] Dropping %d degenerate markers (all-NaN or zero-variance):\n', ...
            sum(degen));
    fprintf('  %s\n', strjoin(markers(degen), ', '));
end
markers = markers(~degen);
assert(~isempty(markers), 'No markers survive NaN/variance filter — check the linked table.');
assert(numel(markers) >= 2, ...
    ['[cycif_plot_embedding] Only %d usable marker column(s) after filtering; ' ...
     'need >=2 for a 2D embedding.'], numel(markers));
nMk = numel(markers);

% --- Filter cells (spike / ctrl) via shared helper ----------------------
T = cycif_filter_cells(T, group_by, exclude_spike, drop_ctrl);

% --- Subsample BEFORE prep (so T and features stay row-aligned) ----------
% Subsampling first means cycif_prep_features operates on the reduced T, and
% its keep_row mask indexes directly back into this same T -- preserving the
% emb_cell_idx / faceting alignment the downstream plotting relies on.
if height(T) > max_pts
    rng(seed);   % reproducible subsample
    T = T(randperm(height(T), max_pts), :);
    fprintf('[cycif_plot_embedding] Subsampled to %d cells (before prep).\n', max_pts);
end

% --- Build embedding features via the shared prep pipeline ---------------
% ENFORCED single feature pipeline (cycif_prep_features): transform -> cell
% completeness filter -> impute -> z-score -> PCA(-drop). This replaces the
% former inline log2/impute/zscore path, which embedded IMPUTED cells and
% collapsed low-completeness rows to identical points (degenerate k-NN graph,
% UMAP spectral-init failure, tSNE blob). Prep drops those cells and PCA-
% decorrelates, fixing the degeneracy at the source. All methods (umap/tsne/
% pca) now embed on the same X that clustering uses.
[Xz, prep_idx, prep_info] = cycif_prep_features(T, ...
    'markers',         markers, ...
    'transform',       char(transform), ...
    'cofactor',        cofactor, ...
    'normalize_cells', char(normalize_cells), ...
    'zscore',          zscore_embed, ...
    'n_pcs',           n_pcs, ...
    'drop_pcs',        drop_pcs, ...
    'impute_nan',      impute_nan, ...
    'min_marker_frac', min_marker_frac, ...
    'keep_cells',      keep_cells, ...
    'verbose',         true);

% Align T to the cells prep kept, so emb_cell_idx + faceting match Xz rows.
T = T(prep_info.keep_row, :);
assert(height(T) == size(Xz,1), ...
    '[cycif_plot_embedding] Row misalignment after prep (%d T vs %d Xz).', height(T), size(Xz,1));
fprintf('[cycif_plot_embedding] Embedding %d cells x %d features (%s%s, n_pcs=%d).\n', ...
    size(Xz,1), size(Xz,2), prep_info.transform, ...
    i_dropstr(prep_info), prep_info.n_pcs);

% --- Set RNG seed --------------------------------------------------------
rng(seed);

% --- Embed ---------------------------------------------------------------
switch method_req
    case 'umap'
        % Reuse a passed UmapHandler, else construct one (validates umap on
        % first use via its own ready() check). No silent PCA fallback: a UMAP
        % failure now errors loudly rather than masquerading as PCA output.
        if isempty(umap_handler)
            umap_handler = UmapHandler;
        end
        assert(isa(umap_handler, 'UmapHandler'), ...
            '[cycif_plot_embedding] ''umap_handler'' must be a UmapHandler instance.');
        fprintf('[cycif_plot_embedding] Running UMAP (n_neighbors=%d, min_dist=%.3g, seed=%d)...\n', ...
                n_neighbors, min_dist, seed);
        emb    = umap_handler.embed(Xz, 'n_neighbors', n_neighbors, ...
                                    'min_dist', min_dist, 'n_components', 2, 'seed', seed);
        method = 'UMAP';

    case 'tsne'
        perp = min(perplexity, floor(height(T)/4));
        assert(perp >= 1, ...
            ['[cycif_plot_embedding] Too few cells (%d) for tSNE after NaN filtering.\n' ...
             'Use ''impute_nan'',''lowerquantile'' to retain all cells.'], height(T));
        if perp < perplexity
            warning('[cycif_plot_embedding] Perplexity capped at %d (n/4) for n=%d cells.', perp, height(T));
        end
        fprintf('[cycif_plot_embedding] Running tSNE (perplexity=%d, seed=%d)...\n', perp, seed);
        emb    = tsne(Xz, 'NumDimensions', 2, 'Perplexity', perp, 'Standardize', false);
        method = 'tSNE';

    case 'pca'
        [~, emb] = pca(Xz);
        emb      = emb(:, 1:2);
        method   = 'PCA';

    otherwise
        error('[cycif_plot_embedding] Unknown method ''%s''. Use ''umap'', ''tsne'', or ''pca''.', method_req);
end
fprintf('[cycif_plot_embedding] %s embedding complete.\n', method);

% --- Cell index for joining emb back to the source table ----------------
% Always carry (xy, cellid). If the table is a POOLED multi-replicate table
% (has 'rep'/'uid' from cycif_pool_replicates), carry those too: (xy,cellid)
% alone is NOT unique across replicates, so a downstream join on (xy,cellid)
% would mis-map cells between reps. cycif_cluster prefers 'uid' when present.
join_cols = {'xy','cellid'};
if ismember('uid', T.Properties.VariableNames); join_cols{end+1} = 'uid'; end
if ismember('rep', T.Properties.VariableNames); join_cols{end+1} = 'rep'; end
emb_cell_idx = T(:, join_cols);

% Expose the exact feature matrix the embedding was built on, so downstream
% clustering (cycif_cluster 'features') uses the SAME cells in the SAME order
% as emb -- eliminating the alignment drift from building X twice.
embed_features = Xz;

% --- Group / facet levels ------------------------------------------------
grp_levels = unique(T.(group_by), 'stable');
nGrp       = numel(grp_levels);
cmap       = i_colormap(nGrp, palette);

use_facet = ~isempty(facet_by);
if use_facet
    facet_levels = unique(T.(facet_by), 'stable');
else
    facet_levels = {''};
end
nFacet = numel(facet_levels);

% --- Figure layout -------------------------------------------------------
legend_frac = 0.15;
panel_frac  = (1 - legend_frac) / nFacet;
panel_h     = 0.72;
panel_bot   = 0.12;

if isempty(fig_pos)
    fig_pos = [80 80 min(480*nFacet + 200, 2400) 500];
end

fig = figure('Color','w', 'Position', fig_pos, ...
             'Name', sprintf('%s | group=%s facet=%s', method, group_by, facet_by));

% --- Shared axis limits (square) across all facets -----------------------
x_pad    = 0.05 * range(emb(:,1));
y_pad    = 0.05 * range(emb(:,2));
xlims    = [min(emb(:,1))-x_pad,  max(emb(:,1))+x_pad];
ylims    = [min(emb(:,2))-y_pad,  max(emb(:,2))+y_pad];
ax_range = max(diff(xlims), diff(ylims));
xc = mean(xlims); yc = mean(ylims);
xlims = [xc - ax_range/2, xc + ax_range/2];
ylims = [yc - ax_range/2, yc + ax_range/2];

for fi = 1:nFacet
    if use_facet
        sel         = T.(facet_by) == facet_levels(fi);
        panel_title = string(facet_levels(fi));
    else
        sel         = true(height(T), 1);
        panel_title = method;
    end

    left = (fi-1) * panel_frac;
    ax   = axes(fig, 'Position', [left, panel_bot, panel_frac*0.88, panel_h]); %#ok<LAXES>
    hold(ax, 'on');

    emb_sel = emb(sel, :);
    grp_sel = T.(group_by)(sel);

    for gi = nGrp:-1:1
        s = grp_sel == grp_levels(gi);
        scatter(ax, emb_sel(s,1), emb_sel(s,2), mksz, cmap(gi,:), 'filled', ...
                'MarkerFaceAlpha', alpha, 'DisplayName', char(grp_levels(gi)));
    end

    title(ax, panel_title, 'Interpreter','none', 'FontSize', 10);
    xlabel(ax, [method ' 1'], 'Interpreter','none');
    if fi == 1, ylabel(ax, [method ' 2'], 'Interpreter','none'); end
    xlim(ax, xlims); ylim(ax, ylims);
    grid(ax, 'off'); box(ax, 'off');
    set(ax, 'XTick',[], 'YTick',[]);
    ax.FontSize        = 9;
    ax.DataAspectRatio = [1 1 1];
    hold(ax, 'off');
end

% --- Dedicated legend axis -----------------------------------------------
lh = gobjects(nGrp, 1);
for gi = 1:nGrp
    lh(gi) = patch(NaN, NaN, cmap(gi,:), 'FaceAlpha', 0.8, ...
                   'EdgeColor', cmap(gi,:), ...
                   'DisplayName', char(grp_levels(gi)));
end
ax_leg = axes(fig, 'Position',[1-legend_frac, 0, legend_frac, 1], 'Visible','off');
leg = legend(ax_leg, lh, 'Interpreter','none', 'FontSize', 8);
leg.Box      = 'off';
leg.Units    = 'normalized';
leg.Position = [1-legend_frac+0.01, 0.5-leg.Position(4)/2, leg.Position(3), leg.Position(4)];

% --- Sgtitle -------------------------------------------------------------
xform_str = char(prep_info.transform);   % actual transform used (from prep), not use_log2
if isfield(prep_info,'dropped_pcs') && ~isempty(prep_info.dropped_pcs)
    xform_str = sprintf('%s|dropPC%s', xform_str, strjoin(string(prep_info.dropped_pcs), ','));
end
spike_str  = 'spikes excluded'; if ~exclude_spike, spike_str = 'spikes included'; end
sgtitle(fig, sprintf('%s of CyCIF markers  |  group: %s  |  facet: %s  |  %s  |  %s', ...
        method, group_by, facet_by, xform_str, spike_str), ...
        'FontWeight','bold', 'FontSize', 10, 'Interpreter','none');

% --- Save as SVG ---------------------------------------------------------
if ~isempty(p.Results.save_path)
    fpath = p.Results.save_path;
    if ~strcmpi(fpath(end-3:end), '.svg'), fpath = [fpath '.svg']; end
    fig.Units             = 'inches';
    fig_sz                = fig.Position(3:4);
    fig.PaperUnits        = 'inches';
    fig.PaperSize         = fig_sz;
    fig.PaperPosition     = [0, 0, fig_sz];
    fig.PaperPositionMode = 'manual';
    print(fig, fpath, '-dsvg', '-painters');
    fprintf('Saved SVG: %s\n', fpath);
end
end


% =========================================================================
%  PRIVATE HELPERS
% =========================================================================

function markers = i_auto_markers(T)
    cols    = T.Properties.VariableNames;
    is_num  = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
    is_id   = ismember(cols, {'xy','cellid','ptx_is_spike'});
    markers = cols(is_num & ~is_id);
end

function v = local_nzvar(col)
% Variance of the non-NaN entries; 0 if fewer than 2 valid values.
c = col(~isnan(col));
if numel(c) < 2; v = 0; else; v = var(c); end
end

function s = i_dropstr(prep_info)
% Compact ' drop PCk' tag for the log line; '' when nothing dropped.
if isfield(prep_info,'dropped_pcs') && ~isempty(prep_info.dropped_pcs)
    s = sprintf(' drop PC%s', strjoin(string(prep_info.dropped_pcs), ','));
else
    s = '';
end
end

function keep = local_resolve_keeplist(keep_cells, T)
% Resolve an external cell keep-list into a logical mask over the rows of T.
% Accepts: logical [height(T)x1]; numeric row indices into T; or a table with
% 'xy' and 'cellid' columns (matched against T.xy / T.cellid). The table form
% is the natural bridge from QC output, e.g.:
%   keep_tbl = QC_cell_table(QC.survived_good >= thr, {'xy','cellid'})
n = height(T);
if islogical(keep_cells)
    assert(numel(keep_cells)==n, ...
        'keep_cells logical length (%d) must equal height(T) (%d).', numel(keep_cells), n);
    keep = keep_cells(:);
elseif isnumeric(keep_cells)
    keep = false(n,1);
    valid = keep_cells(keep_cells>=1 & keep_cells<=n);
    keep(round(valid)) = true;
elseif istable(keep_cells)
    assert(all(ismember({'xy','cellid'}, keep_cells.Properties.VariableNames)), ...
        'keep_cells table must have ''xy'' and ''cellid'' columns.');
    key_T    = string(T.xy) + "_" + string(T.cellid);
    key_keep = string(keep_cells.xy) + "_" + string(keep_cells.cellid);
    keep = ismember(key_T, key_keep);
else
    error('keep_cells must be logical, numeric indices, or an xy+cellid table.');
end
end

function cmap = i_colormap(n, palette)
if nargin < 2 || isempty(palette), palette = 'bright'; end

switch lower(palette)
    case 'bright'
        base = [ ...
            0.900, 0.100, 0.100;
            0.100, 0.500, 0.900;
            0.100, 0.850, 0.300;
            0.000, 0.620, 0.200;
            0.000, 0.390, 0.100;
            0.900, 0.100, 0.550;
            1.000, 0.700, 0.000;
            0.960, 0.480, 0.000;
            0.980, 0.860, 0.100];
        if n <= size(base,1), cmap = base(1:n,:);
        else, cmap = lines(n); warning('i_colormap: n=%d exceeds palette, using lines().', n); end

    case 'pal1'
        base = [ ...
            0.259, 0.263, 0.282;
            1.000, 0.435, 0.761;
            0.647, 0.000, 0.380;
            0.463, 0.337, 0.875;
            0.106, 0.106, 0.106;
            0.776, 0.855, 0.992];
        if n <= size(base,1), cmap = base(1:n,:);
        else, cmap = lines(n); warning('i_colormap: n=%d exceeds palette, using lines().', n); end

    case 'pal2'
        base = [ ...
            0.945, 0.431, 0.745;
            0.898, 0.702, 0.847;
            0.059, 0.420, 0.996;
            0.875, 0.918, 0.996;
            0.710, 0.620, 0.800;
            0.118, 0.122, 0.098;
            0.647, 0.647, 0.647];
        if n <= size(base,1), cmap = base(1:n,:);
        else, cmap = lines(n); warning('i_colormap: n=%d exceeds palette, using lines().', n); end

    case 'parula'
        raw = parula(max(n,2));
        idx = round(linspace(1, size(raw,1), n));
        cmap = raw(idx,:);

    case 'lines'
        cmap = lines(n);

    otherwise
        try
            raw  = feval(lower(palette), max(n,2));
            idx  = round(linspace(1, size(raw,1), n));
            cmap = raw(idx,:);
        catch
            warning('i_colormap: unknown palette ''%s'', using lines().', palette);
            cmap = lines(n);
        end
end
end