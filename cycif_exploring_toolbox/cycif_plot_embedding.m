function fig = cycif_plot_embedding(T, markers, group_by, facet_by, varargin)
% cycif_plot_embedding  Dimensionality reduction of CyCIF markers.
%                       Supports UMAP, tSNE, and PCA. Embedding is computed
%                       once on the full (subsampled) dataset; facet panels
%                       share the same coordinate space.
%
% INPUT
%   T         — table from cycif_build_table
%   markers   — {1×M} marker names (columns in T); [] = auto-detect
%   group_by  — T column to color embedding by      (default 'tx1_label')
%   facet_by  — T column to split into subplots     (default '' = single panel)
%
% NAME-VALUE OPTIONS
%   'method'       string   'umap' | 'tsne' | 'pca'          (default 'umap')
%   'use_log2'     logical  log2(max(x,1)) before embedding   (default true)
%   'exclude_spike' logical drop ptx_is_spike wells          (default true)
%   'max_pts'      scalar   random subsample cap             (default Inf = no subsample)
%   'n_neighbors'  scalar   UMAP n_neighbors                 (default 15)
%   'min_dist'     scalar   UMAP min_dist                    (default 0.1)
%   'alpha'        scalar   scatter marker alpha             (default 0.4)
%   'mksz'         scalar   scatter marker size              (default 8)
%   'max_nan_frac' scalar   drop markers with > this frac NaN (default 1.0 = keep all)
%   'fig_pos'      [x y w h]                                  (default auto)
%
% OUTPUT
%   fig  — figure handle
%
% USAGE
%   fig = cycif_plot_embedding(T, [], 'tx1_label', 'ptx_pool_label');
%   fig = cycif_plot_embedding(T, [], 'tx1_label', '', 'use_log2', false);
%   fig = cycif_plot_embedding(T, [], 'tx1_label', '', 'exclude_spike', false);
%   fig = cycif_plot_embedding(T, [], 'tx1_label', '', 'method', 'tsne');
%   fig = cycif_plot_embedding(T, [], 'tx1_label', '', 'method', 'pca');

if nargin < 3 || isempty(group_by), group_by = 'tx1_label'; end
if nargin < 4,                       facet_by = '';          end

% --- Parse name-value options ---------------------------------------------
p = inputParser();
addParameter(p, 'method',        'tsne');
addParameter(p, 'use_log2',      true);
addParameter(p, 'exclude_spike', true);
addParameter(p, 'max_pts',       Inf);
addParameter(p, 'n_neighbors',   15);
addParameter(p, 'min_dist',      0.1);
addParameter(p, 'alpha',         0.4);
addParameter(p, 'mksz',          8);
addParameter(p, 'max_nan_frac',  1.0);
addParameter(p, 'fig_pos',       []);
parse(p, varargin{:});

method_req     = lower(p.Results.method);
use_log2       = p.Results.use_log2;
exclude_spike  = p.Results.exclude_spike;
max_pts        = p.Results.max_pts;
n_neighbors    = p.Results.n_neighbors;
min_dist       = p.Results.min_dist;
alpha          = p.Results.alpha;
mksz           = p.Results.mksz;
max_nan_frac   = p.Results.max_nan_frac;
fig_pos        = p.Results.fig_pos;

if isempty(markers), markers = i_auto_markers(T); end

% --- Drop markers exceeding NaN threshold ----------------------------
% Heterogeneous panels mean condition-specific antibodies (e.g. iNOS,
% CD206) are NaN in wells that didn't receive that stain. Requiring
% all markers to be present collapses the embedding to a biased subset.
nan_frac = cellfun(@(m) mean(isnan(T.(m))), markers);
keep     = nan_frac <= max_nan_frac;
if any(~keep)
    fprintf('[cycif_plot_embedding] Dropping %d markers with >%.0f%% NaN:\n', ...
            sum(~keep), 100*max_nan_frac);
    fprintf('  %s\n', markers{~keep});
end
markers = markers(keep);
assert(~isempty(markers), 'No markers survive NaN filter — lower max_nan_frac.');
nMk = numel(markers);

% --- Filter spike wells ---------------------------------------------------
if exclude_spike && ismember('ptx_is_spike', T.Properties.VariableNames)
    n_spike = sum(T.ptx_is_spike);
    if n_spike > 0
        warning('[cycif_plot_embedding] Excluding %d spike-treatment cells (ptx_is_spike=true). Pass ''exclude_spike'',false to include.', n_spike);
    end
    T = T(~T.ptx_is_spike, :);
end

% --- Build marker matrix ------------------------------------------------
Xmat = table2array(T(:, markers));
if use_log2, Xmat = log2(max(Xmat, 1)); end

% --- Per-condition median imputation ------------------------------------
% ASSUMPTION: NaN intensities for condition-specific antibodies (e.g. iNOS
% in M2 wells, CD206 in M1 wells) are filled with the median of cells in
% the same tx1_label group that WERE stained. This preserves all cells in
% the embedding but is circular: condition identity is used to impute values
% that partly drive condition separation. Treat the embedding as exploratory
% only — do not use for classification or statistical testing.
cond_labels  = T.tx1_label;
grp_imp      = unique(cond_labels, 'stable');
n_imputed    = 0;
for m = 1:nMk
    col = Xmat(:, m);
    if ~any(isnan(col)), continue; end
    for g = 1:numel(grp_imp)
        in_grp  = cond_labels == grp_imp(g);
        is_nan  = isnan(col);
        fill_idx = in_grp & is_nan;
        if ~any(fill_idx), continue; end
        grp_med = median(col(in_grp & ~is_nan), 'omitnan');
        if isnan(grp_med)
            % Entire group unstained — fall back to global median
            grp_med = median(col(~is_nan), 'omitnan');
        end
        if isnan(grp_med), grp_med = 0; end   % all-NaN marker
        col(fill_idx)  = grp_med;
        n_imputed      = n_imputed + sum(fill_idx);
    end
    Xmat(:, m) = col;
end
fprintf('[cycif_plot_embedding] ASSUMPTION — per-condition median imputation:\n');
fprintf('  %d cell-marker values imputed across %d markers.\n', n_imputed, nMk);
fprintf('  Condition label (tx1_label) was used to compute imputed medians.\n');
fprintf('  Embedding is exploratory only — do not use for classification.\n');

% --- Subsample -----------------------------------------------------------
if height(T) > max_pts
    idx  = randperm(height(T), max_pts);
    T    = T(idx, :);
    Xmat = Xmat(idx, :);
    fprintf('[cycif_plot_embedding] Subsampled to %d cells.\n', max_pts);
end

% --- Z-score per marker --------------------------------------------------
mu   = mean(Xmat);
sig  = std(Xmat);
sig(sig == 0) = 1;   % guard against zero-variance markers
Xz   = (Xmat - mu) ./ sig;

% --- Embed ---------------------------------------------------------------
switch method_req
    case 'umap'
        try
            emb    = run_umap(Xz, 'n_neighbors', n_neighbors, 'min_dist', min_dist, ...
                              'verbose', 'none');
            method = 'UMAP';
        catch ME
            warning('[cycif_plot_embedding] run_umap failed — falling back to PCA.\n  Error: %s', ME.message);
            fprintf('[cycif_plot_embedding] Error identifier: %s\n', ME.identifier);
            [~, emb] = pca(Xz); emb = emb(:,1:2);
            method = 'PCA (umap fallback)';
        end

    case 'tsne'
        % tsne is built into MATLAB R2017b+
        % Perplexity auto-set to min(30, n/4) to avoid warnings on small n
        perp = min(30, floor(height(T)/4));
        fprintf('[cycif_plot_embedding] Running tSNE (perplexity=%d)...\n', perp);
        emb    = tsne(Xz, 'NumDimensions', 2, 'Perplexity', perp, ...
                      'Standardize', false);   % already z-scored above
        method = 'tSNE';

    case 'pca'
        [~, emb] = pca(Xz);
        emb      = emb(:, 1:2);
        method   = 'PCA';

    otherwise
        error('[cycif_plot_embedding] Unknown method ''%s''. Use ''umap'', ''tsne'', or ''pca''.', method_req);
end
fprintf('[cycif_plot_embedding] %s embedding complete.\n', method);

% --- Group / facet levels ------------------------------------------------
grp_levels = unique(T.(group_by), 'stable');
nGrp       = numel(grp_levels);
cmap       = lines(nGrp);

use_facet = ~isempty(facet_by);
if use_facet
    facet_levels = unique(T.(facet_by), 'stable');
else
    facet_levels = {''};
end
nFacet = numel(facet_levels);

% --- Figure --------------------------------------------------------------
if isempty(fig_pos)
    fig_pos = [80 80 min(520*nFacet, 2400) 480];
end

fig = figure('Color','w', 'Position', fig_pos, ...
             'Name', sprintf('%s | group=%s facet=%s', method, group_by, facet_by));

for fi = 1:nFacet
    if use_facet
        sel         = T.(facet_by) == facet_levels(fi);
        panel_title = string(facet_levels(fi));
    else
        sel         = true(height(T), 1);
        panel_title = method;
    end

    ax = subplot(1, nFacet, fi);
    hold(ax, 'on');

    emb_sel = emb(sel, :);
    grp_sel = T.(group_by)(sel);

    % Plot in reverse order so first group sits on top
    for gi = nGrp:-1:1
        s = grp_sel == grp_levels(gi);
        scatter(ax, emb_sel(s,1), emb_sel(s,2), mksz, cmap(gi,:), 'filled', ...
                'MarkerFaceAlpha', alpha, 'DisplayName', char(grp_levels(gi)));
    end

    title(ax, panel_title, 'Interpreter','none', 'FontSize', 10);
    xlabel(ax, [method ' 1']);
    if fi == 1, ylabel(ax, [method ' 2']); end
    axis(ax, 'equal');
    grid(ax, 'on'); box(ax, 'off');
    ax.FontSize = 9;

    % Legend on last panel only
    if fi == nFacet
        legend(ax, 'Location','eastoutside', 'Interpreter','none', 'FontSize', 8);
    end

    hold(ax, 'off');
end

% Subtitle notes transform and spike status
xform_str = 'log_2'; if ~use_log2, xform_str = 'raw'; end
spike_str  = 'spikes excluded'; if ~exclude_spike, spike_str = 'spikes included'; end

sgtitle(fig, sprintf('%s of CyCIF markers  |  group: %s  |  facet: %s  |  %s  |  %s', ...
        method, group_by, facet_by, xform_str, spike_str), ...
        'FontWeight','bold', 'FontSize', 10);
end


% =========================================================================
%  PRIVATE HELPERS
% =========================================================================

function markers = i_auto_markers(T)
% Detect marker columns as any numeric column that is not xy or cellid.
% Using type detection rather than a hardcoded exclusion list so that
% dynamic tx* columns added by cycif_build_table are always excluded.
    cols    = T.Properties.VariableNames;
    is_num  = cellfun(@(c) isnumeric(T.(c)) || islogical(T.(c)), cols);
    is_id   = ismember(cols, {'xy','cellid','ptx_is_spike'});
    markers = cols(is_num & ~is_id);
end