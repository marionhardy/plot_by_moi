function [T, cluster_labels] = cycif_cluster(T, emb, varargin)
% cycif_cluster  Cluster CyCIF cells on a 2D embedding (tSNE/PCA/UMAP),
%                visualize results, and add cluster labels to T.
%
% INPUT
%   T   — table from cycif_build_table (full, pre-filter)
%   emb — [nEmb x 2] embedding coordinates (from cycif_plot_embedding)
%         rows must correspond to T rows after spike filtering (no subsampling).
%
% NAME-VALUE OPTIONS
%   'method'        string    'kmeans' | 'dbscan' | 'louvain' | 'leiden' | 'both'
%                             (default 'both'). louvain = self-contained graph
%                             clustering; leiden = faster Python-bridge version
%                             (needs leidenalg; pass py_handler to reuse an env).
%   'py_handler'    handle    optional PythonHandler/UmapHandler for leiden
%   'min_cluster_size' scalar clusters with fewer cells are relabeled "small"
%                             (shown as one bracketed annotation row, not merged
%                             into neighbors). Cleans singleton sprays from an
%                             over-split resolution without coarsening good
%                             clusters. Default 1 (off).
%   'features'      matrix    [nEmb x nFeat] marker/PCA features to CLUSTER on
%                             (e.g. X from cycif_prep_features). [] = cluster on
%                             emb (legacy). The scatter always uses emb. Prefer
%                             passing features: the 2-D embedding has meaningless
%                             inter-cluster distances, so clustering on it is
%                             unreliable (that is what this whole refactor fixes).
%   'graph_k'       scalar    louvain kNN-graph neighbors    (default 15)
%   'gamma'         scalar    louvain resolution (>1 finer)  (default 1.0)
%   'seed'          scalar    RNG for kmeans/louvain         (default 0)
%   'k'             scalar    number of clusters for k-means  (default 8)
%   'epsilon'       scalar    DBSCAN neighborhood radius      (default auto)
%   'eps_scale'     scalar    multiplier on median k-NN dist for auto epsilon
%                             (default 1.5). Raise (e.g. 2-3) to merge sparse,
%                             strongly-polarized cells into clusters instead of
%                             dropping them to noise. Ignored if 'epsilon' set.
%   'zscore_rows'   logical   z-score marker columns across cluster-means
%                             (default true). false = plot mean log2 intensity
%                             directly (full magnitude, sequential colormap).
%   'clim_pct'      [lo hi]   robust CLim from data percentiles, e.g. [2 98]
%                             ([] = full range). Tames outlier-driven flattening.
%   'min_pts'       scalar    DBSCAN min neighbors            (default 10)
%   'emb_cell_idx'  table     [nEmb x 2] with xy + cellid columns, the 3rd
%                             output of cycif_plot_embedding. When supplied,
%                             cluster labels are placed back into T by MATCHING
%                             (xy, cellid) -- robust to any row filtering,
%                             subsampling, or reordering the embedding applied.
%                             This is REQUIRED in practice: cycif_plot_embedding
%                             applies a per-cell completeness filter
%                             (min_marker_frac / keep_cells, drops cells measured
%                             in < ~50% of markers before imputation) that
%                             cycif_cluster has no way to reconstruct -- so the
%                             row counts will not match without this. STRONGLY
%                             PREFERRED for every call.
%                             (default [] = fall back to spike-filter matching,
%                             which only works if NO completeness filter ran)
%   'exclude_spike' logical   drop ptx_is_spike wells         (default true)
%                             ONLY used for the fallback row-matching when
%                             emb_cell_idx is not supplied; must then match the
%                             setting used in cycif_plot_embedding. Ignored when
%                             emb_cell_idx is given.
%   'group_by'      string    T column for condition panel     (default 'tx1_label')
%   'colormap'      string    cluster colormap                 (default 'bright')
%   'center_zero'   logical   pin heatmap white to 0 via symmetric CLim
%                             (default true). false = imagesc autoscale.
%   'alpha'         scalar    scatter alpha                    (default 0.4)
%   'mksz'          scalar    scatter marker size              (default 8)
%   'fig_pos'       [x y w h] figure position                 (default auto)
%   'save_path'     string    base path for SVG export         (default '' = no save)
%
% OUTPUT
%   T              — input table with added 'cluster' column (categorical)
%                    cells not in emb receive cluster = 'unassigned'
%   cluster_labels — [nEmb x 1] string array, one label per emb row
%
% USAGE
%   [fig, emb, idx] = cycif_plot_embedding(T, 'method','tsne', 'seed',0);
%   [T, labels] = cycif_cluster(T, emb, 'emb_cell_idx', idx);   % PREFERRED
%   [T, labels] = cycif_cluster(T, emb, 'method','kmeans', 'k',9, 'emb_cell_idx', idx);
%   [T, labels] = cycif_cluster(T, emb, 'method','dbscan', 'epsilon',3);   % fallback

p = inputParser();
addParameter(p, 'exclude_spike', true);
addParameter(p, 'emb_cell_idx',  []);
addParameter(p, 'features',   []);   % [nEmb x nFeat] marker/PCA features to CLUSTER on
                                     % (e.g. cycif_prep_features X); [] = cluster on emb.
                                     % emb is always used for the scatter regardless.
addParameter(p, 'method',    'both');
addParameter(p, 'graph_k',   15);    % louvain/leiden: neighbors in kNN graph
addParameter(p, 'gamma',     1.0);   % louvain/leiden: resolution (>1 finer)
addParameter(p, 'py_handler', []);   % optional PythonHandler/UmapHandler for leiden
addParameter(p, 'min_cluster_size', 1);   % clusters below this -> "small" annotation row
addParameter(p, 'seed',      0);     % louvain/kmeans RNG
addParameter(p, 'k',         7);
addParameter(p, 'epsilon',   []);
addParameter(p, 'eps_scale', 1.5);   % multiplier on median k-NN dist for auto epsilon;
                                     % raise to merge sparse polarized cells into clusters
addParameter(p, 'min_pts',   10);
addParameter(p, 'group_by',  'tx1_label');
addParameter(p, 'colormap',  'bright');
addParameter(p, 'center_zero', true);   % pin heatmap white to 0 (symmetric CLim)
addParameter(p, 'zscore_rows', true);   % z-score marker columns across cluster-means;
                                        % false = plot mean log2 intensity directly
addParameter(p, 'clim_pct',    []);     % [lo hi] percentiles for robust CLim, e.g. [2 98];
                                        % [] = full range (symmetric if center_zero)
addParameter(p, 'alpha',     0.4);
addParameter(p, 'mksz',      8);
addParameter(p, 'fig_pos',   []);
addParameter(p, 'save_path', '');
parse(p, varargin{:});

exclude_spike = p.Results.exclude_spike;
emb_cell_idx  = p.Results.emb_cell_idx;
features      = p.Results.features;
method    = lower(p.Results.method);
graph_k   = p.Results.graph_k;
gamma     = p.Results.gamma;
py_handler = p.Results.py_handler;
min_cluster_size = p.Results.min_cluster_size;
seed      = p.Results.seed;
k         = p.Results.k;
epsilon   = p.Results.epsilon;
eps_scale = p.Results.eps_scale;
min_pts   = p.Results.min_pts;
group_by  = p.Results.group_by;
palette   = p.Results.colormap;
center_zero = p.Results.center_zero;
zscore_rows = p.Results.zscore_rows;
clim_pct    = p.Results.clim_pct;
alpha     = p.Results.alpha;
mksz      = p.Results.mksz;
fig_pos   = p.Results.fig_pos;
save_path = p.Results.save_path;

nEmb = size(emb, 1);

% --- Clustering space ----------------------------------------------------
% Cluster on marker/PCA features when provided (correct: the embedding is a
% lossy 2-D projection with meaningless inter-cluster distances). Fall back to
% clustering on emb for back-compat when features are not given. The scatter
% always uses emb regardless.
if ~isempty(features)
    assert(size(features,1) == nEmb, 'cycif_cluster:featRows', ...
        'features has %d rows but emb has %d; they must be row-aligned.', ...
        size(features,1), nEmb);
    Xc = features;
    fprintf('[cycif_cluster] Clustering on features [%d x %d] (embedding used for viz only).\n', ...
            size(Xc,1), size(Xc,2));
else
    Xc = emb;
    fprintf('[cycif_cluster] Clustering on embedding [%d x %d] (no features passed).\n', ...
            size(Xc,1), size(Xc,2));
end

% --- Auto epsilon for DBSCAN ---------------------------------------------
if isempty(epsilon) && ismember(method, {'dbscan','both'})
    knn_d   = pdist2(Xc, Xc, 'euclidean', 'Smallest', min_pts+1);
    epsilon = median(knn_d(end,:)) * eps_scale;
    fprintf('[cycif_cluster] DBSCAN auto epsilon = %.3f (%.2gx median %d-NN dist)\n', ...
            epsilon, eps_scale, min_pts);
end

% --- Run clustering ------------------------------------------------------
lou_labels = [];
switch method
    case 'kmeans'
        km_labels = i_run_kmeans(Xc, k, seed);
        db_labels = [];
    case 'dbscan'
        db_labels = i_run_dbscan(Xc, epsilon, min_pts);
        km_labels = [];
    case 'louvain'
        lou_labels = i_run_louvain(Xc, graph_k, gamma, seed);
        km_labels = []; db_labels = [];
    case 'leiden'
        lou_labels = i_run_leiden(Xc, graph_k, gamma, seed, py_handler);
        km_labels = []; db_labels = [];
    case 'both'
        km_labels = i_run_kmeans(Xc, k, seed);
        db_labels = i_run_dbscan(Xc, epsilon, min_pts);
    otherwise
        error('[cycif_cluster] Unknown method ''%s'' (kmeans|dbscan|louvain|both).', method);
end

if ~isempty(lou_labels)
    cluster_labels = lou_labels;
    if strcmp(method,'leiden'), primary_method = 'Leiden'; else, primary_method = 'Louvain'; end
elseif ~isempty(km_labels)
    cluster_labels = km_labels;
    primary_method = 'k-means';
else
    cluster_labels = db_labels;
    primary_method = 'DBSCAN';
end

fprintf('[cycif_cluster] Primary labels from %s: %d unique clusters\n', ...
        primary_method, numel(unique(cluster_labels)));

% --- Minimum cluster size ------------------------------------------------
% Relabel clusters below min_cluster_size to "small" (a non-cluster label,
% handled like noise/unassigned: excluded from the cluster row block, shown as
% a bracketed annotation row). This cleans up singleton/fragment sprays from an
% over-split resolution WITHOUT globally coarsening the good clusters -- prefer
% this over lowering gamma when the large clusters are already well-resolved.
% NOTE: 'small' cells are pooled into ONE annotation row; they are not merged
% into neighboring clusters (no arbitrary reassignment).
if min_cluster_size > 1
    cl = string(cluster_labels);
    u  = unique(cl);
    cnt = arrayfun(@(x) sum(cl == x), u);
    tooSmall = u(cnt < min_cluster_size);
    if ~isempty(tooSmall)
        nCells = sum(ismember(cl, tooSmall));
        cl(ismember(cl, tooSmall)) = "small";
        cluster_labels = cl;
        fprintf('[cycif_cluster] min_cluster_size=%d: %d clusters (%d cells) relabeled "small".\n', ...
                min_cluster_size, numel(tooSmall), nCells);
        fprintf('[cycif_cluster] %d clusters remain above threshold.\n', ...
                numel(unique(cl(cl ~= "small"))));
    end
end

% --- Add cluster column to T ---------------------------------------------
% Two placement strategies:
%   (A) PREFERRED -- join by (xy, cellid) when emb_cell_idx is supplied. This
%       is robust to whatever row filtering / subsampling / reordering the
%       embedding applied, because it matches on cell identity rather than
%       reconstructing the filter. In particular cycif_plot_embedding applies a
%       per-cell COMPLETENESS filter (min_marker_frac, default ~0.5, or an
%       explicit keep_cells list) that drops low-followability cells before
%       imputation -- cluster cannot reconstruct that, so this join is the only
%       correct way to align the two.
%   (B) FALLBACK -- reconstruct the spike-only filter by row order. Only valid
%       when the embedding applied exactly the spike filter and nothing else
%       (no completeness filter, no subsample). Errors with a clear message
%       otherwise.
T.cluster = repmat("unassigned", height(T), 1);

if ~isempty(emb_cell_idx)
    % --- (A) identity join ---
    % Prefer 'uid' when BOTH emb_cell_idx and T carry it: for a pooled
    % multi-replicate table (from cycif_pool_replicates) (xy,cellid) repeats
    % across replicates, so joining on (xy,cellid) alone would map every
    % replicate's cell to the first replicate's row -- silently corrupting
    % labels across reps. uid (rep_xy_cellid) is globally unique. Fall back to
    % (xy,cellid) for single-experiment tables that have no uid.
    idx_vars = emb_cell_idx.Properties.VariableNames;
    if ismember('uid', idx_vars) && ismember('uid', T.Properties.VariableNames)
        join_keys = {'uid'};
    else
        join_keys = {'xy','cellid'};
    end
    assert(istable(emb_cell_idx) && all(ismember(join_keys, idx_vars)), ...
        '[cycif_cluster] emb_cell_idx must be a table with %s columns (3rd output of cycif_plot_embedding).', ...
        strjoin(join_keys, '+'));
    assert(height(emb_cell_idx) == nEmb, ...
        ['[cycif_cluster] emb_cell_idx has %d rows but emb has %d. They must be the ' ...
         'matched pair returned together by cycif_plot_embedding.'], ...
        height(emb_cell_idx), nEmb);
    assert(all(ismember(join_keys, T.Properties.VariableNames)), ...
        '[cycif_cluster] T must contain %s to join embedding rows back.', strjoin(join_keys, '+'));

    % Map each emb row to its row in T via the join key(s).
    [tf, locT] = ismember(emb_cell_idx(:, join_keys), T(:, join_keys));
    if ~all(tf)
        warning(['[cycif_cluster] %d of %d embedding cells did not match a (%s) ' ...
                 'row in T -- those labels are dropped. Is this the same T passed to ' ...
                 'cycif_plot_embedding?'], sum(~tf), nEmb, strjoin(join_keys, ','));
    end
    T.cluster(locT(tf)) = cluster_labels(tf);
    n_assigned = sum(tf);

else
    % --- (B) fallback: spike-filter row-order reconstruction ---
    if exclude_spike && ismember('ptx_is_spike', T.Properties.VariableNames)
        emb_rows = find(~T.ptx_is_spike);
        row_src  = 'after spike filter';
    else
        emb_rows = (1:height(T))';
        row_src  = 'unfiltered (exclude_spike off, or no ptx_is_spike column)';
    end

    assert(numel(emb_rows) == nEmb, ...
        ['[cycif_cluster] Row mismatch: T has %d rows (%s) but emb has %d rows.\n' ...
         'The %d-row gap is rows the EMBEDDING dropped that this fallback cannot\n' ...
         'reconstruct -- most likely cycif_plot_embedding''s per-cell COMPLETENESS\n' ...
         'filter (min_marker_frac / keep_cells drops low-followability cells), or a\n' ...
         'max_pts subsample. PASS THE 3RD OUTPUT of cycif_plot_embedding:\n' ...
         '  [fig, emb, idx] = cycif_plot_embedding(T, ...);\n' ...
         '  [T, labels] = cycif_cluster(T, emb, ''emb_cell_idx'', idx);'], ...
        numel(emb_rows), row_src, nEmb, abs(numel(emb_rows) - nEmb));

    T.cluster(emb_rows) = cluster_labels;
    n_assigned = numel(emb_rows);
end

T.cluster = categorical(T.cluster);
fprintf('[cycif_cluster] Cluster column added to T (%d/%d cells assigned).\n', ...
        n_assigned, height(T));

% --- Markers for heatmap -------------------------------------------------
markers = i_auto_markers(T);

% --- Visualize -----------------------------------------------------------
if strcmp(method, 'both')
    i_plot_comparison(emb, km_labels, db_labels, T, group_by, ...
                      markers, palette, alpha, mksz, fig_pos, save_path, center_zero, zscore_rows, clim_pct);
else
    i_plot_single(emb, cluster_labels, primary_method, T, group_by, ...
                  markers, palette, alpha, mksz, fig_pos, save_path, center_zero, zscore_rows, clim_pct);
end
end


% =========================================================================
%  CLUSTERING HELPERS
% =========================================================================

function labels = i_run_kmeans(X, k, seed)
    if nargin >= 3 && ~isempty(seed), rng(seed); end
    fprintf('[cycif_cluster] Running k-means (k=%d)...\n', k);
    idx    = kmeans(X, k, 'Replicates', 5, 'MaxIter', 500);
    labels = string(idx);
    fprintf('[cycif_cluster] k-means done. Cluster sizes: ');
    for c = 1:k, fprintf('%d:%d  ', c, sum(idx==c)); end
    fprintf('\n');
end

function labels = i_run_louvain(X, graph_k, gamma, seed)
    % Graph-based clustering on feature space via self-contained Louvain.
    % No noise label -- Louvain assigns every cell to a community.
    fprintf('[cycif_cluster] Running Louvain (graph_k=%d, gamma=%.3g)...\n', graph_k, gamma);
    idx    = cycif_louvain(X, 'graph_k', graph_k, 'gamma', gamma, 'seed', seed, 'verbose', true);
    labels = string(idx);
end

function labels = i_run_leiden(X, graph_k, gamma, seed, py_handler)
    % Graph-based clustering via Leiden (Python bridge). Same graph as Louvain,
    % faster + better connectivity guarantees. No noise label.
    fprintf('[cycif_cluster] Running Leiden (graph_k=%d, gamma=%.3g)...\n', graph_k, gamma);
    args = {'graph_k', graph_k, 'gamma', gamma, 'seed', seed, 'verbose', true};
    if ~isempty(py_handler), args = [args, {'python_handler', py_handler}]; end
    idx    = cycif_leiden(X, args{:});
    labels = string(idx);
end

function labels = i_run_dbscan(X, epsilon, min_pts)
    fprintf('[cycif_cluster] Running DBSCAN (epsilon=%.3f, min_pts=%d)...\n', ...
            epsilon, min_pts);
    idx    = dbscan(X, epsilon, min_pts);
    labels = string(idx);
    labels(idx == -1) = "noise";
    nClust = numel(unique(idx(idx > 0)));
    nNoise = sum(idx == -1);
    fprintf('[cycif_cluster] DBSCAN done: %d clusters, %d noise points (%.1f%%)\n', ...
            nClust, nNoise, 100*nNoise/numel(idx));
end


% =========================================================================
%  PLOTTING HELPERS
% =========================================================================

function i_plot_comparison(emb, km_labels, db_labels, T, group_by, ...
                            markers, palette, alpha, mksz, fig_pos, save_path, center_zero, zscore_rows, clim_pct)
    [xlims, ylims] = i_emb_limits(emb);
    km_uniq  = i_sort_labels(unique(km_labels, 'stable'));
    db_uniq  = i_sort_labels(unique(db_labels, 'stable'));
    cmap_km  = i_colormap(numel(km_uniq), palette);
    cmap_db  = i_colormap(numel(db_uniq), palette);

    % Condition labels from T (rows matching emb)
    T_emb = T(T.cluster ~= 'unassigned', :);
    if height(T_emb) == size(emb,1)
        cond_labels = string(T_emb.(group_by));
    else
        warning('[cycif_cluster] T assigned rows != emb rows — condition panel skipped.');
        cond_labels = repmat("?", size(emb,1), 1);
    end
    grp_uniq = i_sort_labels(unique(cond_labels, 'stable'));
    cmap_grp = i_colormap(numel(grp_uniq), palette);

    if isempty(fig_pos), fig_pos = [40 40 1600 500]; end
    fig1 = figure('Color','w', 'Position', fig_pos, ...
                  'Name', 'CyCIF Clustering Comparison');

    % Explicitly parent all subplots to fig1
    ax1 = subplot(1,3,1, 'Parent', fig1);
    i_scatter_clusters(ax1, emb, km_labels, km_uniq, cmap_km, xlims, ylims, alpha, mksz, 'k-means');

    ax2 = subplot(1,3,2, 'Parent', fig1);
    i_scatter_clusters(ax2, emb, db_labels, db_uniq, cmap_db, xlims, ylims, alpha, mksz, 'DBSCAN');

    ax3 = subplot(1,3,3, 'Parent', fig1);
    i_scatter_clusters(ax3, emb, cond_labels, grp_uniq, cmap_grp, xlims, ylims, alpha, mksz, group_by);

    sgtitle(fig1, 'CyCIF Clustering: k-means vs DBSCAN vs condition', ...
            'FontWeight','bold', 'Interpreter','none');
    i_save_svg(fig1, save_path, 'clusters');

    i_plot_heatmap(T, km_labels, km_uniq, markers, 'k-means', fig_pos, save_path, center_zero, zscore_rows, clim_pct);
    i_plot_heatmap(T, db_labels, db_uniq, markers, 'DBSCAN',  fig_pos, save_path, center_zero, zscore_rows, clim_pct);
end

function i_plot_single(emb, labels, method_name, T, group_by, ...
                       markers, palette, alpha, mksz, fig_pos, save_path, center_zero, zscore_rows, clim_pct)
    [xlims, ylims] = i_emb_limits(emb);
    uniq     = i_sort_labels(unique(labels, 'stable'));
    cmap     = i_colormap(numel(uniq), palette);

    if isempty(fig_pos), fig_pos = [40 40 1100 500]; end
    fig1 = figure('Color','w', 'Position', fig_pos, ...
                  'Name', ['CyCIF Clustering — ' method_name]);

    % Explicitly parent all subplots to fig1
    ax1 = subplot(1,2,1, 'Parent', fig1);
    i_scatter_clusters(ax1, emb, labels, uniq, cmap, xlims, ylims, alpha, mksz, method_name);

    % Condition panel
    T_emb = T(T.cluster ~= 'unassigned', :);
    if height(T_emb) == size(emb,1)
        cond_lbl = string(T_emb.(group_by));
    else
        cond_lbl = repmat("?", size(emb,1), 1);
    end
    grp_uniq = i_sort_labels(unique(cond_lbl, 'stable'));
    cmap_grp = i_colormap(numel(grp_uniq), palette);

    ax2 = subplot(1,2,2, 'Parent', fig1);
    i_scatter_clusters(ax2, emb, cond_lbl, grp_uniq, cmap_grp, xlims, ylims, alpha, mksz, group_by);

    sgtitle(fig1, ['CyCIF Clustering — ' method_name], ...
            'FontWeight','bold', 'Interpreter','none');
    i_save_svg(fig1, save_path, 'clusters');

    i_plot_heatmap(T, labels, uniq, markers, method_name, fig_pos, save_path, center_zero, zscore_rows, clim_pct);
end

function i_scatter_clusters(ax, emb, labels, uniq, cmap, xlims, ylims, alpha, mksz, ttl)
    hold(ax, 'on');
    for ci = numel(uniq):-1:1
        s = labels == uniq(ci);
        scatter(ax, emb(s,1), emb(s,2), mksz, cmap(ci,:), 'filled', ...
                'MarkerFaceAlpha', alpha, 'DisplayName', char(uniq(ci)));
    end
    title(ax, ttl, 'Interpreter','none', 'FontSize', 10);
    xlim(ax, xlims); ylim(ax, ylims);
    ax.DataAspectRatio = [1 1 1];
    grid(ax,'off'); box(ax,'off');
    set(ax, 'XTick',[], 'YTick',[]);
    if numel(uniq) <= 20
        leg = legend(ax, 'Location','eastoutside', 'Interpreter','none', 'FontSize',7);
        leg.Box = 'off';
    else
        legend(ax, 'off');
    end
    hold(ax, 'off');
end

function i_plot_heatmap(T, labels, uniq, markers, method_name, fig_pos, save_path, center_zero, zscore_rows, clim_pct)
% Hierarchical clustering on both rows (clusters) and columns (markers).
% Dendrograms are drawn alongside the heatmap.
    nC  = numel(uniq);
    nMk = numel(markers);

    % --- Build mean log2 matrix [nC x nMk] ---
    hmat   = NaN(nC, nMk);
    ncells = zeros(nC, 1);
    for ci = 1:nC
        rows = string(T.cluster) == uniq(ci);
        ncells(ci) = sum(rows);
        for mi = 1:nMk
            v = T.(markers{mi})(rows);
            v = v(~isnan(v));
            if ~isempty(v), hmat(ci,mi) = mean(log2(max(v,1))); end
        end
    end

    % Z-score per marker column (across cluster-means), or plot raw mean log2.
    % NOTE: z-scoring is across the nC cluster rows, so it shows each cluster's
    % marker level RELATIVE to other clusters -- absolute intensity is lost, and
    % dynamic range is bounded by the (often skewed) spread of cluster-means.
    % zscore_rows=false plots mean log2 intensity directly (full magnitude).
    if nargin < 9 || isempty(zscore_rows), zscore_rows = true; end
    if zscore_rows
        mu_col  = mean(hmat, 'omitnan');
        std_col = std(hmat,  'omitnan'); std_col(std_col==0) = 1;
        hmat_z  = (hmat - mu_col) ./ std_col;
    else
        hmat_z  = hmat;   % raw mean log2 intensity
    end
    hmat_z(isnan(hmat_z)) = 0;   % replace NaN with 0 for linkage

    % --- Hierarchical clustering on rows and columns ---
    if nC > 1
        Z_row    = linkage(hmat_z,   'ward', 'euclidean');
        row_ord  = optimalleaforder(Z_row, pdist(hmat_z));
    else
        Z_row   = []; row_ord = 1;
    end
    if nMk > 1
        Z_col    = linkage(hmat_z',  'ward', 'euclidean');
        col_ord  = optimalleaforder(Z_col, pdist(hmat_z'));
    else
        Z_col   = []; col_ord = 1;
    end

    hmat_plot  = hmat_z(row_ord, col_ord);
    % Tag non-cluster rows (noise/unassigned/non-positive IDs) so they read as
    % annotations, not clusters. i_sort_labels already places them last, but the
    % row dendrogram may reorder; tag by content to stay correct either way.
    row_labels = arrayfun(@(i) i_row_label(uniq(row_ord(i)), ncells(row_ord(i))), ...
                          (1:nC)', 'UniformOutput', false);
    col_labels = markers(col_ord);

    % --- Figure layout: col dendrogram | heatmap | colorbar ---
    fig_w = max(700, nMk*45 + 180);
    fig_h = max(350, nC*35  + 150);
    if isempty(fig_pos), fig_pos = [40 40 900 400]; end
    hfig  = figure('Color','w', ...
                   'Position', [fig_pos(1)+50, fig_pos(2)+50, fig_w, fig_h], ...
                   'Name', ['Cluster Heatmap — ' method_name]);

    % Normalized positions: [left bottom width height]
    dend_col_pos = [0.10, 0.82, 0.72, 0.12];   % column dendrogram (top)
    dend_row_pos = [0.02, 0.10, 0.08, 0.70];   % row dendrogram (left)
    heat_pos     = [0.10, 0.10, 0.72, 0.70];   % heatmap

    % Column dendrogram
    if ~isempty(Z_col)
        ax_dcol = axes(hfig, 'Position', dend_col_pos);
        [~, ~, col_perm] = dendrogram(ax_dcol, Z_col, 0, ...
                                      'Reorder', col_ord, 'Orientation','top');
        ax_dcol.XTick = []; ax_dcol.YTick = [];
        ax_dcol.Box = 'off'; ax_dcol.XColor = 'none';
        ax_dcol.Color = 'none';
    end

    % Row dendrogram
    if ~isempty(Z_row)
        ax_drow = axes(hfig, 'Position', dend_row_pos);
        dendrogram(ax_drow, Z_row, 0, ...
                   'Reorder', row_ord, 'Orientation','right');
        ax_drow.XTick = []; ax_drow.YTick = [];
        ax_drow.Box = 'off'; ax_drow.YColor = 'none';
        ax_drow.Color = 'none';
        % Flip Y so row 1 is at top (matches imagesc)
        ax_drow.YDir = 'reverse';
    end

    % Heatmap
    ax = axes(hfig, 'Position', heat_pos);
    imagesc(ax, hmat_plot);

    % Colormap: diverging blue-white-red for z-scores (0 = white midpoint);
    % sequential for raw intensity, where a diverging map would imply a
    % meaningless zero-crossing.
    n2  = 128;
    if zscore_rows
        rr  = [linspace(0.17,1,n2), ones(1,n2)];
        gg  = [linspace(0.51,1,n2), linspace(1,0.30,n2)];
        bb  = [linspace(0.73,1,n2), linspace(1,0.17,n2)];
        colormap(ax, [rr', gg', bb']);
    else
        colormap(ax, parula(2*n2));
    end

    % --- Color limits ---
    % Priority: clim_pct (robust) > center_zero (symmetric about 0) > autoscale.
    % center_zero only makes sense for z-scores; in raw-intensity mode 0 is not a
    % meaningful midpoint, so it is ignored there.
    if nargin < 8 || isempty(center_zero), center_zero = true; end
    if nargin < 10, clim_pct = []; end
    vals = hmat_plot(:); vals = vals(isfinite(vals));

    if ~isempty(clim_pct)
        % Robust limits from data percentiles (guards outlier-driven range).
        lo = prctile(vals, clim_pct(1));
        hi = prctile(vals, clim_pct(2));
        if zscore_rows && center_zero
            cmax = max(abs([lo hi])); if cmax==0 || ~isfinite(cmax), cmax = 1; end
            clim(ax, [-cmax, cmax]);          % keep white on 0
        else
            if hi <= lo, hi = lo + 1; end
            clim(ax, [lo, hi]);
        end
    elseif zscore_rows && center_zero
        cmax = max(abs(vals)); if isempty(cmax) || cmax==0 || ~isfinite(cmax), cmax = 1; end
        clim(ax, [-cmax, cmax]);
    end   % else: leave imagesc autoscale (full [min max])

    cb = colorbar(ax);
    if zscore_rows
        cb.Label.String  = 'z-score (log_2 intensity)';
    else
        cb.Label.String  = 'mean log_2 intensity';
    end
    cb.Label.Interpreter = 'none';

    ax.XTick               = 1:nMk;
    ax.XTickLabel          = col_labels;
    ax.XTickLabelRotation  = 45;
    ax.YTick               = 1:nC;
    ax.YTickLabel          = row_labels;
    ax.FontSize            = 9;
    ax.TickLabelInterpreter = 'none';
    axis(ax, 'tight');
    title(ax, ['Marker means per cluster — ' method_name], ...
          'Interpreter','none', 'FontWeight','bold');

    i_save_svg(hfig, save_path, ['heatmap_' lower(strrep(method_name,' ','_'))]);
end


% =========================================================================
%  UTILITY HELPERS
% =========================================================================

function [xlims, ylims] = i_emb_limits(emb)
    xp = 0.05 * range(emb(:,1));
    yp = 0.05 * range(emb(:,2));
    xlims = [min(emb(:,1))-xp, max(emb(:,1))+xp];
    ylims = [min(emb(:,2))-yp, max(emb(:,2))+yp];
    ar    = max(diff(xlims), diff(ylims));
    xc    = mean(xlims); yc = mean(ylims);
    xlims = [xc-ar/2, xc+ar/2];
    ylims = [yc-ar/2, yc+ar/2];
end

function s = i_row_label(lbl, n)
% Row label with count; brackets non-cluster rows (noise/unassigned/<=0) so
% they are visually distinct from true cluster rows.
    lbl = string(lbl);
    num = str2double(lbl);
    isSpecial = lbl == "noise" || lbl == "unassigned" || lbl == "small" || (~isnan(num) && num <= 0);
    if isSpecial
        s = sprintf('[%s] (n=%d)', lbl, n);   % bracketed = not a cluster
    else
        s = sprintf('%s (n=%d)', lbl, n);
    end
end

function sorted = i_sort_labels(labels)
% Order true cluster rows first (ascending numeric), then non-cluster labels
% (noise / unassigned / non-positive IDs) as annotated bottom rows. This keeps
% DBSCAN noise and any stray '0'/negative label out of the CLUSTER block while
% still showing them for inspection.
    labels = string(labels);
    isSpecialName = labels == "noise" | labels == "unassigned" | labels == "small";
    num           = str2double(labels);
    isNonPosNum   = ~isnan(num) & num <= 0;      % 0 or negative IDs are not clusters
    isSpecial     = isSpecialName | isNonPosNum;

    clusters = labels(~isSpecial);
    special  = labels( isSpecial);
    [~, ix]  = sort(str2double(clusters));
    % de-duplicate special while preserving a stable order (noise, unassigned, 0, ...)
    special  = unique(special, 'stable');
    sorted   = [clusters(ix); special];
end

function i_save_svg(fig, save_path, suffix)
    if isempty(save_path), return; end
    fpath = [save_path '_' suffix '.svg'];
    fig.Units             = 'inches';
    fig_sz                = fig.Position(3:4);
    fig.PaperUnits        = 'inches';
    fig.PaperSize         = fig_sz;
    fig.PaperPosition     = [0, 0, fig_sz];
    fig.PaperPositionMode = 'manual';
    print(fig, fpath, '-dsvg', '-painters');
    fprintf('Saved SVG: %s\n', fpath);
end

function markers = i_auto_markers(T)
    cols   = T.Properties.VariableNames;
    is_num = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
    is_id  = ismember(cols, {'xy','cellid','ptx_is_spike','cluster'});
    markers = cols(is_num & ~is_id);
end

function cmap = i_colormap(n, palette)
if nargin < 2 || isempty(palette), palette = 'bright'; end
switch lower(palette)
    case 'bright'
        base = [0.900,0.100,0.100; 0.100,0.500,0.900; 0.100,0.850,0.300;
                0.000,0.620,0.200; 0.000,0.390,0.100; 0.900,0.100,0.550;
                1.000,0.700,0.000; 0.960,0.480,0.000; 0.980,0.860,0.100;
                0.400,0.400,0.400; 0.600,0.200,0.800; 0.200,0.800,0.800];
        if n <= size(base,1), cmap = base(1:n,:);
        else, cmap = lines(n); end
    case 'parula'
        raw = parula(max(n,2));
        cmap = raw(round(linspace(1,size(raw,1),n)),:);
    otherwise
        try
            raw  = feval(lower(palette), max(n,2));
            cmap = raw(round(linspace(1,size(raw,1),n)),:);
        catch
            cmap = lines(n);
        end
end
end