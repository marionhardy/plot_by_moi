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
%   'method'        string    'kmeans' | 'dbscan' | 'both'   (default 'both')
%   'k'             scalar    number of clusters for k-means  (default 8)
%   'epsilon'       scalar    DBSCAN neighborhood radius      (default auto)
%   'min_pts'       scalar    DBSCAN min neighbors            (default 10)
%   'exclude_spike' logical   drop ptx_is_spike wells         (default true)
%                             must match setting in cycif_plot_embedding
%   'group_by'      string    T column for condition panel     (default 'tx1_label')
%   'colormap'      string    cluster colormap                 (default 'bright')
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
%   [fig, emb] = cycif_plot_embedding(T, 'method','tsne', 'seed',0);
%   [T, labels] = cycif_cluster(T, emb);
%   [T, labels] = cycif_cluster(T, emb, 'method','kmeans', 'k',9);
%   [T, labels] = cycif_cluster(T, emb, 'method','dbscan', 'epsilon',3);

p = inputParser();
addParameter(p, 'exclude_spike', true);
addParameter(p, 'method',    'both');
addParameter(p, 'k',         8);
addParameter(p, 'epsilon',   []);
addParameter(p, 'min_pts',   10);
addParameter(p, 'group_by',  'tx1_label');
addParameter(p, 'colormap',  'bright');
addParameter(p, 'alpha',     0.4);
addParameter(p, 'mksz',      8);
addParameter(p, 'fig_pos',   []);
addParameter(p, 'save_path', '');
parse(p, varargin{:});

exclude_spike = p.Results.exclude_spike;
method    = lower(p.Results.method);
k         = p.Results.k;
epsilon   = p.Results.epsilon;
min_pts   = p.Results.min_pts;
group_by  = p.Results.group_by;
palette   = p.Results.colormap;
alpha     = p.Results.alpha;
mksz      = p.Results.mksz;
fig_pos   = p.Results.fig_pos;
save_path = p.Results.save_path;

nEmb = size(emb, 1);

% --- Auto epsilon for DBSCAN ---------------------------------------------
if isempty(epsilon) && ismember(method, {'dbscan','both'})
    knn_d   = pdist2(emb, emb, 'euclidean', 'Smallest', min_pts+1);
    epsilon = median(knn_d(end,:)) * 1.5;
    fprintf('[cycif_cluster] DBSCAN auto epsilon = %.3f (1.5x median %d-NN dist)\n', ...
            epsilon, min_pts);
end

% --- Run clustering ------------------------------------------------------
switch method
    case 'kmeans'
        km_labels = i_run_kmeans(emb, k);
        db_labels = [];
    case 'dbscan'
        db_labels = i_run_dbscan(emb, epsilon, min_pts);
        km_labels = [];
    case 'both'
        km_labels = i_run_kmeans(emb, k);
        db_labels = i_run_dbscan(emb, epsilon, min_pts);
    otherwise
        error('[cycif_cluster] Unknown method ''%s''.', method);
end

if ~isempty(km_labels)
    cluster_labels = km_labels;
    primary_method = 'k-means';
else
    cluster_labels = db_labels;
    primary_method = 'DBSCAN';
end

fprintf('[cycif_cluster] Primary labels from %s: %d unique clusters\n', ...
        primary_method, numel(unique(cluster_labels)));

% --- Add cluster column to T by row-order matching -----------------------
T.cluster = repmat("unassigned", height(T), 1);

if exclude_spike && ismember('ptx_is_spike', T.Properties.VariableNames)
    emb_rows = find(~T.ptx_is_spike);
else
    emb_rows = (1:height(T))';
end

assert(numel(emb_rows) == nEmb, ...
    ['[cycif_cluster] Row mismatch: T has %d rows after spike filter but emb has %d rows.\n' ...
     'Ensure exclude_spike matches the setting in cycif_plot_embedding,\n' ...
     'and that no subsampling was applied (max_pts=Inf).'], ...
    numel(emb_rows), nEmb);

T.cluster(emb_rows) = cluster_labels;
T.cluster = categorical(T.cluster);
fprintf('[cycif_cluster] Cluster column added to T (%d/%d cells assigned).\n', ...
        sum(T.cluster ~= 'unassigned'), height(T));

% --- Markers for heatmap -------------------------------------------------
markers = i_auto_markers(T);

% --- Visualize -----------------------------------------------------------
if strcmp(method, 'both')
    i_plot_comparison(emb, km_labels, db_labels, T, group_by, ...
                      markers, palette, alpha, mksz, fig_pos, save_path);
else
    i_plot_single(emb, cluster_labels, primary_method, T, group_by, ...
                  markers, palette, alpha, mksz, fig_pos, save_path);
end
end


% =========================================================================
%  CLUSTERING HELPERS
% =========================================================================

function labels = i_run_kmeans(emb, k)
    fprintf('[cycif_cluster] Running k-means (k=%d)...\n', k);
    idx    = kmeans(emb, k, 'Replicates', 5, 'MaxIter', 500);
    labels = string(idx);
    fprintf('[cycif_cluster] k-means done. Cluster sizes: ');
    for c = 1:k, fprintf('%d:%d  ', c, sum(idx==c)); end
    fprintf('\n');
end

function labels = i_run_dbscan(emb, epsilon, min_pts)
    fprintf('[cycif_cluster] Running DBSCAN (epsilon=%.3f, min_pts=%d)...\n', ...
            epsilon, min_pts);
    idx    = dbscan(emb, epsilon, min_pts);
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
                            markers, palette, alpha, mksz, fig_pos, save_path)
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

    i_plot_heatmap(T, km_labels, km_uniq, markers, 'k-means', fig_pos, save_path);
    i_plot_heatmap(T, db_labels, db_uniq, markers, 'DBSCAN',  fig_pos, save_path);
end

function i_plot_single(emb, labels, method_name, T, group_by, ...
                       markers, palette, alpha, mksz, fig_pos, save_path)
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

    i_plot_heatmap(T, labels, uniq, markers, method_name, fig_pos, save_path);
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

function i_plot_heatmap(T, labels, uniq, markers, method_name, fig_pos, save_path)
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

    % Z-score per marker column
    mu_col  = mean(hmat, 'omitnan');
    std_col = std(hmat,  'omitnan'); std_col(std_col==0) = 1;
    hmat_z  = (hmat - mu_col) ./ std_col;
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
    row_labels = arrayfun(@(i) sprintf('%s (n=%d)', uniq(row_ord(i)), ncells(row_ord(i))), ...
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

    % Diverging blue-white-red colormap
    n2  = 128;
    rr  = [linspace(0.17,1,n2), ones(1,n2)];
    gg  = [linspace(0.51,1,n2), linspace(1,0.30,n2)];
    bb  = [linspace(0.73,1,n2), linspace(1,0.17,n2)];
    colormap(ax, [rr', gg', bb']);

    cb = colorbar(ax);
    cb.Label.String      = 'z-score (log_2 intensity)';
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

function sorted = i_sort_labels(labels)
    labels  = string(labels);
    numeric = labels(labels ~= "noise" & labels ~= "unassigned");
    special = labels(labels == "noise" | labels == "unassigned");
    [~, ix] = sort(str2double(numeric));
    sorted  = [numeric(ix); special];
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