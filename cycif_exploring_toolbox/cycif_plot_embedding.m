function [fig, emb, emb_cell_idx] = cycif_plot_embedding(T, varargin)
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
%
% OUTPUT
%   fig          — figure handle
%   emb          — [nCells x 2] embedding coordinates
%                  rows correspond to T after spike filtering and subsampling
%   emb_cell_idx — table [nCells x 2] with xy and cellid columns
%                  use this to join emb rows back to dataloc or original T
%
% USAGE
%   cycif_plot_embedding(T);
%   [fig, emb_tsne, idx] = cycif_plot_embedding(T, 'method','tsne', 'seed',0);
%   cycif_plot_embedding(T, 'group_by','tx1_label', 'facet_by','ptx_pool_label');

p = inputParser();
addParameter(p, 'markers',       []);
addParameter(p, 'group_by',      'tx1_label');
addParameter(p, 'facet_by',      '');
addParameter(p, 'method',        'tsne');
addParameter(p, 'colormap',      'bright');
addParameter(p, 'use_log2',      true);
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
addParameter(p, 'perplexity',    30);
addParameter(p, 'seed',          0);
addParameter(p, 'fig_pos',       []);
addParameter(p, 'save_path',     '');
parse(p, varargin{:});

markers        = p.Results.markers;
group_by       = p.Results.group_by;
facet_by       = p.Results.facet_by;
method_req     = lower(p.Results.method);
palette        = p.Results.colormap;
use_log2       = p.Results.use_log2;
exclude_spike  = p.Results.exclude_spike;
drop_ctrl      = p.Results.drop_ctrl;
max_pts        = p.Results.max_pts;
n_neighbors    = p.Results.n_neighbors;
min_dist       = p.Results.min_dist;
alpha          = p.Results.alpha;
mksz           = p.Results.mksz;
max_nan_frac   = p.Results.max_nan_frac;
min_marker_frac = p.Results.min_marker_frac;
keep_cells     = p.Results.keep_cells;
impute_nan     = lower(p.Results.impute_nan);
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

% --- Build marker matrix -------------------------------------------------
Xmat = table2array(T(:, markers));
if use_log2
    nz = ~isnan(Xmat);
    Xmat(nz) = log2(max(Xmat(nz), 1));
end

% --- Cell completeness / followability filter (BEFORE imputation) --------
% Heavily-NaN cells (poor followability: measured/well-registered in few
% rounds) otherwise get imputed to identical per-marker constants, collapsing
% into one giant identical point. That breaks tSNE's perplexity binary search
% (zero-distance neighbors) and forms an artifactual central blob. Drop those
% cells here, before imputation, so only followable cells are embedded.
%
% Two mechanisms (keep_cells overrides min_marker_frac when provided):
%   keep_cells       explicit cell keep-list. Either logical [height(T)x1], or
%                    numeric row indices into T, or a table with xy+cellid
%                    (e.g. built from QC.survived_good >= threshold).
%   min_marker_frac  keep cells with at least this fraction of markers measured
%                    (non-NaN). Default 0.5.
frac_measured = mean(~isnan(Xmat), 2);
nCellsIn = numel(frac_measured);

% Completeness distribution (per-cell fraction of markers measured), so the
% threshold choice is evidence-based.
fprintf('[cycif_plot_embedding] --- Cell completeness (fraction of %d markers measured) ---\n', nMk);
bands = [0 0.25 0.50 0.75 1.0001];
blab  = {'  0-25%','  25-50%',' 50-75%','75-100%'};
for b = 1:numel(bands)-1
    inb = frac_measured >= bands(b) & frac_measured < bands(b+1);
    fprintf('    %s : %6d cells (%.1f%%)\n', blab{b}, sum(inb), 100*mean(inb));
end
fprintf('    fully measured (100%%): %6d cells (%.1f%%)\n', ...
        sum(frac_measured>=1), 100*mean(frac_measured>=1));
fprintf('    median completeness: %.0f%% | mean: %.0f%%\n', ...
        100*median(frac_measured), 100*mean(frac_measured));

if ~isempty(keep_cells)
    keep_row = local_resolve_keeplist(keep_cells, T);
    fprintf('[cycif_plot_embedding] External keep-list applied: %d of %d cells retained (%d dropped, %.1f%%).\n', ...
            sum(keep_row), nCellsIn, sum(~keep_row), 100*mean(~keep_row));
    fprintf('    completeness of KEPT cells:     median %.0f%% (min %.0f%%)\n', ...
            100*median(frac_measured(keep_row)), 100*min(frac_measured(keep_row)));
    if any(~keep_row)
        fprintf('    completeness of DROPPED cells:  median %.0f%% (max %.0f%%)\n', ...
                100*median(frac_measured(~keep_row)), 100*max(frac_measured(~keep_row)));
    end
else
    keep_row = frac_measured >= min_marker_frac;
    fprintf('[cycif_plot_embedding] Completeness filter (>=%.0f%% markers): %d of %d retained (%d dropped, %.1f%%).\n', ...
            100*min_marker_frac, sum(keep_row), nCellsIn, sum(~keep_row), 100*mean(~keep_row));
end

% Duplicate-collapse warning: report how many of the DROPPED cells would have
% imputed to (near-)identical rows -- the tSNE-breaking population.
if any(~keep_row)
    dmin = min(frac_measured(~keep_row)); dmed = median(frac_measured(~keep_row));
    fprintf('    (dropped cells are the low-completeness tail that imputes to identical rows;\n');
    fprintf('     dropped completeness range %.0f%%-%.0f%%, would otherwise collapse the embedding.)\n', ...
            100*dmin, 100*max(frac_measured(~keep_row)));
end

% Per-marker retained coverage: after the cell filter, how well is each marker
% actually measured? Flags markers that are still mostly imputed.
kept_cov = mean(~isnan(Xmat(keep_row,:)), 1);
[cov_sorted, ord] = sort(kept_cov);
nshow = min(5, nMk);
fprintf('    lowest-coverage markers among kept cells:\n');
for i = 1:nshow
    fprintf('       %-16s %.0f%% measured\n', markers{ord(i)}, 100*cov_sorted(i));
end

assert(sum(keep_row) >= 2, ...
    ['[cycif_plot_embedding] Only %d cells pass the completeness/keep filter; ' ...
     'lower ''min_marker_frac'' or widen ''keep_cells''.'], sum(keep_row));
T    = T(keep_row, :);
Xmat = Xmat(keep_row, :);
fprintf('[cycif_plot_embedding] Embedding %d cells x %d markers.\n', height(T), nMk);

% --- NaN imputation for embedding ---------------------------------------
% Unstained wells have NaN — they must be filled before embedding.
% Strategy is controlled by 'impute_nan' parameter:
%   'lowerquantile' — global 25th percentile per marker (default).
%                     Places unstained cells at the low end of the
%                     distribution without assuming zero expression.
%   'median'        — per-condition median. Circular but neutral.
%   'none'          — drop all cells with any NaN (may lose many cells).
n_imputed = 0;
switch impute_nan
    case 'lowerquantile'
        for m = 1:nMk
            col = Xmat(:, m);
            if ~any(isnan(col)), continue; end
            fill_val  = quantile(col(~isnan(col)), 0.25);
            if isnan(fill_val), fill_val = 0; end
            fill_idx  = isnan(col);
            col(fill_idx) = fill_val;
            n_imputed = n_imputed + sum(fill_idx);
            Xmat(:, m) = col;
        end
        fprintf('[cycif_plot_embedding] ASSUMPTION — lower-quartile imputation:\n');
        fprintf('  %d NaN values filled with per-marker global 25th percentile.\n', n_imputed);
        fprintf('  Unstained cells placed at low end of each marker distribution.\n');
        fprintf('  Embedding is exploratory only — do not use for classification.\n');

    case 'median'
        cond_labels = T.tx1_label;
        grp_imp     = unique(cond_labels, 'stable');
        for m = 1:nMk
            col = Xmat(:, m);
            if ~any(isnan(col)), continue; end
            for g = 1:numel(grp_imp)
                in_grp   = cond_labels == grp_imp(g);
                is_nan   = isnan(col);
                fill_idx = in_grp & is_nan;
                if ~any(fill_idx), continue; end
                grp_med = median(col(in_grp & ~is_nan), 'omitnan');
                if isnan(grp_med), grp_med = median(col(~is_nan), 'omitnan'); end
                if isnan(grp_med), grp_med = 0; end
                col(fill_idx) = grp_med;
                n_imputed = n_imputed + sum(fill_idx);
            end
            Xmat(:, m) = col;
        end
        fprintf('[cycif_plot_embedding] ASSUMPTION — per-condition median imputation:\n');
        fprintf('  %d NaN values filled. Condition label used — circular assumption.\n', n_imputed);
        fprintf('  Embedding is exploratory only — do not use for classification.\n');

    case 'none'
        valid = all(~isnan(Xmat), 2);
        n_dropped = sum(~valid);
        if n_dropped > 0
            fprintf('[cycif_plot_embedding] Dropping %d cells with any NaN (impute_nan=none).\n', n_dropped);
            T    = T(valid, :);
            Xmat = Xmat(valid, :);
        end

    otherwise
        error('[cycif_plot_embedding] Unknown impute_nan ''%s''. Use ''lowerquantile'', ''median'', or ''none''.', impute_nan);
end

% --- Subsample -----------------------------------------------------------
if height(T) > max_pts
    idx  = randperm(height(T), max_pts);
    T    = T(idx, :);
    Xmat = Xmat(idx, :);
    fprintf('[cycif_plot_embedding] Subsampled to %d cells.\n', max_pts);
end

% --- Z-score per marker --------------------------------------------------
mu  = mean(Xmat);
sig = std(Xmat);
sig(sig == 0) = 1;
Xz  = (Xmat - mu) ./ sig;

% --- Set RNG seed --------------------------------------------------------
rng(seed);

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

% --- Cell index for joining emb back to dataloc --------------------------
emb_cell_idx = T(:, {'xy','cellid'});

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
xform_str = 'log_2'; if ~use_log2, xform_str = 'raw'; end
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