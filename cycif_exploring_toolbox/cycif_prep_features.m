function [X, cell_idx, info] = cycif_prep_features(T, varargin)
% cycif_prep_features  Shared feature preparation for CyCIF embedding & clustering.
%
%   Single source of truth for turning a CyCIF cell table into a numeric feature
%   matrix for distance-based methods (tSNE / UMAP / k-means / Louvain). Both the
%   embedding and the clustering call this so they operate on IDENTICAL features.
%
%   Pipeline:  select markers -> transform -> completeness filter -> impute NaN
%              -> z-score -> (optional) PCA denoise.
%
%   IMPORTANT: every transform here is confined to the DISTANCE METRIC used for
%   embedding/clustering. It never alters reported intensities -- downstream
%   readouts (heatmaps, fold-changes, dot-plots) stay on LINEAR intensity, read
%   straight from T. This mirrors the CyTOF/CyCIF standard: arcsinh for the
%   neighbour graph, raw for quantitation.
%
% Usage examples:
%   [X, idx]        = cycif_prep_features(T);                        % arcsinh + 15 PCs
%   [X, idx, info]  = cycif_prep_features(T, 'transform','log2', 'n_pcs',0);
%   [X, idx, info]  = cycif_prep_features(T, 'n_pcs',15, 'plot_scree',true);  % tune n_pcs
%   [X, idx, info]  = cycif_prep_features(T, 'n_pcs',7, 'drop_pcs',[1]);      % remove amplitude PC
%   [X, idx, info]  = cycif_prep_features(T, 'cofactor',150, 'zscore',true);
%   % feed the SAME X to both:
%   emb    = tsne(X);                       % visualization
%   labels = kmeans(X, 8);                  % clustering
%
% INPUT
%   T  -- table from cycif_build_table (cells x variables)
%
% NAME-VALUE OPTIONS
%   'markers'      {1xM} cell   marker column names; [] = auto-detect  (default [])
%   'transform'    string       'arcsinh'|'log2'|'linear'             (default 'arcsinh')
%                               arcsinh: asinh(x/cofactor) -- linear near 0
%                                        (background), log-like for large x;
%                                        the CyCIF/CyTOF default.
%                               log2   : log2(max(x,1)) -- matches legacy embedding.
%                               linear : raw intensity (NOT recommended for a
%                                        Euclidean metric; loudest markers dominate).
%   'cofactor'     scalar       arcsinh cofactor (knee location)      (default 0.4)
%                               Set empirically from the background floor via
%                               cycif_feature_hist. 0.4 matches this dataset's
%                               per-marker background median (~0.43) AND the
%                               visual knee -- data is fold-over-background
%                               (background ~0.2-0.9), NOT raw fluorescence, so
%                               the classic 150 cofactor is wrong here. Re-check
%                               cycif_feature_hist for a different normalization.
%   'zscore'       logical      per-marker z-score after transform    (default true)
%   'normalize_cells' string    per-cell amplitude normalization BEFORE transform:
%                               'none'   -- off (default)
%                               'mean'   -- divide each cell by its mean marker level
%                               'median' -- divide each cell by its median marker level
%                               ONLY VALID for strictly-POSITIVE data (ratio /
%                               fold-over-background / raw intensity). ERRORS on
%                               background-SUBTRACTED data (values straddling 0) --
%                               use drop_pcs=[1] there instead. Induces
%                               compositional closure; interpret PCs as relative
%                               enrichment. Readouts from T stay linear.
%   'drop_pcs'     vector       PC indices to EXCLUDE from X after PCA (default []).
%                               The correct amplitude-axis control for signed
%                               (background-subtracted) data: if PC1 is a uniform-
%                               sign amplitude/brightness axis, drop_pcs=[1] removes
%                               it without any positivity assumption. info.explained
%                               and .coeff reflect the RETAINED PCs only.
%   'n_pcs'        scalar       PCA components to retain; 0 = no PCA   (default 15)
%                               PCA denoises + decorrelates before the neighbour
%                               graph. Capped at min(n_pcs, nMarkers, nCells-1).
%   'impute_nan'   string       'lowerquantile'|'median'|'none'       (default 'lowerquantile')
%                               NaN fill BEFORE z-score/PCA. 'none' drops any row
%                               with a remaining NaN.
%   'min_marker_frac' scalar    per-cell completeness floor (frac measured) (default 0.9)
%   'keep_cells'   varies       explicit keep-list (logical | idx | xy+cellid
%                               table); overrides min_marker_frac         (default [])
%   'plot_scree'   logical      show PCA scree (variance per PC + cumulative,
%                               n_pcs cut marked) to guide n_pcs choice  (default false)
%   'plot_corr'    logical      show marker correlation heatmap of the PCA-input
%                               features (diagnoses amplitude vs block structure)
%                               (default false)
%   'verbose'      logical      print completeness / PCA diagnostics   (default true)
%
% OUTPUT
%   X         -- [nKept x D]  feature matrix. D = n_pcs (PCA on) or nMarkers (off).
%   cell_idx  -- table [nKept x K] identity keys for joining X rows back to T.
%                Always includes xy, cellid; includes uid when T has it. This is
%                the SAME contract cycif_cluster's emb_cell_idx join expects, so
%                the third output of an embedding built on X aligns 1:1.
%   info      -- struct:
%                  .markers       markers used
%                  .transform     transform name
%                  .cofactor      arcsinh cofactor (NaN if not arcsinh)
%                  .zscore        logical
%                  .n_pcs         PCs actually retained (post-cap; 0 if off)
%                  .explained     [n_pcs x 1] % variance per retained PC ([] if off)
%                  .cum_explained cumulative % variance at n_pcs ([] if off)
%                  .keep_row      logical mask into original T rows
%                  .mu,.sig       z-score center/scale (per marker; [] if zscore off)
%                  .coeff         PCA loadings [nMarkers x n_pcs] ([] if off)
%
% See also: cycif_feature_hist (cofactor/transform diagnostic), cycif_plot_embedding.

% ---- parse ----
p = inputParser();
addParameter(p, 'markers',          []);
addParameter(p, 'transform',        'arcsinh');
addParameter(p, 'cofactor',         0.4);
addParameter(p, 'zscore',           true);
addParameter(p, 'normalize_cells',  'none');   % 'none'|'mean'|'median': per-cell
                                               % amplitude normalization BEFORE
                                               % transform, to remove the technical
                                               % brightness/size axis (see header)
addParameter(p, 'n_pcs',            15);
addParameter(p, 'drop_pcs',         []);   % PC indices to EXCLUDE from X after PCA,
                                           % e.g. [1] to remove an amplitude axis
                                           % (correct amplitude control for signed data)
addParameter(p, 'impute_nan',       'lowerquantile');
addParameter(p, 'min_marker_frac',  0.9);
addParameter(p, 'keep_cells',       []);
addParameter(p, 'plot_scree',       false);   % opt-in PCA scree plot for n_pcs selection
addParameter(p, 'plot_corr',        false);   % opt-in marker correlation heatmap
addParameter(p, 'verbose',          true);
parse(p, varargin{:});

markers    = p.Results.markers;
transform  = lower(string(p.Results.transform));
cofactor   = p.Results.cofactor;
do_zscore  = logical(p.Results.zscore);
norm_cells = lower(string(p.Results.normalize_cells));
n_pcs      = p.Results.n_pcs;
drop_pcs   = p.Results.drop_pcs(:)';
impute_nan = lower(string(p.Results.impute_nan));
min_frac   = p.Results.min_marker_frac;
keep_cells = p.Results.keep_cells;
plot_scree = logical(p.Results.plot_scree);
plot_corr  = logical(p.Results.plot_corr);
verbose    = logical(p.Results.verbose);

assert(ismember(transform, ["arcsinh","log2","linear"]), ...
    'cycif_prep_features:transform', ...
    '''transform'' must be arcsinh | log2 | linear (got ''%s'').', transform);
assert(ismember(impute_nan, ["lowerquantile","median","none"]), ...
    'cycif_prep_features:impute', ...
    '''impute_nan'' must be lowerquantile | median | none (got ''%s'').', impute_nan);

% ---- markers ----
if isempty(markers)
    markers = i_auto_markers(T);
end
markers = cellstr(markers);
nMk = numel(markers);
assert(nMk >= 2, 'cycif_prep_features:markers', 'Need >=2 marker columns (found %d).', nMk);

Xraw = table2array(T(:, markers));   % LINEAR intensities, cells x markers

% ---- per-cell amplitude normalization (BEFORE transform, on linear) ----
% Divides each cell by its own central marker level, removing the technical
% brightness/size axis (the uniform-positive PC1 = 55% variance seen in this
% data). Operates on linear intensities so the ratio is meaningful.
%
% CAVEAT (compositional closure): after per-cell normalization the markers are
% no longer independent -- scaling each cell to a common level forces a mild
% negative dependence among markers (if one rises as a fraction, others fall).
% This can induce SPURIOUS anti-correlation. It is the accepted trade for
% removing a dominant amplitude confound; interpret resulting PC loadings as
% RELATIVE enrichment (compositional), not absolute intensity. Readouts from T
% remain linear and are unaffected.
switch norm_cells
    case "none"
        % leave Xraw as-is
    case {"mean","median"}
        % VALIDITY GUARD: per-cell division is only meaningful for strictly
        % POSITIVE data (ratio / fold-over-background / raw intensity). On
        % background-SUBTRACTED data (values straddling 0), the per-cell mean
        % sits near zero and division is nonsensical (sign flips, blow-ups).
        % Detect and refuse rather than silently produce garbage.
        neg_frac = mean(Xraw(:) < 0, 'omitnan');
        if neg_frac > 0.01
            error('cycif_prep_features:signedNorm', ...
                ['normalize_cells=''%s'' is invalid: %.1f%% of marker values are ' ...
                 'negative, so this data is background-SUBTRACTED, not positive ' ...
                 'ratio data. Per-cell division requires positive values. Use ' ...
                 'drop_pcs=[1] to remove the amplitude axis instead.'], ...
                 norm_cells, 100*neg_frac);
        end
        if norm_cells == "mean"
            cell_ref = mean(Xraw, 2, 'omitnan');
        else
            cell_ref = median(Xraw, 2, 'omitnan');
        end
        cell_ref(cell_ref <= 0 | isnan(cell_ref)) = NaN;   % guard div-by-0
        Xraw = Xraw ./ cell_ref;                            % broadcast per row
        if verbose
            nbad = sum(isnan(cell_ref));
            fprintf('[cycif_prep_features] per-cell %s normalization applied%s.\n', ...
                norm_cells, i_badcellstr(nbad));
        end
    otherwise
        error('cycif_prep_features:normCells', ...
            '''normalize_cells'' must be none | mean | median (got ''%s'').', norm_cells);
end

% ---- transform (metric only) ----
nz = ~isnan(Xraw);
Xt = Xraw;
switch transform
    case "arcsinh"
        assert(cofactor > 0, 'cycif_prep_features:cofactor', 'cofactor must be > 0.');
        Xt(nz) = asinh(Xraw(nz) ./ cofactor);
    case "log2"
        Xt(nz) = log2(max(Xraw(nz), 1));
    case "linear"
        % as-is
end

% ---- completeness / keep filter (BEFORE impute) ----
% Rationale (inherited from cycif_plot_embedding): low-completeness cells impute
% to near-identical rows, collapsing into an artifactual point that breaks the
% neighbour graph. Drop them here so only followable cells reach embed/cluster.
frac_measured = mean(nz, 2);
nIn = height(T);
if ~isempty(keep_cells)
    keep_row = i_resolve_keeplist(keep_cells, T);
else
    keep_row = frac_measured >= min_frac;
end
assert(sum(keep_row) >= 3, 'cycif_prep_features:tooFewCells', ...
    'Only %d cells pass the completeness/keep filter; lower min_marker_frac.', sum(keep_row));

if verbose
    fprintf('[cycif_prep_features] transform=%s%s | zscore=%d | n_pcs=%d\n', ...
        transform, i_cofstr(transform, cofactor), do_zscore, n_pcs);
    fprintf('[cycif_prep_features] completeness filter (>=%.0f%% of %d markers): %d/%d kept (%.1f%% dropped)\n', ...
        100*min_frac, nMk, sum(keep_row), nIn, 100*mean(~keep_row));
end

Tk = T(keep_row, :);
Xt = Xt(keep_row, :);

% ---- NaN imputation (on transformed scale) ----
switch impute_nan
    case "lowerquantile"
        for m = 1:nMk
            col = Xt(:, m); bad = isnan(col);
            if any(bad)
                fillv = quantile(col(~bad), 0.25);
                if isnan(fillv), fillv = 0; end
                col(bad) = fillv; Xt(:, m) = col;
            end
        end
    case "median"
        for m = 1:nMk
            col = Xt(:, m); bad = isnan(col);
            if any(bad)
                fillv = median(col(~bad), 'omitnan');
                if isnan(fillv), fillv = 0; end
                col(bad) = fillv; Xt(:, m) = col;
            end
        end
    case "none"
        good = all(~isnan(Xt), 2);
        if verbose && any(~good)
            fprintf('[cycif_prep_features] impute=none: dropping %d rows with residual NaN.\n', sum(~good));
        end
        Tk = Tk(good, :); Xt = Xt(good, :);
        keep_row(keep_row) = good;   % keep mask consistent with final rows
end

% ---- z-score (per marker) ----
if do_zscore
    mu  = mean(Xt, 1, 'omitnan');
    sig = std(Xt, 0, 1, 'omitnan'); sig(sig == 0) = 1;
    Xz  = (Xt - mu) ./ sig;
else
    mu = []; sig = [];
    Xz = Xt;
end
Xz(isnan(Xz)) = 0;   % safety: any residual NaN -> 0

% ---- opt-in marker correlation heatmap (on the features PCA sees) --------
% Shows the covariation structure that drives the PCs. A broadly warm
% (positive) matrix => a shared amplitude axis (uniform-sign PC1); structured
% blocks => real co-regulation. Diverging map, white at 0.
if plot_corr
    Rho = corr(Xz, 'rows','pairwise');
    figure('Color','w', 'Name','cycif_prep_features marker correlation');
    imagesc(Rho, [-1 1]); axis image; colorbar;
    n2 = 128;   % diverging blue-white-red (white=0), matches pipeline convention
    dm = [ [linspace(0.17,1,n2) ones(1,n2)]', ...
           [linspace(0.51,1,n2) linspace(1,0.30,n2)]', ...
           [linspace(0.73,1,n2) linspace(1,0.17,n2)]' ];
    colormap(gca, dm);
    set(gca, 'XTick',1:nMk, 'XTickLabel',markers, 'YTick',1:nMk, 'YTickLabel',markers, ...
             'TickLabelInterpreter','none'); xtickangle(90);
    title(sprintf('Marker correlation (arcsinh %.3g, z-scored)', cofactor));
end

% ---- PCA denoise (optional) ----
explained = []; cum_explained = []; coeff = [];
if n_pcs > 0
    n_keep = min([n_pcs, nMk, size(Xz,1)-1]);
    [coeff_full, score, ~, ~, expl] = pca(Xz);   % pca centers internally
    n_keep = min(n_keep, size(score,2));

    % --- Why dropping a PC (e.g. PC1) is legitimate, not data-mangling ------
    % PCs are an ORTHOGONAL basis: each captures a variance direction
    % independent of the others. Removing one leaves the remaining PCs exactly
    % unchanged (unlike dropping a raw marker, which would distort every
    % correlated marker). So excluding a nuisance PC is a clean projection onto
    % the complementary subspace, not a reweighting of the rest.
    %
    % Dropping PC1 here is justified specifically because PC1 was diagnosed as a
    % TECHNICAL AMPLITUDE axis, not biology:
    %   - 13/13 uniform-sign loadings of near-equal magnitude => every marker
    %     rises/falls together = a global brightness/size/staining-depth scale,
    %     not a polarization contrast (which would show M1-vs-M2 sign structure).
    %   - It captured ~55% of variance precisely because that scale is shared
    %     across all channels, swamping the sparse per-marker signal.
    %   - The data is background-SUBTRACTED (additive offset removed in
    %     ct_cycif_measure), so this amplitude confound is additive/correlated
    %     and survives subtraction; it is NOT the multiplicative confound that
    %     per-cell division fixes -- hence drop-the-PC, not normalize.
    % Projecting it out lets the neighbour graph / clustering see relative marker
    % PROFILE (what a cell is) instead of overall intensity (how bright it is).
    % Caveat: confirm via loadings before dropping -- only drop a PC shown to be
    % a nuisance axis; do not drop by rank reflexively.
    keep_pc = setdiff(1:n_keep, drop_pcs, 'stable');   % PCs to retain in X
    dropped = intersect(1:n_keep, drop_pcs);
    X = score(:, keep_pc);
    coeff = coeff_full(:, keep_pc);
    explained = expl(keep_pc);
    cum_explained = sum(explained);
    if verbose
        fprintf('[cycif_prep_features] PCA: kept %d of first %d PCs (%.1f%% variance)', ...
            numel(keep_pc), n_keep, cum_explained);
        if ~isempty(dropped)
            fprintf(' | DROPPED PC%s (%.1f%% variance removed)', ...
                strjoin(string(dropped),','), sum(expl(dropped)));
        end
        fprintf('\n');
    end
    n_pcs_out = numel(keep_pc);

    % Opt-in PCA selection plot: full scree with the retained-PC cut marked.
    % Uses the FULL explained vector (all PCs) so the flattening point is visible.
    if plot_scree
        cum_all = cumsum(expl);
        figure('Color','w', 'Name','cycif_prep_features PCA scree');
        yyaxis left
        bar(expl, 'FaceColor',[0.7 0.8 0.9], 'EdgeColor','none'); ylabel('% variance / PC');
        yyaxis right
        plot(cum_all, '-o', 'LineWidth',1.5, 'MarkerSize',4); ylabel('cumulative % variance');
        ylim([0 100]); yline(90, '--', '90%');
        xline(n_keep + 0.5, 'r-', sprintf('n\\_pcs=%d (%.1f%%)', n_keep, cum_all(n_keep)), ...
              'LabelVerticalAlignment','bottom', 'LineWidth',1.2);
        xlabel('principal component');
        title(sprintf('PCA scree: arcsinh(%.3g), z-score=%d', cofactor, do_zscore));
        grid on;
    end
else
    X = Xz;
    n_pcs_out = 0;
    if plot_scree
        warning('cycif_prep_features:noScree', ...
            'plot_scree requested but n_pcs=0 (PCA off) -- nothing to plot.');
    end
end

% ---- identity keys ----
key_vars = {'xy','cellid'};
if ismember('uid', Tk.Properties.VariableNames), key_vars = [key_vars, {'uid'}]; end
assert(all(ismember({'xy','cellid'}, Tk.Properties.VariableNames)), ...
    'cycif_prep_features:keys', 'T must contain xy and cellid to build cell_idx.');
cell_idx = Tk(:, key_vars);

% ---- info ----
info = struct( ...
    'markers',       {markers}, ...
    'transform',     char(transform), ...
    'normalize_cells', char(norm_cells), ...
    'cofactor',      i_cofval(transform, cofactor), ...
    'zscore',        do_zscore, ...
    'n_pcs',         n_pcs_out, ...
    'dropped_pcs',   drop_pcs, ...
    'explained',     explained, ...
    'cum_explained', cum_explained, ...
    'keep_row',      keep_row, ...
    'mu',            mu, ...
    'sig',           sig, ...
    'coeff',         coeff);
end


% ====================================================================
function markers = i_auto_markers(T)
% Numeric, non-ID columns. Matches cycif_plot_embedding's convention.
cols    = T.Properties.VariableNames;
is_num  = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
is_id   = ismember(cols, {'xy','cellid','uid','ptx_is_spike'});
markers = cols(is_num & ~is_id);
end

function keep = i_resolve_keeplist(keep_cells, T)
% logical mask | numeric idx | xy+cellid table -> logical over rows of T.
n = height(T);
if islogical(keep_cells)
    assert(numel(keep_cells)==n, 'keep_cells logical length (%d) must equal height(T) (%d).', numel(keep_cells), n);
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

function s = i_badcellstr(nbad)
if nbad > 0, s = sprintf(' (%d cells with non-positive reference left as NaN)', nbad); else, s = ''; end
end

function s = i_cofstr(transform, cofactor)
if transform == "arcsinh", s = sprintf(' (cofactor=%.4g)', cofactor); else, s = ''; end
end

function v = i_cofval(transform, cofactor)
if transform == "arcsinh", v = cofactor; else, v = NaN; end
end