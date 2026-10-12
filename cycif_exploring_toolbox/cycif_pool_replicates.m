function [T, info] = cycif_pool_replicates(Tlist, varargin)
% CYCIF_POOL_REPLICATES
% Pool per-cell tables from multiple experimental replicates into one table,
% with per-replicate per-marker scaling to remove run-to-run staining batch
% effects BEFORE a shared embedding / clustering / polarization score.
%
% WHY THIS STAGE EXISTS
%   cycif_build_table stores RAW intensities and adds no experiment/replicate
%   provenance. (xy, cellid) is unique only WITHIN one experiment -- XY28 in
%   rep1 and XY28 in rep3 are different cells with the same key. Naive
%   [T1;T2;T3] therefore (a) collides those keys, breaking any later
%   (xy,cellid) join (e.g. emb_cell_idx back into T), and (b) pools raw values
%   that are not comparable across runs when staining drifts run-to-run. This
%   function fixes both: it adds a 'rep' tag + globally-unique 'uid', and
%   z-scores each marker WITHIN each replicate so a per-run intensity shift/
%   scale does not become artefactual structure in the shared feature space.
%
%   Scaling is deliberately a DISTINCT stage between cycif_build_table and the
%   embedding. Because the marker values are now pre-scaled, the embedding
%   must NOT scale them again -- call cycif_plot_embedding with
%   'use_log2',false and 'zscore_embed',false on the pooled table.
%
% METHOD (per this run's chosen configuration)
%   Per replicate, per marker column: z = (x - mean) / std, computed on RAW
%   LINEAR values (no log), ignoring NaN. This standardizes each replicate's
%   marker distribution to mean 0 / sd 1, removing a per-run location+scale
%   shift. NOTE: mean/std on heavy-tailed raw intensities are tail-sensitive;
%   inspect the replicate-colored embedding to confirm reps intermix (i.e. the
%   batch effect was a location/scale shift this removes) rather than forming
%   separate islands (which would indicate a distribution-SHAPE difference that
%   needs quantile normalization or a dedicated integration method instead).
%
% INPUTS
%   Tlist   cell array {T1, T2, T3, ...} of per-replicate tables from
%           cycif_build_table. All must share the SAME marker columns.
%
% PAIRED OPTIONS
%   'rep_labels'  cellstr  one label per table (default {'rep1','rep2',...}).
%   'markers'     cellstr  marker columns to scale ([] = auto: all columns
%                          ending in _nuc/_cyt/_cell). Metadata columns are
%                          never scaled.                          (default [])
%   'scale'       logical  apply per-replicate z-score (default true). false =
%                          pool + tag only, no scaling (raw pooled table).
%   'verbose'     logical  print per-replicate scaling stats     (default true)
%
% OUTPUTS
%   T     pooled table: all input rows, in input order, with two added columns
%         'rep'  (string, the replicate label) and 'uid' (string,
%         rep_xy_cellid -- globally unique). Marker columns are z-scored
%         per replicate when scale=true.
%   info  struct: .rep_labels, .n_per_rep, .markers_scaled,
%                 .mu [nRep x nMk], .sig [nRep x nMk]  (the scaling constants
%                 actually applied, for audit / back-transform).
%
% DOWNSTREAM USAGE
%   T = cycif_pool_replicates({T1,T2,T3});
%   % shared embedding on the pre-scaled pool -- disable the embedding's own
%   % scaling so values are not scaled twice:
%   [fig, emb, idx] = cycif_plot_embedding(T, 'method','tsne', ...
%        'use_log2',false, 'zscore_embed',false, ...
%        'group_by','rep');                 % <- color by rep to judge batch
%   [T, labels] = cycif_cluster(T, emb, 'emb_cell_idx', idx, 'method','dbscan');

% ---- parse ----
p = struct('rep_labels', [], 'markers', [], 'scale', true, 'verbose', true);
for k = 1:2:numel(varargin)
    key = lower(varargin{k});
    if ~isfield(p, key); error('cycif_pool_replicates: unknown option ''%s''', varargin{k}); end
    p.(key) = varargin{k+1};
end
assert(iscell(Tlist) && ~isempty(Tlist), ...
    'cycif_pool_replicates: Tlist must be a non-empty cell array of tables.');
nRep = numel(Tlist);

% ---- replicate labels ----
if isempty(p.rep_labels)
    rep_labels = arrayfun(@(i) sprintf('rep%d', i), 1:nRep, 'uni', 0);
else
    rep_labels = cellstr(p.rep_labels);
    assert(numel(rep_labels) == nRep, ...
        'cycif_pool_replicates: rep_labels has %d entries but %d tables given.', ...
        numel(rep_labels), nRep);
end

% ---- verify identical variable sets across replicates ----
vars1 = Tlist{1}.Properties.VariableNames;
for i = 2:nRep
    vi = Tlist{i}.Properties.VariableNames;
    if ~isequal(sort(vars1), sort(vi))
        only1 = setdiff(vars1, vi);
        only2 = setdiff(vi, vars1);
        error(['cycif_pool_replicates: replicate %d has a different column set than replicate 1.\n' ...
               '  only in rep1: %s\n  only in rep%d: %s\n' ...
               'Reconcile the marker/metadata columns before pooling (a marker excluded in\n' ...
               'one replicate but not others is the usual cause).'], ...
               i, strjoin(only1, ', '), i, strjoin(only2, ', '));
    end
end

% ---- resolve marker columns to scale ----
if isempty(p.markers)
    markers = vars1(endsWith(vars1, {'_nuc','_cyt','_cell'}));
else
    markers = cellstr(p.markers);
    missing = markers(~ismember(markers, vars1));
    assert(isempty(missing), ...
        'cycif_pool_replicates: requested marker(s) not in table: %s', strjoin(missing, ', '));
end
nMk = numel(markers);
assert(nMk >= 1, 'cycif_pool_replicates: no marker columns found to scale.');

% ---- per-replicate scaling + tagging ----
mu_all  = nan(nRep, nMk);
sig_all = nan(nRep, nMk);
n_per   = zeros(1, nRep);
out = cell(nRep, 1);

for i = 1:nRep
    Ti = Tlist{i};
    n_per(i) = height(Ti);

    if p.scale
        for m = 1:nMk
            col = Ti.(markers{m});
            mu  = mean(col, 'omitnan');
            sg  = std(col, 'omitnan');
            if isnan(sg) || sg == 0; sg = 1; end   % degenerate: leave centered only
            mu_all(i,m)  = mu;
            sig_all(i,m) = sg;
            Ti.(markers{m}) = (col - mu) ./ sg;
        end
    end

    % provenance + globally-unique key. (xy,cellid) alone repeats across reps.
    Ti.rep = repmat(string(rep_labels{i}), height(Ti), 1);
    Ti.uid = Ti.rep + "_" + string(Ti.xy) + "_" + string(Ti.cellid);
    out{i} = Ti;

    if p.verbose
        if p.scale
            fprintf('[cycif_pool_replicates] %s: %d cells, %d markers z-scored (raw linear).\n', ...
                    rep_labels{i}, height(Ti), nMk);
        else
            fprintf('[cycif_pool_replicates] %s: %d cells, no scaling (scale=false).\n', ...
                    rep_labels{i}, height(Ti));
        end
    end
end

% ---- concatenate ----
T = vertcat(out{:});

% uid uniqueness guard: if it is NOT unique, a downstream (rep,xy,cellid) or
% uid join would silently mis-map. Surface it here rather than let it corrupt
% cluster label placement.
if numel(unique(T.uid)) ~= height(T)
    n_dup = height(T) - numel(unique(T.uid));
    warning(['cycif_pool_replicates: %d duplicate uid(s) after pooling. (xy,cellid) is not ' ...
             'unique within a replicate for this data -- downstream joins on uid may mis-map. ' ...
             'Check that cellid is unique per (rep,xy).'], n_dup);
end

info = struct('rep_labels', {rep_labels}, 'n_per_rep', n_per, ...
              'markers_scaled', {markers}, 'mu', mu_all, 'sig', sig_all, ...
              'scaled', p.scale);

if p.verbose
    fprintf('[cycif_pool_replicates] Pooled %d replicates -> %d cells x %d markers.\n', ...
            nRep, height(T), nMk);
    if p.scale
        fprintf('  Marker values are per-replicate z-scored: call cycif_plot_embedding with\n');
        fprintf('  ''use_log2'',false and ''zscore_embed'',false to avoid double-scaling.\n');
        fprintf('  Color the embedding by ''rep'' to confirm replicates intermix (batch removed).\n');
    end
end
end