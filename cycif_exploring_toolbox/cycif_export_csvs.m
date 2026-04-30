function out_dir = cycif_export_csvs(T, out_dir, varargin)
% cycif_export_csvs  Export CyCIF table and optional embeddings to CSVs
%                    for downstream conversion to AnnData (.h5ad).
%
% INPUT
%   T        — table from cycif_build_table
%   out_dir  — folder to write CSVs into (created if missing)
%
% NAME-VALUE OPTIONS
%   'use_log2'    logical    log2(max(x,1)) transform markers      (default true)
%   'impute_nan'  string     NaN fill strategy before export:
%                            'none'          — keep NaN (default)
%                            'lowerquantile' — global 25th pct per marker
%                            'median'        — global median per marker
%   'embeddings'  struct     pre-computed embeddings to export.
%                            Fields are method names, values are [nCells x 2].
%                            e.g. struct('tsne', emb_tsne, 'pca', emb_pca)
%                            Pass [] to skip embedding export.
%
% OUTPUT
%   out_dir  — path to folder containing:
%                X.csv          cells x markers intensity matrix
%                obs.csv        cell metadata (conditions, xy, cellid)
%                var.csv        marker names and index
%                emb_<name>.csv one file per embedding in 'embeddings'
%
% USAGE
%   out_dir = cycif_export_csvs(T, 'C:\path\to\export');
%   out_dir = cycif_export_csvs(T, 'C:\path\to\export', ...
%               'embeddings', struct('tsne', emb_tsne, 'pca', emb_pca));
%   out_dir = cycif_export_csvs(T, 'C:\path\to\export', ...
%               'impute_nan','lowerquantile', 'embeddings', struct('tsne', emb_tsne));

p = inputParser();
addParameter(p, 'use_log2',   true);
addParameter(p, 'impute_nan', 'none');
addParameter(p, 'embeddings', []);
parse(p, varargin{:});

use_log2   = p.Results.use_log2;
impute_nan = lower(p.Results.impute_nan);
embeddings = p.Results.embeddings;

% --- Create output directory ----------------------------------------------
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fprintf('[cycif_export_csvs] Writing to: %s\n', out_dir);

% --- Identify marker columns (numeric, not identity/flag) -----------------
cols   = T.Properties.VariableNames;
is_num = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
is_id  = ismember(cols, {'xy','cellid','ptx_is_spike'});
markers = cols(is_num & ~is_id);
nMk     = numel(markers);
fprintf('[cycif_export_csvs] %d markers: %s\n', nMk, strjoin(markers, ', '));

% --- X matrix: cells x markers -------------------------------------------
Xmat = table2array(T(:, markers));
if use_log2
    nz = ~isnan(Xmat);
    Xmat(nz) = log2(max(Xmat(nz), 1));
    fprintf('[cycif_export_csvs] Applied log2 transform (NaN preserved).\n');
end

% --- NaN imputation ------------------------------------------------------
switch impute_nan
    case 'none'
        fprintf('[cycif_export_csvs] NaN preserved in X.csv (impute_nan=none).\n');
    case 'lowerquantile'
        n_imp = 0;
        for m = 1:nMk
            col = Xmat(:,m); nan_idx = isnan(col);
            if ~any(nan_idx), continue; end
            q25 = quantile(col(~nan_idx), 0.25);
            if isnan(q25), q25 = 0; end
            col(nan_idx) = q25; n_imp = n_imp + sum(nan_idx);
            Xmat(:,m) = col;
        end
        fprintf('[cycif_export_csvs] Imputed %d NaN with per-marker 25th percentile.\n', n_imp);
    case 'median'
        n_imp = 0;
        for m = 1:nMk
            col = Xmat(:,m); nan_idx = isnan(col);
            if ~any(nan_idx), continue; end
            med = median(col(~nan_idx), 'omitnan');
            if isnan(med), med = 0; end
            col(nan_idx) = med; n_imp = n_imp + sum(nan_idx);
            Xmat(:,m) = col;
        end
        fprintf('[cycif_export_csvs] Imputed %d NaN with per-marker global median.\n', n_imp);
    otherwise
        error('[cycif_export_csvs] Unknown impute_nan ''%s''. Use ''none'', ''lowerquantile'', or ''median''.', impute_nan);
end

% Write with marker names as header
Xmat_table = array2table(Xmat, 'VariableNames', markers);
writetable(Xmat_table, fullfile(out_dir, 'X.csv'));
fprintf('[cycif_export_csvs] Wrote X.csv [%d cells x %d markers]\n', size(Xmat,1), nMk);

% --- obs: cell metadata ---------------------------------------------------
meta_cols = cols(~is_num | is_id);
% Always include xy and cellid first
priority  = {'xy','cellid'};
rest      = setdiff(meta_cols, priority, 'stable');
obs_cols  = [priority, rest];
obs_cols  = obs_cols(ismember(obs_cols, cols));

Tobs = T(:, obs_cols);
% Convert ALL columns to string for clean, unambiguous CSV export.
% categorical, logical, and numeric columns are all cast to string.
% This avoids writetable silently mangling categoricals in some MATLAB versions.
for c = 1:width(Tobs)
    vname = Tobs.Properties.VariableNames{c};
    col   = Tobs.(vname);
    if iscategorical(col)
        Tobs.(vname) = string(col);      % categorical -> string
    elseif isnumeric(col) || islogical(col)
        Tobs.(vname) = string(col);      % numeric/logical -> string
    elseif ~isstring(col) && ischar(col)
        Tobs.(vname) = string(col);      % char -> string
    end
    % string columns pass through unchanged
end
writetable(Tobs, fullfile(out_dir, 'obs.csv'));
fprintf('[cycif_export_csvs] Wrote obs.csv [%d columns]: %s\n', ...
        width(Tobs), strjoin(Tobs.Properties.VariableNames, ', '));

% --- var: marker index ----------------------------------------------------
var_table = table(markers(:), 'VariableNames', {'marker_name'});
writetable(var_table, fullfile(out_dir, 'var.csv'));
fprintf('[cycif_export_csvs] Wrote var.csv\n');

% --- embeddings -----------------------------------------------------------
if ~isempty(embeddings)
    emb_names = fieldnames(embeddings);
    for ei = 1:numel(emb_names)
        name = emb_names{ei};
        emb  = embeddings.(name);
        assert(size(emb,1) == height(T), ...
            'Embedding ''%s'' has %d rows but T has %d — must match.', ...
            name, size(emb,1), height(T));
        col1 = sprintf('%s_1', name);
        col2 = sprintf('%s_2', name);
        emb_table = array2table(emb(:,1:2), 'VariableNames', {col1, col2});
        fname = fullfile(out_dir, sprintf('emb_%s.csv', name));
        writetable(emb_table, fname);
        fprintf('[cycif_export_csvs] Wrote emb_%s.csv\n', name);
    end
end

fprintf('[cycif_export_csvs] Done.\n');
end