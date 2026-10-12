function cycif_to_h5ad(T, out_path, varargin)
% cycif_to_h5ad  Export CyCIF table to cellxgene-compatible .h5ad file.
%                Orchestrates cycif_export_csvs.m + cycif_to_h5ad.py.
%
% INPUT
%   T         — table from cycif_build_table
%   out_path  — full path for output .h5ad file
%
% NAME-VALUE OPTIONS
%   'use_log2'     logical   log2(max(x,1)) transform            (default false)
%   'impute_nan'   string    NaN fill strategy: 'none'|'lowerquantile'|'median'
%                            (default 'none') — passed to cycif_export_csvs
%   'embeddings'   struct    pre-computed COORDINATE embeddings [nCells x 2].
%                            e.g. struct('tsne', emb_tsne, 'pca', emb_pca)
%                            Land in adata.obsm as 'X_tsne', 'X_pca'.
%   'annotations'  struct    per-cell LABEL vectors [nCells x 1] (cluster ids,
%                            leiden, etc). Each field becomes an obs column and
%                            (if named in the Python cat list) a cellxgene
%                            categorical filter. e.g. struct('leiden', emb)
%   'pyexe'        string    path to Python executable
%                            (default: cellpose-matlab env)
%   'keep_csvs'    logical   keep intermediate CSV folder         (default false)
%
% USAGE
%   % Basic
%   cycif_to_h5ad(T, 'C:\path\to\output\cycif.h5ad');
%
%   % With coordinate embeddings + leiden cluster annotation
%   cycif_to_h5ad(T, 'C:\path\to\output\cycif.h5ad', ...
%                 'embeddings',  struct('tsne', emb_tsne, 'pca', emb_pca), ...
%                 'annotations', struct('leiden', leiden_labels), ...
%                 'impute_nan',  'lowerquantile');
%
% NOTE
%   embeddings = 2-D coordinates (go to obsm). annotations = 1-D labels
%   (go to obs). Do not pass cluster labels through 'embeddings' — that path
%   forces exactly 2 columns and will error or misread a label vector.

p = inputParser();
addParameter(p, 'use_log2',    false);
addParameter(p, 'impute_nan',  'none');
addParameter(p, 'embeddings',  []);
addParameter(p, 'annotations', []);
addParameter(p, 'pyexe',       'C:\Users\mhardy\AppData\Local\anaconda3\envs\cellpose-matlab\python.exe');
addParameter(p, 'keep_csvs',   false);
parse(p, varargin{:});

pyexe     = p.Results.pyexe;
keep_csvs = p.Results.keep_csvs;

% Locate the Python script alongside this .m file (or on path)
py_script = which('cycif_to_h5ad.py');
assert(~isempty(py_script), ...
    'cycif_to_h5ad.py not found on MATLAB path — add its folder with addpath().');

% --- Step 0: Merge annotations into T as obs columns ----------------------
% Label vectors ride out as ordinary obs columns; the Python side casts
% known category names (e.g. 'cluster','leiden') to categorical for cellxgene.
annot = p.Results.annotations;
if ~isempty(annot)
    afields = fieldnames(annot);
    for ai = 1:numel(afields)
        aname = afields{ai};
        avec  = annot.(aname);
        assert(numel(avec) == height(T), ...
            'Annotation ''%s'' has %d entries but T has %d rows — must match.', ...
            aname, numel(avec), height(T));
        % Store as string so it exports as a clean categorical, never as a
        % marker (cycif_export_csvs routes numeric non-id cols into X).
        T.(aname) = string(avec(:));
        fprintf('[cycif_to_h5ad] Added annotation ''%s'' to obs (%d cells).\n', ...
                aname, numel(avec));
    end
end

% --- Step 1: Export CSVs to a temp folder ---------------------------------
csv_dir = [tempname '_cycif_csvs'];
cycif_export_csvs(T, csv_dir, ...
    'use_log2',   p.Results.use_log2, ...
    'impute_nan', p.Results.impute_nan, ...
    'embeddings', p.Results.embeddings);

% --- Step 1b: Set OpenMP env vars to avoid MATLAB/Python DLL conflict ------
% MATLAB already loads libiomp5md.dll; numpy loads it again on import.
% KMP_DUPLICATE_LIB_OK allows both to coexist (same workaround as celltracer_v3).
setenv('KMP_DUPLICATE_LIB_OK',  'TRUE');
setenv('OMP_NUM_THREADS',        '1');
setenv('MKL_NUM_THREADS',        '1');
setenv('OPENBLAS_NUM_THREADS',   '1');

% --- Step 2: Check anndata is available in the Python env -----------------
[st, out] = system(sprintf('"%s" -c "import anndata; print(anndata.__version__)"', pyexe));
if st ~= 0
    fprintf('[cycif_to_h5ad] anndata not found — installing...\n');
    system(sprintf('"%s" -m pip install anndata', pyexe));
end

% --- Step 3: Call Python conversion script --------------------------------
cmd = sprintf('"%s" "%s" "%s" "%s"', pyexe, py_script, csv_dir, out_path);
fprintf('[cycif_to_h5ad] Running: %s\n', cmd);
[st, out] = system(cmd);
disp(out);
assert(st == 0, '[cycif_to_h5ad] Python conversion failed — see output above.');

fprintf('[cycif_to_h5ad] h5ad written to: %s\n', out_path);

% --- Step 4: Clean up temp CSVs -------------------------------------------
if ~keep_csvs
    rmdir(csv_dir, 's');
    fprintf('[cycif_to_h5ad] Temp CSVs removed.\n');
else
    fprintf('[cycif_to_h5ad] Temp CSVs kept at: %s\n', csv_dir);
end
end