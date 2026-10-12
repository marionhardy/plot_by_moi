%UmapHandler
%   PythonHandler subclass exposing UMAP dimensionality reduction to MATLAB.
%
%   Built against the deployed PythonHandler (the MatCellPose stub), whose
%   surface is preparePython / checkPython / getPythonExe with NO global
%   lock and NO abstract contract. This class therefore owns readiness
%   testing (umap import) and the pyrun marshalling itself.
%
% Usage examples:
%   uh  = UmapHandler;                 % uses pyenv / PATH python by default
%   uh.PythonExe = "C:\Python314\python.exe";   % or pin an interpreter
%   uh.ready                           % true if env runs AND umap imports
%   % --- typical: embed the prepped feature matrix ---
%   [X, cell_idx] = cycif_prep_features(T, 'transform','arcsinh', 'n_pcs',7, 'drop_pcs',[1]);
%   emb = uh.embed(X, 'n_neighbors',15, 'min_dist',0.1, 'seed',0);
%
% Requires umap-learn in the selected Python env. On Python 3.14 this needs
% numba>=0.63 / llvmlite>=0.46 (cp314 wheels, available on PyPI). Install once
% at a system shell:  <python> -m pip install umap-learn
%
% See also: PythonHandler, cycif_prep_features, cycif_plot_embedding.

classdef UmapHandler < PythonHandler

    properties
        % UMAP import status, cached after a successful ready() check.
        umapReady (1,1) logical = false;
    end

    methods
        function uh = UmapHandler(varargin)
            uh@PythonHandler(varargin{:});
        end

        % ready  Verify the Python env runs and umap-learn imports.
        %   tf = uh.ready()  -> logical; caches result in uh.umapReady.
        function tf = ready(uh)
            uh.preparePython();     % pick interpreter (env var / pyenv / PATH)
            uh.checkPython();       % errors if the interpreter won't run
            try
                py.importlib.import_module('umap');
                uh.umapReady = true;
            catch ME
                uh.umapReady = false;
                warning('UmapHandler:noUmap', ...
                    ['umap-learn not importable in the active Python env.\n' ...
                     'Install with:  "%s" -m pip install umap-learn\nError: %s'], ...
                    uh.getPythonExe(), ME.message);
            end
            tf = uh.umapReady;
        end

        % embed  Run UMAP on a feature matrix -> 2-D (or n-D) coordinates.
        %
        %   emb = uh.embed(X, Name,Value...)
        %
        % INPUT
        %   X   [nObs x nFeat] double, finite. Feed cycif_prep_features output
        %       so embedding + clustering share identical features.
        % NAME-VALUE (passed to umap.UMAP)
        %   'n_neighbors'  (default 15)   'min_dist'     (default 0.1)
        %   'n_components' (default 2)    'metric'       (default 'euclidean')
        %   'seed'         (default 0)    random_state for reproducibility
        % OUTPUT
        %   emb  [nObs x n_components] double, row-aligned to X.
        function emb = embed(uh, X, varargin)
            q = inputParser();
            addRequired(q, 'X', @(v) isnumeric(v) && ismatrix(v) && ~isempty(v));
            addParameter(q, 'n_neighbors', 15);
            addParameter(q, 'min_dist',    0.1);
            addParameter(q, 'n_components',2);
            addParameter(q, 'metric',      'euclidean');
            addParameter(q, 'seed',        0);
            parse(q, X, varargin{:});
            nn   = double(q.Results.n_neighbors);
            md   = double(q.Results.min_dist);
            ncmp = double(q.Results.n_components);
            met  = char(q.Results.metric);
            seed = double(q.Results.seed);

            % --- validate in MATLAB (fail early, clear message) ---
            X = double(X);
            assert(all(isfinite(X(:))), 'UmapHandler:nonFinite', ...
                'X contains NaN/Inf. cycif_prep_features imputes these before embedding.');
            nObs = size(X,1);
            assert(nObs > nn, 'UmapHandler:tooFewObs', ...
                'n_neighbors (%d) must be < number of observations (%d).', nn, nObs);

            % --- ensure umap is ready ---
            if ~uh.umapReady
                assert(uh.ready(), 'UmapHandler:notReady', ...
                    'umap-learn is not ready; see the warning above for the install command.');
            end

            % --- run UMAP through pyrun; marshal X in, embedding out ---
            code = [ ...
                "import numpy as np, umap"                                      , ...
                "reducer = umap.UMAP(n_neighbors=int(nn), min_dist=float(md)," , ...
                "                    n_components=int(ncmp), metric=str(met)," , ...
                "                    random_state=int(seed))"                    , ...
                "emb = reducer.fit_transform(np.asarray(X, dtype=np.float64))" , ...
                "emb = np.ascontiguousarray(emb, dtype=np.float64)"           ];
            pyEmb = pyrun(code, "emb", ...
                X = X, nn = nn, md = md, ncmp = ncmp, met = met, seed = seed);

            try
                emb = double(pyEmb);
            catch
                emb = double(pyEmb.tolist());
            end
            emb = reshape(emb, nObs, ncmp);
        end
    end
end