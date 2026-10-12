function [labels, info] = cycif_leiden(X, varargin)
% cycif_leiden  Leiden community detection via the Python bridge (leidenalg +
%               igraph). Fast C backend; strictly improves on Louvain
%               (guarantees well-connected communities).
%
%   Builds the SAME k-NN graph as cycif_louvain (via cycif_knn_graph), then
%   runs Leiden with RBConfigurationVertexPartition so the 'gamma' resolution
%   knob matches the Louvain interface. Any difference vs Louvain is therefore
%   the algorithm, not the graph.
%
% Usage examples:
%   labels = cycif_leiden(X);                              % 15-NN, gamma 1.0
%   labels = cycif_leiden(X, 'graph_k',15, 'gamma',1.0, 'seed',0);
%   [labels, info] = cycif_leiden(X, 'gamma',1.4, 'python_handler',uh);
%
% INPUT
%   X  [nObs x nFeat] double, finite. Feed cycif_prep_features output.
%
% NAME-VALUE
%   'graph_k'   scalar   kNN neighbors                          (default 15)
%   'gamma'     scalar   resolution (>1 finer/more communities) (default 1.0)
%   'weighting' string   'gaussian' | 'binary'                  (default 'gaussian')
%   'mutual'    logical  mutual-kNN edges only                  (default false)
%   'n_iter'    scalar   Leiden iterations (-1 = to convergence) (default -1)
%   'seed'      scalar   RNG seed                               (default 0)
%   'python_handler' PythonHandler  instance to prepare the env (default []
%                       constructs a bare PythonHandler; reuse your UmapHandler
%                       to share the already-prepared interpreter)
%   'verbose'   logical                                         (default true)
%
% OUTPUT
%   labels  [nObs x 1] double community IDs 1..K, size-ranked (1 = largest).
%   info    struct: .k, .gamma, .n_communities, .modularity (Leiden quality),
%                   .sizes, .graph (from cycif_knn_graph)
%
% REQUIRES  leidenalg + igraph (python-igraph) in the active Python env:
%   "<python>" -m pip install leidenalg igraph
%
% See also: cycif_knn_graph, cycif_louvain, cycif_prep_features, UmapHandler.

p = inputParser();
addRequired(p, 'X', @(v) isnumeric(v) && ismatrix(v) && ~isempty(v));
addParameter(p, 'graph_k',   15);
addParameter(p, 'gamma',     1.0);
addParameter(p, 'weighting', 'gaussian');
addParameter(p, 'mutual',    false);
addParameter(p, 'n_iter',    -1);
addParameter(p, 'seed',      0);
addParameter(p, 'python_handler', []);
addParameter(p, 'verbose',   true);
parse(p, X, varargin{:});
k        = round(p.Results.graph_k);
gamma    = p.Results.gamma;
weighting= char(p.Results.weighting);
mutual   = logical(p.Results.mutual);
n_iter   = round(p.Results.n_iter);
seed     = round(p.Results.seed);
ph       = p.Results.python_handler;
verbose  = logical(p.Results.verbose);

X = double(X);
n = size(X,1);

% --- 1. Shared graph (identical to cycif_louvain) -----------------------
[A, gInfo] = cycif_knn_graph(X, 'graph_k', k, 'weighting', weighting, ...
                             'mutual', mutual, 'verbose', verbose);

% --- 2. Prepare Python env ----------------------------------------------
% Reuse a passed handler (e.g. your UmapHandler) so the interpreter is shared/
% already prepared; else construct a bare PythonHandler to resolve one.
if isempty(ph)
    ph = PythonHandler;
end
ph.preparePython();
ph.checkPython();
% Verify leidenalg present (clear message + install command on failure).
try
    py.importlib.import_module('leidenalg');
    py.importlib.import_module('igraph');
catch ME
    error('cycif_leiden:noLeiden', ...
        ['leidenalg/igraph not importable in the active Python env.\n' ...
         'Install:  "%s" -m pip install leidenalg igraph\nError: %s'], ...
         ph.getPythonExe(), ME.message);
end

% --- 3. Edge list (upper triangle; igraph is undirected) ----------------
% Extract i<j edges + weights from the symmetric sparse A.
[si, sj, sw] = find(triu(A, 1));
if verbose
    fprintf('[cycif_leiden] passing %d edges to igraph/leidenalg (gamma=%.3g)...\n', ...
        numel(si), gamma);
end

% --- 4. Run Leiden in Python --------------------------------------------
% Pass 0-based edges (Python) and weights; build igraph, run leidenalg with
% RBConfigurationVertexPartition (resolution_parameter = gamma). Return the
% membership vector (0-based) and the partition quality.
code = [ ...
    "import igraph as ig, leidenalg as la, numpy as np"                    , ...
    "ei = np.asarray(si, dtype=np.int64).ravel() - 1"                      , ... % ->0-based
    "ej = np.asarray(sj, dtype=np.int64).ravel() - 1"                      , ...
    "wt = np.asarray(sw, dtype=np.float64).ravel()"                        , ...
    "edges = list(zip(ei.tolist(), ej.tolist()))"                          , ...
    "g = ig.Graph(n=int(nn), edges=edges, edge_attrs={'weight': wt.tolist()})", ...
    "part = la.find_partition(g, la.RBConfigurationVertexPartition,"        , ...
    "        weights='weight', resolution_parameter=float(gamma),"          , ...
    "        n_iterations=int(niter), seed=int(seed))"                       , ...
    "membership = np.asarray(part.membership, dtype=np.float64)"           , ...
    "quality = float(part.quality())"                                      ];
[pyMemb, pyQual] = pyrun(code, ["membership","quality"], ...
    si = si, sj = sj, sw = sw, nn = n, gamma = gamma, niter = n_iter, seed = seed);

try
    memb = double(pyMemb);    % 0-based community per node
catch
    memb = double(pyMemb.tolist());
end
Q    = double(pyQual);
memb = memb(:);
assert(numel(memb) == n, 'cycif_leiden:size', ...
    'Leiden returned %d labels for %d cells.', numel(memb), n);

% --- 5. Size-rank relabel (1 = largest), matching cycif_louvain ----------
[labels, sizes] = i_rank_labels(memb + 1);   % +1 -> 1-based before ranking

info = struct('k', k, 'gamma', gamma, 'n_communities', numel(sizes), ...
    'modularity', Q, 'sizes', sizes, 'graph', gInfo);
if verbose
    fprintf('[cycif_leiden] %d communities, quality=%.4f. Sizes: %s\n', ...
        numel(sizes), Q, mat2str(sizes(:)'));
end
end


% ======================================================================
function [labels, sizes] = i_rank_labels(c)
% Relabel so community 1 is largest, 2 next, etc.
u = unique(c);
cnt = arrayfun(@(x) sum(c==x), u);
[sizes, ord] = sort(cnt, 'descend');
remap = zeros(max(u),1);
remap(u(ord)) = 1:numel(u);
labels = remap(c);
end