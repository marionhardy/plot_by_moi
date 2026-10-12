function [A, info] = cycif_knn_graph(X, varargin)
% cycif_knn_graph  Build a symmetric weighted k-NN similarity graph.
%
%   Single source of truth for the graph used by BOTH cycif_louvain and
%   cycif_leiden, so the two community-detection methods operate on the
%   identical graph and any difference in results is the ALGORITHM, not the
%   graph construction.
%
% Usage examples:
%   A = cycif_knn_graph(X);                            % 15-NN, gaussian, union
%   A = cycif_knn_graph(X, 'graph_k',10, 'mutual',true);
%   [A, info] = cycif_knn_graph(X, 'weighting','binary');
%
% INPUT
%   X  [nObs x nFeat] double, finite (cells x features).
%
% NAME-VALUE
%   'graph_k'   scalar  neighbors per node                    (default 15)
%   'weighting' string  'gaussian' | 'binary'                 (default 'gaussian')
%                       gaussian: exp(-d^2/sigma^2), sigma = median kNN dist
%   'mutual'    logical mutual-kNN (reciprocal edges only)    (default false)
%                       false = symmetric union
%   'verbose'   logical                                       (default true)
%
% OUTPUT
%   A     [nObs x nObs] sparse symmetric adjacency, zero diagonal.
%   info  struct: .k, .weighting, .mutual, .sigma, .n_edges, .mean_degree
%
% See also: cycif_louvain, cycif_leiden.

p = inputParser();
addRequired(p, 'X', @(v) isnumeric(v) && ismatrix(v) && ~isempty(v));
addParameter(p, 'graph_k',   15);
addParameter(p, 'weighting', 'gaussian');
addParameter(p, 'mutual',    false);
addParameter(p, 'verbose',   true);
parse(p, X, varargin{:});
k         = round(p.Results.graph_k);
weighting = lower(string(p.Results.weighting));
mutual    = logical(p.Results.mutual);
verbose   = logical(p.Results.verbose);

X = double(X);
assert(all(isfinite(X(:))), 'cycif_knn_graph:nonFinite', 'X has NaN/Inf.');
n = size(X,1);
assert(n > k, 'cycif_knn_graph:tooFewObs', 'graph_k (%d) must be < nObs (%d).', k, n);

if verbose, fprintf('[cycif_knn_graph] building %d-NN graph on %d cells...\n', k, n); end
[nnIdx, nnDist] = i_knn(X, k+1);       % +1: first neighbor is self
nnIdx  = nnIdx(:, 2:end);              % drop self column
nnDist = nnDist(:, 2:end);

% edge weights
sigma = NaN;
switch weighting
    case "binary"
        w = ones(n, k);
    case "gaussian"
        sigma = median(nnDist(:));
        if sigma <= 0 || ~isfinite(sigma), sigma = 1; end
        w = exp(-(nnDist.^2) ./ (sigma^2));
    otherwise
        error('cycif_knn_graph:weighting', 'weighting must be gaussian | binary.');
end

% assemble sparse adjacency (symmetrize)
ii = repmat((1:n)', 1, k);
A  = sparse(ii(:), nnIdx(:), w(:), n, n);
if mutual
    A = min(A, A.');      % reciprocal edges only
else
    A = max(A, A.');      % symmetric union
end
A = A - diag(diag(A));    % no self loops

nE = nnz(A)/2;
info = struct('k', k, 'weighting', char(weighting), 'mutual', mutual, ...
    'sigma', sigma, 'n_edges', nE, 'mean_degree', 2*nE/n);
if verbose
    if mutual, mtag = ' (mutual)'; else, mtag = ''; end
    fprintf('[cycif_knn_graph] %d edges, mean degree %.1f%s.\n', nE, 2*nE/n, mtag);
end
end


% ======================================================================
function [idx, d] = i_knn(X, k)
% k nearest neighbors. knnsearch (Stats TB) if present, else chunked pdist2.
if exist('knnsearch', 'file')
    [idx, d] = knnsearch(X, X, 'K', k);
    return;
end
n = size(X,1);
idx = zeros(n, k); d = zeros(n, k);
chunk = max(1, floor(2e7 / max(1,n)));
for s = 1:chunk:n
    e = min(s+chunk-1, n);
    D = pdist2(X(s:e,:), X);
    [ds, is] = sort(D, 2, 'ascend');
    idx(s:e,:) = is(:, 1:k);
    d(s:e,:)   = ds(:, 1:k);
end
end