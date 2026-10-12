function [labels, info] = cycif_louvain(X, varargin)
% cycif_louvain  Graph-based clustering via k-NN graph + greedy Louvain
%                modularity optimization. Self-contained (no toolboxes).
%
%   Clusters rows of X (cells x features) by:
%     1. building a mutual/symmetric k-NN similarity graph,
%     2. optimizing modularity with the Louvain method (local moving +
%        aggregation, iterated over levels).
%
%   Designed for marker/PCA feature spaces where k-means' convex-blob and
%   fixed-k assumptions fail and DBSCAN's single global density does not fit a
%   variable-density landscape. Louvain needs no cluster count -- granularity
%   is set by graph_k and the resolution gamma.
%
% Usage examples:
%   labels = cycif_louvain(X);                              % defaults
%   labels = cycif_louvain(X, 'graph_k',15, 'gamma',1.0);
%   [labels, info] = cycif_louvain(X, 'gamma',1.4, 'seed',0);  % finer clusters
%
% INPUT
%   X  [nObs x nFeat] double, finite. Feed cycif_prep_features output (arcsinh/
%      PCA/PC-dropped), NOT the 2-D embedding.
%
% NAME-VALUE
%   'graph_k'  scalar   neighbors per node in the k-NN graph      (default 15)
%   'gamma'    scalar   resolution: >1 smaller/more clusters,     (default 1.0)
%                       <1 larger/fewer; 1 = classic modularity
%   'weighting' string  'binary' | 'gaussian' edge weights        (default 'gaussian')
%                       gaussian: exp(-d^2 / sigma^2), sigma = median kNN dist
%   'mutual'   logical  keep only mutual k-NN edges (sparser,     (default false)
%                       cleaner communities); false = symmetrized union
%   'max_levels' scalar Louvain aggregation levels cap            (default 10)
%   'seed'     scalar   RNG seed (node visit order is randomized) (default 0)
%   'verbose'  logical  print progress                            (default true)
%
% OUTPUT
%   labels  [nObs x 1] double community IDs, 1..K (contiguous, size-ranked so
%           1 is the largest community). No noise label -- Louvain assigns all.
%   info    struct: .k, .gamma, .modularity (final Q), .n_communities,
%                   .n_levels, .sizes (per-community counts)
%
% NOTES
% - Complexity is ~O(nObs * graph_k) memory for the sparse graph; fine for
%   ~50-100k cells. For >few hundred k, subsample or raise sparsity via mutual.
% - Modularity Q typically 0.3-0.7 for well-structured data; very low Q (<0.1)
%   means little community structure at this gamma / in this feature space.
%
% See also: cycif_prep_features, cycif_cluster.

p = inputParser();
addRequired(p, 'X', @(v) isnumeric(v) && ismatrix(v) && ~isempty(v));
addParameter(p, 'graph_k',    15);
addParameter(p, 'gamma',      1.0);
addParameter(p, 'weighting',  'gaussian');
addParameter(p, 'mutual',     false);
addParameter(p, 'max_levels', 10);
addParameter(p, 'seed',       0);
addParameter(p, 'verbose',    true);
parse(p, X, varargin{:});
k        = round(p.Results.graph_k);
gamma    = p.Results.gamma;
weighting= lower(string(p.Results.weighting));
mutual   = logical(p.Results.mutual);
maxLvl   = p.Results.max_levels;
seed     = p.Results.seed;
verbose  = logical(p.Results.verbose);

X = double(X);
assert(all(isfinite(X(:))), 'cycif_louvain:nonFinite', 'X has NaN/Inf.');
n = size(X,1);
assert(n > k, 'cycif_louvain:tooFewObs', 'graph_k (%d) must be < nObs (%d).', k, n);
rng(seed);

% ---- 1. k-NN graph (shared builder: same graph as cycif_leiden) --------
[A, gInfo] = cycif_knn_graph(X, 'graph_k', k, 'weighting', char(weighting), ...
                             'mutual', mutual, 'verbose', verbose);

% ---- 2. Louvain modularity optimization --------------------------------
if verbose, fprintf('[cycif_louvain] optimizing modularity (gamma=%.3g)...\n', gamma); end
[comm, Q, nLvl] = i_louvain(A, gamma, maxLvl, verbose);

% ---- relabel: contiguous, size-ranked (1 = largest) --------------------
[labels, sizes] = i_rank_labels(comm);

info = struct('k', k, 'gamma', gamma, 'modularity', Q, ...
    'n_communities', numel(sizes), 'n_levels', nLvl, 'sizes', sizes, ...
    'graph', gInfo);
if verbose
    fprintf('[cycif_louvain] %d communities, Q=%.4f, %d level(s). Sizes: %s\n', ...
        numel(sizes), Q, nLvl, mat2str(sizes(:)'));
end
end



function [comm, Q, level] = i_louvain(A, gamma, maxLvl, verbose)
% Greedy Louvain: local moving + aggregation, iterated across levels.
n0 = size(A,1);
comm = (1:n0)';          % final membership over ORIGINAL nodes
Acur = A;                % adjacency at current aggregation level
mapToOrig = (1:n0)';     % current-node -> original-node group id (via comm)
level = 0;

while level < maxLvl
    level = level + 1;
    [c, improved] = i_local_moving(Acur, gamma);
    % propagate this level's assignment down to original nodes
    comm = c(comm);      % comm currently maps orig->cur; c maps cur->newcur
    if ~improved || numel(unique(c)) == numel(c)
        break;           % no merges this level -> converged
    end
    % aggregate: build next-level graph where nodes are communities
    Acur = i_aggregate(Acur, c);
    if size(Acur,1) == 1, break; end
end

Q = i_modularity(A, comm, gamma);
end


function [c, improved] = i_local_moving(A, gamma)
% One pass (repeated to convergence) of moving nodes to the neighbor
% community giving the largest modularity gain. Resolution gamma scales the
% null-model term.
n = size(A,1);
deg = full(sum(A,2));            % weighted degree
m2  = sum(deg);                  % 2m (total edge weight *2)
if m2 == 0, c = (1:n)'; improved = false; return; end

% --- Precompute neighbor lists ONCE (avoids O(n) sparse row-slicing per
% node). A is symmetric, so columns == rows; column access is native to CSC.
% Store each node's neighbor indices + weights in cell arrays.
[nbAll, ~, wAll] = find(A);      % all edges (col-major: nbAll grouped by column)
colPtr = cumsum([1; full(sum(A~=0, 1)')]);   % start index per node in the find() order
nbList = cell(n,1); wList = cell(n,1);
for j = 1:n
    r = colPtr(j):(colPtr(j+1)-1);
    nbList{j} = nbAll(r);
    wList{j}  = wAll(r);
end

c = (1:n)';                      % start: each node its own community
sigTot = deg;                    % community total degree

improved = false;
anyMove = true;
order = randperm(n);            % randomized visit order (seeded upstream)
guard = 0;
while anyMove && guard < 100
    guard = guard + 1;
    anyMove = false;
    for oi = 1:n
        u = order(oi);
        cu = c(u);
        du = deg(u);
        sigTot(cu) = sigTot(cu) - du;      % remove u from its community

        nb = nbList{u};                    % O(degree) neighbor lookup
        if isempty(nb)
            sigTot(cu) = sigTot(cu) + du;
            continue;
        end
        wUto = wList{u};
        cn   = c(nb);
        % aggregate edge weight to each candidate community
        [uc, ~, g] = unique(cn);
        kin = accumarray(g, wUto(:), [numel(uc) 1]);   % sum weight u->comm

        % Ensure the original community cu is a candidate even if u has no
        % edges back to it (kin=0), so "stay" is evaluated on equal footing.
        if ~any(uc == cu)
            uc  = [uc; cu];              %#ok<AGROW>
            kin = [kin; 0];              %#ok<AGROW>
        end

        % modularity gain: dQ ~ kin - gamma * sigTot(comm) * du / m2
        gain = kin - gamma .* sigTot(uc) .* (du / m2);
        [bestGain, bi] = max(gain);
        bestC = uc(bi);
        if bestGain <= 0, bestC = cu; end   % all moves worsen Q -> stay

        c(u) = bestC;
        sigTot(bestC) = sigTot(bestC) + du;
        if bestC ~= cu
            anyMove = true; improved = true;
        end
    end
end
% relabel communities contiguous 1..K
[~, ~, c] = unique(c);
end


function Ag = i_aggregate(A, c)
% Collapse communities in c into super-nodes; edge weights summed.
K = max(c);
S = sparse(1:numel(c), c, 1, numel(c), K);   % node -> community indicator
Ag = S.' * A * S;                            % community-community weights
Ag = (Ag + Ag.')/2;                          % keep symmetric
end


function Q = i_modularity(A, c, gamma)
% Generalized modularity with resolution gamma.
deg = full(sum(A,2));
m2  = sum(deg);
if m2 == 0, Q = 0; return; end
K = max(c);
S = sparse(1:numel(c), c, 1, numel(c), K);
eIn  = full(sum((S.' * A * S) .* speye(K), 'all'));   % within-community weight
dTot = full(S.' * deg);                                % community total degree
Q = eIn/m2 - gamma * sum((dTot./m2).^2);
end


function [labels, sizes] = i_rank_labels(c)
% Relabel so community 1 is largest, 2 next, etc.
u = unique(c);
cnt = arrayfun(@(x) sum(c==x), u);
[sizes, ord] = sort(cnt, 'descend');
remap = zeros(max(u),1);
remap(u(ord)) = 1:numel(u);
labels = remap(c);
end