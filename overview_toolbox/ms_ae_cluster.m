function [clustID, info] = ms_ae_cluster(Lat, varargin)
%CT_AE_CLUSTER  Cluster autoencoder latent traces in the FULL latent space.
%
%   Clusters the AE latent `Lat` [nCells x latentDim] directly (not a 2-D
%   embedding), because the latent is the distance-faithful space the AE
%   optimized and there are no interpretable markers to backstop an
%   embedding artefact (unlike CyCIF). A 2-D UMAP/tSNE, if used, is for
%   VISUALISATION ONLY — colour it by these labels, never cluster it.
%
%   Runs k-means and/or DBSCAN (auto-epsilon), reports which, and returns
%   labels aligned row-for-row with Lat so they drop straight into the
%   pooled provenance / batch-diagnostic slot (clustID).
%
% USAGE:
%   clustID = ms_ae_cluster(Lat);                       % both methods, kmeans primary
%   clustID = ms_ae_cluster(Lat,'Method','dbscan');     % density-based (blob+outliers)
%   [clustID,info] = ms_ae_cluster(Lat,'Method','both','K',6,'Standardize',true);
%
% INPUT:
%   Lat            [nCells x latentDim] AE latent (e.g. from tsc.encode)
%
% NAME-VALUE:
%   'Method'       'kmeans' | 'dbscan' | 'both'     [both]
%   'K'            k for k-means                     [8]
%   'Epsilon'      DBSCAN radius ([]=auto)           []
%   'MinPts'       DBSCAN min neighbours             [max(10, latentDim+1)]
%   'Standardize'  z-score each latent dim first     [true]
%   'Primary'      which labels to RETURN as clustID: 'kmeans'|'dbscan' [kmeans]
%   'Verbose'      print progress                    [true]
%
% OUTPUT:
%   clustID  [nCells x 1] contiguous integer labels (1..K); DBSCAN noise = 0
%   info     struct: .km, .db (raw labels), .epsilon, .minpts, .K,
%            .silhouette_km, .method, .standardized
%
% NOTE ON HIGH-D: in latentDim>~10 Euclidean distances concentrate, which
% weakens DBSCAN. MinPts defaults scale with dimension; if DBSCAN returns
% ~1 cluster + heavy noise, prefer k-means or reduce with PCA first.

% ---- parse -------------------------------------------------------------
q.Method='both'; q.K=8; q.Epsilon=[]; q.MinPts=[]; q.Standardize=true;
q.Primary='kmeans'; q.Verbose=true;
q.GMMMaxComp=6; q.GMMVarKeep=0.95;   % ADDED: GMM sweep 1..GMMMaxComp; fit in PCs covering GMMVarKeep variance
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'method',      q.Method=lower(varargin{k+1}); k=k+2;
        case 'k',           q.K=varargin{k+1};             k=k+2;
        case 'epsilon',     q.Epsilon=varargin{k+1};       k=k+2;
        case 'minpts',      q.MinPts=varargin{k+1};        k=k+2;
        case 'standardize', q.Standardize=logical(varargin{k+1}); k=k+2;
        case 'primary',     q.Primary=lower(varargin{k+1}); k=k+2;
        case 'gmmmaxcomp',  q.GMMMaxComp=varargin{k+1};    k=k+2;   % ADDED
        case 'gmmvarkeep',  q.GMMVarKeep=varargin{k+1};    k=k+2;   % ADDED
        case 'verbose',     q.Verbose=logical(varargin{k+1}); k=k+2;
        otherwise, error('ms_ae_cluster: unknown param "%s".', varargin{k});
    end
end

[nCells, dLat] = size(Lat);
assert(nCells>1 && dLat>=1, 'Lat must be [nCells x latentDim].');
if isempty(q.MinPts); q.MinPts = max(10, dLat+1); end   % scale with dimension

% ---- standardize (equal weight per latent dim) -------------------------
X = Lat;
if q.Standardize
    mu = mean(X,1); sd = std(X,0,1); sd(sd==0) = 1;
    X  = (X - mu) ./ sd;
end

info = struct('km',[], 'db',[], 'epsilon',[], 'minpts',q.MinPts, ...
              'K',q.K, 'silhouette_km',NaN, 'method',q.Method, ...
              'standardized',q.Standardize, ...
              'gmm',[], 'gmm_bic',[], 'gmm_nopt',NaN, 'gmm_npc',NaN);   % ADDED: GMM outputs

runKM  = ismember(q.Method,{'kmeans','both'});
runDB  = ismember(q.Method,{'dbscan','both'});
runGMM = ismember(q.Method,{'gmm','both'});   % ADDED

% ---- k-means -----------------------------------------------------------
if runKM
    if q.Verbose; fprintf('ms_ae_cluster: k-means (k=%d) on %d-D latent...\n', q.K, dLat); end
    km = kmeans(X, q.K, 'Replicates', 5, 'MaxIter', 500);
    info.km = km;
    if q.Verbose
        s = silhouette(X, km); info.silhouette_km = mean(s);
        fprintf('  k-means sizes:'); for c=1:q.K; fprintf(' %d:%d',c,sum(km==c)); end
        fprintf(' | mean silhouette = %.3f\n', info.silhouette_km);
    end
end

% ---- DBSCAN (auto-epsilon) ---------------------------------------------
if runDB
    eps_use = q.Epsilon;
    if isempty(eps_use)
        % k-distance heuristic (same idea as cycif_cluster) in full latent D
        knn_d   = pdist2(X, X, 'euclidean', 'Smallest', q.MinPts+1);
        eps_use = median(knn_d(end,:)) * 1.5;
        if q.Verbose
            fprintf('ms_ae_cluster: DBSCAN auto epsilon = %.3f (1.5x median %d-NN dist)\n', ...
                eps_use, q.MinPts);
        end
    end
    db = dbscan(X, eps_use, q.MinPts);       % -1 = noise
    info.db = db; info.epsilon = eps_use;
    if q.Verbose
        nCl = numel(unique(db(db>0))); nNo = sum(db==-1);
        fprintf('  DBSCAN: %d clusters, %d noise (%.1f%%)\n', nCl, nNo, 100*nNo/nCells);
        if nCl <= 1
            fprintf(2,'  !! DBSCAN found <=1 cluster — high-D distance concentration.\n');
            fprintf(2,'     Prefer k-means, or PCA-reduce Lat before DBSCAN.\n');
        end
    end
end

% ---- GMM + BIC (continuum vs #states, in low-D PCA space) --------------
% ADDED: fit GMMs with 1..GMMMaxComp components on the PCs that carry
% GMMVarKeep of the variance. BIC minimum picks the model: 1 component = NO
% real structure (continuum) -> do not cluster; report the axis instead.
if runGMM
    [~, scPC, ~, ~, ex] = pca(X);
    npc = find(cumsum(ex)/sum(ex) >= q.GMMVarKeep, 1, 'first');
    npc = max(npc, 2);                       % keep >=2 so 2-D structure is representable
    Y   = scPC(:, 1:npc);
    info.gmm_npc = npc;
    if q.Verbose
        fprintf('ms_ae_cluster: GMM+BIC on %d PCs (%.0f%% variance), 1..%d comps...\n', ...
            npc, 100*sum(ex(1:npc))/sum(ex), q.GMMMaxComp);
    end
    bic = nan(1,q.GMMMaxComp); gm = cell(1,q.GMMMaxComp);
    for g = 1:q.GMMMaxComp
        try
            gm{g} = fitgmdist(Y, g, 'Replicates',5, 'RegularizationValue',1e-6, ...
                              'Options',statset('MaxIter',500));
            bic(g) = gm{g}.BIC;
        catch
            bic(g) = Inf;                    % ill-conditioned fit -> skip
        end
    end
    [~, nopt] = min(bic);
    info.gmm_bic = bic; info.gmm_nopt = nopt;
    if ~isinf(bic(nopt)); info.gmm = cluster(gm{nopt}, Y); end  % hard assignment
    if q.Verbose
        fprintf('  BIC by #components:'); fprintf(' %d:%.0f', [1:q.GMMMaxComp; bic]); fprintf('\n');
        fprintf('  BIC-optimal = %d component(s)', nopt);
        if nopt==1
            fprintf(2, '  -> CONTINUUM: one Gaussian best. Do NOT cluster; report the PC axis.\n');
        else
            fprintf(' -> %d states supported.\n', nopt);
        end
    end
end

% ---- choose primary labels to return -----------------------------------
switch q.Primary
    case 'kmeans'
        assert(runKM,'Primary=kmeans but Method excluded k-means.');
        raw = info.km;
    case 'dbscan'
        assert(runDB,'Primary=dbscan but Method excluded DBSCAN.');
        raw = info.db; raw(raw==-1) = 0;      % noise -> 0
    case 'gmm'                                                  % ADDED
        assert(runGMM && ~isempty(info.gmm),'Primary=gmm but GMM did not run/fit.'); % ADDED
        raw = info.gmm;                                         % ADDED: BIC-optimal GMM assignment
    otherwise, error('Primary must be ''kmeans'', ''dbscan'', or ''gmm''.'); % CHANGED: gmm now valid
end

% contiguous integer labels (noise stays 0 if present)
if any(raw==0)
    pos = raw>0; [~,~,c] = unique(raw(pos));
    clustID = zeros(nCells,1); clustID(pos) = c;   % 0 = noise, 1..K = clusters
else
    [~,~,clustID] = unique(raw);
end
end
