function S = ms_embedding(S, varargin)
%CT_TSC_EMBEDDING  Embed the AE latent (S.Lat) to 2-D and store in S.cellT.
%
%   TSC analog of cycif_plot_embedding, adapted: we embed the LEARNED LATENT
%   (S.Lat), NOT raw traces — the autoencoder already did representation
%   learning. Latent dims are on comparable scales, so standardization is
%   OFF by default (unlike CyCIF markers). Coordinates are written back into
%   S.cellT (umap1/umap2 or tsne1/tsne2) so they persist for reuse (colouring,
%   interactive trace pull-out). Row order is preserved (aligned to data/cycifT).
%
% USAGE:
%   S = ms_embedding(S, 'Method','umap-py');        % real UMAP via Python
%   S = ms_embedding(S, 'Method','umap');           % run_umap (MATLAB, if installed)
%   S = ms_embedding(S, 'Method','pca');            % just refresh PC1/PC2
%   S = ms_embedding(S, 'Method','umap-py','NNeighbors',30,'MinDist',0.1);
%
% NAME-VALUE:
%   'Method'      'umap-py' | 'umap' | 'tsne' | 'pca'   ['umap']
%                 'umap-py' = umap-learn via active pyenv (recommended if
%                 run_umap MATLAB isn't installed). Writes umap1/umap2.
%                 A PCA fallback on any failure is stored as pca1/pca2 and
%                 flagged in S.audit.embed_ok=false (honest labelling).
%   'Standardize' z-score latent dims first      [false]
%   'NNeighbors'  UMAP n_neighbors               [30]
%   'MinDist'     UMAP min_dist                  [0.1]
%   'Perplexity'  tSNE perplexity                [30]
%   'Seed'        RNG seed                       [0]
%   'Verbose'     [true]
%
% OUTPUT: S with S.cellT.(umap1/umap2 | tsne1/tsne2) added (PC1/PC2 refreshed
%   if Method='pca'), and S.audit.embed_method / .embed_ok recorded.

q.Method='umap'; q.Standardize=false; q.NNeighbors=30; q.MinDist=0.1;
q.Perplexity=30; q.Seed=0; q.Source='lat'; q.Prefix=''; q.Verbose=true;
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'method',      q.Method=lower(varargin{k+1}); k=k+2;
        case 'standardize', q.Standardize=logical(varargin{k+1}); k=k+2;
        case 'nneighbors',  q.NNeighbors=varargin{k+1}; k=k+2;
        case 'mindist',     q.MinDist=varargin{k+1}; k=k+2;
        case 'perplexity',  q.Perplexity=varargin{k+1}; k=k+2;
        case 'seed',        q.Seed=varargin{k+1}; k=k+2;
        case 'source',      q.Source=lower(varargin{k+1}); k=k+2;   % 'lat' | 'cycif'
        case 'prefix',      q.Prefix=varargin{k+1}; k=k+2;          % column-name prefix
        case 'verbose',     q.Verbose=logical(varargin{k+1}); k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

% --- pick the matrix to embed + default column prefix ---
switch q.Source
    case 'lat'
        X = S.Lat;
        if isempty(q.Prefix), q.Prefix=''; end          % umap1/umap2, PC1/PC2
        % CyCIF markers are z-scored by default (very different ranges);
        % latent is not (learned, comparable scales). Flip the default here:
    case 'cycif'
        X = table2array(S.cycifT);
        if isempty(q.Prefix), q.Prefix='cy'; end        % cyumap1/cyumap2, cyPC1/cyPC2
        if ~any(strcmpi(varargin(1:2:end),'standardize')), q.Standardize=true; end % markers: z-score by default
    otherwise, error('Source must be ''lat'' or ''cycif''.');
end
assert(size(X,1)==height(S.cellT), 'source rows != cellT rows — S is misaligned.');
if q.Standardize
    mu=mean(X); sd=std(X); sd(sd==0)=1; X=(X-mu)./sd;
end
rng(q.Seed);

P = q.Prefix;   % column-name prefix ('' for latent, 'cy' for cycif)
switch q.Method
    case 'umap'
        try
            emb = run_umap(X, 'n_neighbors',q.NNeighbors, 'min_dist',q.MinDist, 'verbose','none');
            S.cellT.([P 'umap1'])=emb(:,1); S.cellT.([P 'umap2'])=emb(:,2);
            S.audit.embed_method='umap'; S.audit.embed_ok=true;
        catch ME
            warning(['run_umap (MATLAB) failed (%s). Falling back to PCA, stored as ' ...
                '%spca1/%spca2 (NOT umap). For real UMAP use Method=''umap-py''.'], ME.message, P, P);
            [~,sc]=pca(X); S.cellT.([P 'pca1'])=sc(:,1); S.cellT.([P 'pca2'])=sc(:,2);
            S.audit.embed_method='pca(umap-fallback)'; S.audit.embed_ok=false;
        end
    case 'umap-py'
        try
            L = py.numpy.array(X);
            umod = py.importlib.import_module('umap');
            reducer = umod.UMAP(pyargs('n_neighbors',int32(q.NNeighbors), ...
                                       'min_dist',double(q.MinDist), ...
                                       'random_state',int32(q.Seed)));
            emb = double(reducer.fit_transform(L));
            S.cellT.([P 'umap1'])=emb(:,1); S.cellT.([P 'umap2'])=emb(:,2);
            S.audit.embed_method='umap-py'; S.audit.embed_ok=true;
        catch ME
            warning(['umap-py failed (%s). Is umap-learn in the active pyenv? ' ...
                'Falling back to PCA as %spca1/%spca2.'], ME.message, P, P);
            [~,sc]=pca(X); S.cellT.([P 'pca1'])=sc(:,1); S.cellT.([P 'pca2'])=sc(:,2);
            S.audit.embed_method='pca(umap-py-fallback)'; S.audit.embed_ok=false;
        end
    case 'tsne'
        perp=min(q.Perplexity, floor(size(X,1)/4));
        emb=tsne(X,'NumDimensions',2,'Perplexity',perp,'Standardize',false);
        S.cellT.([P 'tsne1'])=emb(:,1); S.cellT.([P 'tsne2'])=emb(:,2);
        S.audit.embed_method='tsne'; S.audit.embed_ok=true;
    case 'pca'
        [~,sc]=pca(X);
        if isempty(P)
            S.cellT.PC1=sc(:,1); S.cellT.PC2=sc(:,2);
        else
            S.cellT.([P 'PC1'])=sc(:,1); S.cellT.([P 'PC2'])=sc(:,2);
        end
        S.audit.embed_method='pca'; S.audit.embed_ok=true;
    otherwise, error('Method must be umap-py|umap|tsne|pca');
end
if q.Verbose
    fprintf('ms_embedding: %s on %s [%d x %d]%s -> S.cellT.%s*\n', ...
        q.Method, q.Source, size(X,1), size(X,2), ...
        ternary(q.Standardize,' (z-scored)',''), P);
end
end
function o=ternary(c,a,b); if c, o=a; else, o=b; end; end
