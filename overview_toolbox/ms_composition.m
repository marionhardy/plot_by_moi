function C = ms_composition(S, clustCol, byCol, varargin)
%CT_CLUSTER_COMPOSITION  Stacked barplot of cluster composition by a metadata
%   variable — which treatments/conditions populate each cluster.
%
%   Works for latent OR CyCIF clusters (any label column), grouped by any
%   cellT metadata (tx2_label, condition, rep, ...). Shows population overlap:
%   is a cluster treatment-specific, or shared across conditions?
%
% USAGE:
%   ms_composition(S, 'clust_latent', 'tx2_label');
%   ms_composition(S, 'clust_cycif',  'condition', 'Norm','fraction');
%   ms_composition(S, 'clust_latent', 'clust_cycif');  % cluster x cluster
%
% INPUT:
%   S         object from dl_build_live_cycif_tsc
%   clustCol  cluster label column (bars = clusters)
%   byCol     metadata column to break each bar down by (the stack colours)
%
% NAME-VALUE:
%   'Norm'      'fraction' (each bar sums to 1) | 'count'   ['fraction']
%   'DropNoise' drop cluster label 0/""/"noise"            [true]
%   'Flip'      swap: bars = byCol, stacks = clusters       [false]
%
% OUTPUT:
%   C  struct: .counts [nCluster x nLevel], .clusters, .levels, .frac

q.Norm='fraction'; q.DropNoise=true; q.Flip=false;
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'norm',      q.Norm=lower(varargin{k+1}); k=k+2;
        case 'dropnoise', q.DropNoise=logical(varargin{k+1}); k=k+2;
        case 'flip',      q.Flip=logical(varargin{k+1}); k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

getcol=@(c) S.cellT.(char(c));
cl  = string(getcol(clustCol));
by  = string(getcol(byCol));

% keep cells with both labels present
ok = cl~="" & by~="" & ~ismissing(cl) & ~ismissing(by);
if q.DropNoise, ok = ok & cl~="0" & lower(cl)~="noise"; end
cl=cl(ok); by=by(ok);

clu = unique(cl); lev = unique(by);
counts = zeros(numel(clu), numel(lev));
for i=1:numel(clu)
    for j=1:numel(lev)
        counts(i,j) = sum(cl==clu(i) & by==lev(j));
    end
end

if q.Flip
    counts=counts'; tmp=clu; clu=lev; lev=tmp;
    xlab=char(byCol); leglab=char(clustCol);
else
    xlab=char(clustCol); leglab=char(byCol);
end

frac = counts ./ max(sum(counts,2),1);
Y = counts; if strcmp(q.Norm,'fraction'), Y=frac; end

figure('Color','w','Name',sprintf('composition: %s by %s',char(clustCol),char(byCol)));
bar(Y,'stacked');
set(gca,'XTick',1:numel(clu),'XTickLabel',cellstr(clu),'XTickLabelRotation',45, ...
    'TickLabelInterpreter','none');
xlabel(xlab,'Interpreter','none');
ylabel(ternary(strcmp(q.Norm,'fraction'),'fraction of cluster','cell count'));
legend(cellstr(lev),'Location','bestoutside','Interpreter','none');
title(sprintf('%s composition by %s',char(clustCol),char(byCol)),'Interpreter','none');

C=struct('counts',counts,'frac',frac,'clusters',{cellstr(clu)},'levels',{cellstr(lev)});
end
function o=ternary(c,a,b); if c,o=a; else,o=b; end; end
