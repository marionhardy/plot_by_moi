function P = ms_cluster_profile(S, clustCol, varargin)
%CT_TSC_CLUSTER_PROFILE  Characterise latent clusters by BOTH their live-trace
%   dynamics (mean trace per cluster) AND their CyCIF endpoint markers.
%
%   The unique cross-modal readout of this dataset: for each dynamic cluster,
%   show (LEFT) the mean sensor trajectory that DEFINES it, and (RIGHT) the
%   endpoint protein state it ENDS in. Same cluster order on both sides so the
%   eye maps "this dynamic -> this protein state" directly. This is the
%   trace-space analog of cycif_cluster's marker heatmap, plus the marker
%   heatmap itself, on one figure.
%
% USAGE:
%   P = ms_cluster_profile(S, 'clust_latent');
%   P = ms_cluster_profile(S, clustID_vector);       % pass labels directly
%   P = ms_cluster_profile(S, 'clust_latent', 'Log2Markers',true, 'ZMarkers',true);
%
% INPUT:
%   S         struct from dl_build_live_cycif_tsc (needs .data, .cycifT, .cellT,
%             .meta.markerNames; channels inferred from size(data,3))
%   clustCol  name of a cluster column in S.cellT, OR a [nCells x 1] label vec
%
% NAME-VALUE:
%   'Channels'    cellstr of channel names for data's 3rd dim
%                 [{'RAMPKAR2','HYLIGHT'} if size(data,3)==2, else auto]
%   'Log2Markers' log2(max(x,1)) markers before averaging   [true]
%   'ZMarkers'    z-score each marker across clusters for the heatmap [true]
%   'Band'        trace spread: 'sem' | 'std' | 'none'      ['sem']
%   'DropNoise'   exclude cluster label 0 (DBSCAN noise)    [true]
%   'Verbose'     [true]
%
% OUTPUT P (struct):
%   .clusters      cluster ids plotted (in order)
%   .meanTrace     {nCh} each [nClust x nT] mean trace per cluster
%   .markerMean    [nClust x nMk] mean (log2) marker per cluster
%   .markerZ       [nClust x nMk] z-scored version shown in heatmap
%   .markerNames   1 x nMk
%   .nPerCluster   [nClust x 1]

q.Channels={}; q.Log2Markers=true; q.ZMarkers=true; q.Band='sem';
q.DropNoise=true; q.Verbose=true;
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'channels',    q.Channels=varargin{k+1}; k=k+2;
        case 'log2markers', q.Log2Markers=logical(varargin{k+1}); k=k+2;
        case 'zmarkers',    q.ZMarkers=logical(varargin{k+1}); k=k+2;
        case 'band',        q.Band=lower(varargin{k+1}); k=k+2;
        case 'dropnoise',   q.DropNoise=logical(varargin{k+1}); k=k+2;
        case 'verbose',     q.Verbose=logical(varargin{k+1}); k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

% --- resolve cluster labels ---
if isnumeric(clustCol) || iscategorical(clustCol)
    lab = clustCol(:);
else
    assert(ismember(char(clustCol),S.cellT.Properties.VariableNames), ...
        'clustCol "%s" not in cellT.', char(clustCol));
    lab = S.cellT.(char(clustCol));
end
if iscategorical(lab), lab=double(lab); end
lab = double(lab);

% --- channels for the trace panel ---
nCh = size(S.data,3);
if isempty(q.Channels)
    if isfield(S.meta,'channels') && numel(S.meta.channels)==nCh
        ch = cellstr(S.meta.channels);
    elseif nCh==2
        ch = {'RAMPKAR2','HYLIGHT'};
    else
        ch = arrayfun(@(i)sprintf('ch%d',i),1:nCh,'Un',0);
    end
else
    ch = cellstr(q.Channels);
end

% --- cluster set ---
cl = unique(lab(~isnan(lab)));
if q.DropNoise, cl(cl==0)=[]; end
nClust = numel(cl);
nT = size(S.data,2);

% --- LEFT: mean trace per cluster per channel ---
meanTrace = cell(1,nCh); bandTrace = cell(1,nCh);
for c=1:nCh
    meanTrace{c}=nan(nClust,nT); bandTrace{c}=nan(nClust,nT);
end
nPer=zeros(nClust,1);
for i=1:nClust
    rows = lab==cl(i); nPer(i)=sum(rows);
    for c=1:nCh
        M = S.data(rows,:,c);
        meanTrace{c}(i,:) = mean(M,1,'omitnan');
        switch q.Band
            case 'sem', bandTrace{c}(i,:)=std(M,0,1,'omitnan')./sqrt(max(sum(~isnan(M),1),1));
            case 'std', bandTrace{c}(i,:)=std(M,0,1,'omitnan');
            otherwise,  bandTrace{c}(i,:)=zeros(1,nT);
        end
    end
end

% --- RIGHT: mean marker per cluster ---
mk = S.cycifT.Properties.VariableNames;
Mk = table2array(S.cycifT);
if q.Log2Markers, Mk = log2(max(Mk,1)); end
markerMean=nan(nClust,numel(mk));
for i=1:nClust
    markerMean(i,:) = mean(Mk(lab==cl(i),:),1,'omitnan');
end
markerZ = markerMean;
if q.ZMarkers
    mu=mean(markerMean,1,'omitnan'); sd=std(markerMean,0,1,'omitnan'); sd(sd==0)=1;
    markerZ=(markerMean-mu)./sd;
end

% --- PLOT: traces (left column, one per channel) + heatmap (right) ---
figure('Color','w','Name','TSC cluster profile: dynamics + endpoint');
tl=tiledlayout(nCh, 2, 'TileSpacing','compact','Padding','compact');
cols=lines(nClust);
t=1:nT;
for c=1:nCh
    ax=nexttile(tl,(c-1)*2+1); hold(ax,'on');
    for i=1:nClust
        m=meanTrace{c}(i,:); b=bandTrace{c}(i,:);
        if ~strcmp(q.Band,'none')
            fill(ax,[t fliplr(t)],[m+b fliplr(m-b)],cols(i,:), ...
                'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
        end
        plot(ax,t,m,'Color',cols(i,:),'LineWidth',1.5, ...
            'DisplayName',sprintf('c%d (n=%d)',cl(i),nPer(i)));
    end
    ylabel(ax,ch{c},'Interpreter','none'); xlim(ax,[1 nT]);
    if c==1, title(ax,'Mean trace per cluster'); legend(ax,'Location','best'); end
    if c==nCh, xlabel(ax,'time (frames)'); end
end
% marker heatmap spanning the right column
axh=nexttile(tl,2,[nCh 1]);
imagesc(axh, markerZ);
colormap(axh, parula); cb=colorbar(axh);
cb.Label.String = ternary(q.ZMarkers,'z-scored mean (log2)','mean log2');
set(axh,'XTick',1:numel(mk),'XTickLabel',mk,'XTickLabelRotation',45, ...
    'TickLabelInterpreter','none', ...
    'YTick',1:nClust,'YTickLabel',compose('c%d',cl));
title(axh,'Endpoint CyCIF markers per cluster');
title(tl,'Latent clusters: live-trace dynamics (left) vs endpoint protein state (right)');

P=struct('clusters',cl,'meanTrace',{meanTrace},'markerMean',markerMean, ...
    'markerZ',markerZ,'markerNames',{mk},'nPerCluster',nPer);
if q.Verbose
    fprintf('ms_cluster_profile: %d clusters | %d channels | %d markers.\n', ...
        nClust, nCh, numel(mk));
    fprintf('  cells per cluster:'); fprintf(' c%d:%d',[cl'; nPer']); fprintf('\n');
end
end
function o=ternary(c,a,b); if c, o=a; else, o=b; end; end
