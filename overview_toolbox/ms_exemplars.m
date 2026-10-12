function ex = ms_exemplars_window(S, clustCol, win, varargin)
%MS_EXEMPLARS_WINDOW  Per-cluster exemplar traces: FULL trace next to a ZOOMED
%   window of timepoints you specify, for both sensors. See overall dynamics
%   and fine detail in a region of interest side by side.
%
%   Layout per cluster row: [ch1 full | ch1 window | ch2 full | ch2 window].
%   The window is shaded on the full-trace panels for reference.
%
% USAGE:
%   ms_exemplars_window(S, 'clust_latent', [80 120]);
%   ms_exemplars_window(S, 'clust_latent', [80 120], 'N',10);
%
% INPUT:
%   S         object from dl_build_live_cycif_tsc (.data, .cellT, .Lat)
%   clustCol  cellT column name (or [nCells x 1] label vector)
%   win       [t0 t1] timepoint window to zoom (frame indices)
%
% NAME-VALUE:
%   'N'         traces per cluster                       [10]
%   'Channels'  cellstr for data 3rd dim   [auto/{RAMPKAR2,HYLIGHT}]
%   'Pick'      'random' | 'representative'              ['representative']
%   'DropNoise' exclude label 0/""/"noise"               [true]
%   'Seed'      RNG for 'random'                         [0]
%
% OUTPUT:
%   ex  struct: .clusters, .rows{k} (cellT.row of picked cells), .win

q.N=10; q.Channels={}; q.Pick='representative'; q.DropNoise=true; q.Seed=0;
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'n',         q.N=varargin{k+1}; k=k+2;
        case 'channels',  q.Channels=varargin{k+1}; k=k+2;
        case 'pick',      q.Pick=lower(varargin{k+1}); k=k+2;
        case 'dropnoise', q.DropNoise=logical(varargin{k+1}); k=k+2;
        case 'seed',      q.Seed=varargin{k+1}; k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

nT = size(S.data,2);
assert(numel(win)==2 && win(1)>=1 && win(2)<=nT && win(1)<win(2), ...
    'win must be [t0 t1] within [1 %d].', nT);
wt = win(1):win(2);

if isnumeric(clustCol)||iscategorical(clustCol), lab=clustCol;
else, lab=S.cellT.(char(clustCol)); end
labs = string(lab);
cl = unique(labs);
if q.DropNoise, cl(cl=="0"|cl==""|lower(cl)=="noise")=[]; end
nClust=numel(cl);

nCh=size(S.data,3);
if isempty(q.Channels)
    if isfield(S.meta,'channels')&&numel(S.meta.channels)==nCh, ch=cellstr(S.meta.channels);
    elseif nCh==2, ch={'RAMPKAR2','HYLIGHT'};
    else, ch=arrayfun(@(i)sprintf('ch%d',i),1:nCh,'Un',0); end
else, ch=cellstr(q.Channels); end
t=1:nT;

rng(q.Seed);
cols = lines(nClust);
ex.clusters=cl; ex.rows=cell(1,nClust); ex.win=win;

figure('Color','w','Name',['exemplars + window: ' char(string(clustCol))]);
tl=tiledlayout(nClust, 2*nCh, 'TileSpacing','compact','Padding','compact');

for i=1:nClust
    idx = find(labs==cl(i));
    if strcmp(q.Pick,'representative') && ~isempty(S.Lat)
        c0 = mean(S.Lat(idx,:),1);
        dd = sum((S.Lat(idx,:)-c0).^2,2);
        [~,ord]=sort(dd,'ascend');
        pick = idx(ord(1:min(q.N,numel(idx))));
    else
        pick = idx(randperm(numel(idx), min(q.N,numel(idx))));
    end
    ex.rows{i} = S.cellT.row(pick);

    for c=1:nCh
        M  = S.data(pick,:,c);
        mu = mean(M,1,'omitnan');
        ylo = min(M(:),[],'omitnan'); yhi = max(M(:),[],'omitnan');
        if ~isfinite(ylo)||~isfinite(yhi)||ylo==yhi, ylo=0; yhi=1; end

        axF = nexttile(tl, (i-1)*2*nCh + (c-1)*2 + 1); hold(axF,'on');
        patch(axF, [win(1) win(2) win(2) win(1)], [ylo ylo yhi yhi], ...
            [0.85 0.85 0.85], 'EdgeColor','none', 'HandleVisibility','off');
        for j=1:size(M,1)
            h=plot(axF, t, M(j,:), 'Color',cols(i,:), 'LineWidth',0.5);
            try, h.Color(4)=0.35; catch, end
        end
        plot(axF, t, mu, 'Color',cols(i,:)*0.5, 'LineWidth',1.8);
        xlim(axF,[1 nT]);
        if c==1, ylabel(axF,sprintf('c%s (n=%d)',cl(i),numel(idx)),'Interpreter','none'); end
        if i==1, title(axF,[ch{c} ' (full)'],'Interpreter','none'); end
        if i==nClust, xlabel(axF,'time (frames)'); end

        axW = nexttile(tl, (i-1)*2*nCh + (c-1)*2 + 2); hold(axW,'on');
        for j=1:size(M,1)
            h=plot(axW, wt, M(j,wt), 'Color',cols(i,:), 'LineWidth',0.6);
            try, h.Color(4)=0.4; catch, end
        end
        plot(axW, wt, mu(wt), 'Color',cols(i,:)*0.5, 'LineWidth',1.8);
        xlim(axW,[win(1) win(2)]);
        if i==1, title(axW,sprintf('%s (tp %d-%d)',ch{c},win(1),win(2)),'Interpreter','none'); end
        if i==nClust, xlabel(axW,'time (frames)'); end
    end
end
title(tl, sprintf('%d exemplars/cluster - full + window tp %d-%d (%s pick)', ...
    q.N, win(1), win(2), q.Pick), 'Interpreter','none');
end