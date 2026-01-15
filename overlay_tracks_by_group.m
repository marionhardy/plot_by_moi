function overlay_tracks_from_ids(dataloc, GroupStruct, varargin)
% GroupStruct has numeric vectors: .Oscillatory, .NonOscillatory (global cell IDs)

ip = inputParser; ip.CaseSensitive=false;
addParameter(ip,'channel','ERKTR');
addParameter(ip,'group','osc',@(s) any(strcmpi(s,{'osc','non'})));
addParameter(ip,'ntracks',20,@(x) x>0 && isscalar(x));
addParameter(ip,'smooth',[],@(x) isempty(x)||isscalar(x));
addParameter(ip,'tile',true,@islogical);
addParameter(ip,'tmax_hours',[],@(x) isempty(x)||isscalar(x)); % truncate by hours from start
parse(ip,varargin{:}); p=ip.Results; chan=char(p.channel);

IDs = strcmpi(p.group,'osc') * GroupStruct.Oscillatory + ...
      strcmpi(p.group,'non') * GroupStruct.NonOscillatory;
if isempty(IDs), error('Selected group has no IDs.'); end
IDs = IDs(:);

% frames per hour
looptime = 6; if isfield(dataloc,'movieinfo') && isfield(dataloc.movieinfo,'tsamp'), looptime = dataloc.movieinfo.tsamp; end
tpk = 60/looptime;

goodXY = find(cellfun(@(x) ~isempty(x) && isfield(x,'cellindex') && ~isempty(x.cellindex), dataloc.d));
if p.tile, figure('Color','w'); tl=tiledlayout('flow','TileSpacing','compact','Padding','compact'); end

for xy = goodXY(:)'
    D = dataloc.d{xy};
    if ~isfield(D,'data') || ~isfield(D.data,chan), continue; end
    cellIDs_xy = D.cellindex(:);
    rows = find(ismember(cellIDs_xy, IDs));
    if isempty(rows), continue; end
    if numel(rows) > p.ntracks, rows = rows(randperm(numel(rows), p.ntracks)); end

    X = D.data.(chan);
    lasttp = size(X,2);
    if ~isempty(p.tmax_hours), lasttp = min(lasttp, floor(p.tmax_hours*tpk)); end
    dat = X(rows, 1:lasttp);
    if ~isempty(p.smooth), dat = movmean(dat, p.smooth, 2, 'omitnan'); end

    if p.tile, ax = nexttile(tl); else, figure('Color','w'); ax=gca; end
    hold(ax,'on'); plot(ax, dat','LineWidth',0.8);
    plot(ax, mean(dat,1,'omitnan'),'k-','LineWidth',2);
    title(ax, sprintf('XY %d — %s (%s)', xy, chan, p.group));
    xlabel(ax,'frames'); ylabel(ax, strrep(chan,'_',' '));
end
end
