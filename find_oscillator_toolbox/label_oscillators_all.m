function osc_labels = label_oscillators_all(dataloc, chanName, maxCells, xyInclude, minTrackPts, varargin)
% label_oscillators_all  Manual triage across selected XYs for one channel,
% condition-aware sampling (does NOT bottleneck to the smallest condition),
% with optional min group size threshold.
%
% osc_labels = label_oscillators_all(dataloc, chanName, maxCells, xyInclude, minTrackPts)
% osc_labels = label_oscillators_all(..., 'minPerGroup', 20, 'seed', 1)
%
% Inputs:
%   dataloc      : main dataloc struct with .d cell array and .platemapd.idx
%   chanName     : channel field in S.data (e.g. 'HYLIGHT')
%   maxCells     : max number of cells to label (e.g. 1000). Use Inf for all.
%   xyInclude    : numeric XY list OR logical mask of allowed XYs ([] = all)
%   minTrackPts  : minimum # finite points required for a cell trace
%
% Name-value options:
%   'minPerGroup' : drop conditions with < this many eligible cells (default 20)
%   'seed'        : rng seed for reproducible shuffling (default [])
%   'expar'       : passed to ct_maptx (default {'Cell','pTx','Tx1','Tx2','Tx3','Tx4'})
%   'reduce'      : passed to ct_maptx (default false)
%
% Behavior:
%   - Builds condition groups via ct_maptx(dataloc.platemapd.idx)
%   - Restricts to XYs in xyInclude
%   - Finds eligible cells (>= minTrackPts finite points)
%   - Drops tiny groups (< minPerGroup)
%   - Samples across conditions WITHOUT being limited by smallest group:
%       * gives everyone a baseline,
%       * takes all from small groups,
%       * fills remaining quota from larger groups (no replacement)
%   - Interactive labeling keys:
%       o = oscillating
%       n = not oscillating
%       u = unsure
%       b = back
%       q = quit
%
% Output struct osc_labels:
%   osc_labels.cells(k) has fields:
%       xy, row, cellindex, full, PtxTx, Ptx, Tx
%   osc_labels.labels(k): 1 osc, 0 non, -1 unsure, NaN unlabeled
%   osc_labels.chanName, osc_labels.dt_minutes
%   osc_labels.minPerGroup, osc_labels.droppedGroups

% ---------------- defaults ----------------
p = struct();
p.minpergroup = 20;
p.seed        = [];
p.expar       = {'Cell','pTx','Tx1','Tx2','Tx3','Tx4'};
p.reduce      = false;

% ---------------- parse varargin ----------------
if rem(numel(varargin),2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    p.(lower(varargin{k})) = varargin{k+1};
end

if nargin < 3 || isempty(maxCells), maxCells = inf; end
if nargin < 4, xyInclude = []; end
if nargin < 5 || isempty(minTrackPts), minTrackPts = 1; end

if ~isempty(p.seed)
    rng(p.seed);
end

% ---------------- dt (hours) ----------------
dt_minutes = [];
if isfield(dataloc,'movieinfo') && isfield(dataloc.movieinfo,'tsamp') && ~isempty(dataloc.movieinfo.tsamp)
    ts = unique(dataloc.movieinfo.tsamp);
    dt_minutes = ts(1);
end
if isempty(dt_minutes) || ~isfinite(dt_minutes) || dt_minutes <= 0
    warning('Could not read dataloc.movieinfo.tsamp; defaulting dt=10 min.');
    dt_minutes = 10;
end
dt_hours = dt_minutes / 60;

% ---------------- build XY mask ----------------
nXY = numel(dataloc.d);
if isempty(xyInclude)
    xyMask = true(1, nXY);
else
    if islogical(xyInclude)
        xyMask = xyInclude(:).';
        if numel(xyMask) ~= nXY
            error('xyInclude logical mask must have length numel(dataloc.d).');
        end
    else
        % normalize xyInclude to numeric vector
        if iscell(xyInclude)
            xyInclude = cell2mat(xyInclude);
        end

        xyInclude = xyInclude(:).';        % row vector
        xyInclude = xyInclude(isfinite(xyInclude));
        xyInclude = unique(xyInclude);

        % bounds check
        xyInclude = xyInclude(xyInclude >= 1 & xyInclude <= nXY);

        xyMask = false(1, nXY);
        xyMask(xyInclude) = true;
    end
end

% also require dataloc.d{xy} exists
xyMask = xyMask & arrayfun(@(xy) ~isempty(dataloc.d{xy}), 1:nXY);

% ---------------- map conditions (ct_maptx) ----------------
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'idx')
    error('dataloc.platemapd.idx missing; cannot group by condition.');
end

mappedTxs = ct_maptx(dataloc.platemapd.idx, 'expar', p.expar, 'reduce', p.reduce);
txNames = fieldnames(mappedTxs);
if isempty(txNames)
    warning('ct_maptx returned no conditions.');
    osc_labels = struct('cells',struct([]),'labels',[],'chanName',chanName,'dt_minutes',dt_minutes);
    return;
end

% ---------------- collect eligible cells per condition ----------------
groups = struct('name',{}, 'cells',{}, 'nEligible',{});
Tmax_global = 0;

for iMap = 1:numel(txNames)
    nm = txNames{iMap};

    xyThis = mappedTxs.(nm).xy(:);
    xyThis = xyThis(isfinite(xyThis));
    xyThis = xyThis(xyThis>=1 & xyThis<=nXY);
    xyThis = xyThis(xyMask(xyThis));  % restrict to xyInclude

    if isempty(xyThis)
        continue;
    end

    % Build readable strings once per condition
    tpPerHr = 60/dt_minutes;
    [~, preStr, txStr, fullStr] = ct__buildTxStrings(mappedTxs.(nm).tx, tpPerHr);

    cellList = struct('xy',{}, 'row',{}, 'cellindex',{}, 'full',{}, 'PtxTx',{}, 'Ptx',{}, 'Tx',{});

    for ii = 1:numel(xyThis)
        xy = xyThis(ii);
        S = dataloc.d{xy};

        if isempty(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        if ~isfield(S,'cellindex') || isempty(S.cellindex)
            continue;
        end

        X = S.data.(chanName); % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end

        [nCellsXY, Txy] = size(X);
        Tmax_global = max(Tmax_global, Txy);

        for r = 1:nCellsXY
            tr = X(r,:);
        
            % ---- HARD REQUIREMENT: full-length trace ----
            if numel(tr) < Tmax_global
                continue;   % skip incomplete movies entirely
            end
        
            % ---- still enforce finite-point requirement ----
            effLen = sum(isfinite(tr));
            if effLen < minTrackPts
                continue;
            end
        
            c.xy = xy;
            c.row = r;
            c.cellindex = S.cellindex(r);
        
            c.full  = string(nm);
            c.PtxTx = string(fullStr);
            c.Ptx   = string(preStr);
            c.Tx    = string(txStr);
        
            cellList(end+1) = c; %#ok<AGROW>
        end
    end

    g.name = string(nm);
    g.cells = cellList;
    g.nEligible = numel(cellList);

    groups(end+1) = g; %#ok<AGROW>
end

if isempty(groups)
    warning('No eligible cells found after filtering XYs and minTrackPts.');
    osc_labels = struct('cells',struct([]),'labels',[],'chanName',chanName,'dt_minutes',dt_minutes);
    return;
end

% ---------------- drop tiny groups ----------------
nEligibleAll = [groups.nEligible];
keepGroups = nEligibleAll >= p.minpergroup;
groupsDropped = groups(~keepGroups);
groups = groups(keepGroups);

if isempty(groups)
    warning('All conditions had < minPerGroup=%d eligible cells. Nothing to label.', p.minpergroup);
    osc_labels = struct('cells',struct([]),'labels',[],'chanName',chanName,'dt_minutes',dt_minutes, ...
                        'droppedGroups', {groupsDropped});
    return;
end

fprintf('Condition filtering: kept %d groups (dropped %d groups with < %d eligible cells).\n', ...
    numel(groups), numel(groupsDropped), p.minpergroup);

% ---------------- sample cells across groups (NOT limited by smallest) ----------------
nGroups = numel(groups);

if isinf(maxCells)
    maxCells = sum([groups.nEligible]); % take everything available
end

nAvail = [groups.nEligible];
totalAvail = sum(nAvail);

maxCellsEff = min(maxCells, totalAvail);

fprintf('Sampling across %d groups (avail total=%d). Requested maxCells=%d => using %d.\n', ...
    nGroups, totalAvail, maxCells, maxCellsEff);

% Allocation ("water filling")
basePer = floor(maxCellsEff / nGroups);
alloc = min(nAvail, basePer * ones(1,nGroups));

remaining = maxCellsEff - sum(alloc);
if remaining > 0
    active = find(nAvail > alloc);
    while remaining > 0 && ~isempty(active)
        order = active(randperm(numel(active)));
        for ii = 1:numel(order)
            g = order(ii);
            if alloc(g) < nAvail(g) && remaining > 0
                alloc(g) = alloc(g) + 1;
                remaining = remaining - 1;
            end
        end
        active = find(nAvail > alloc);
    end
end

fprintf('Per-group allocation (min=%d, max=%d), total=%d\n', min(alloc), max(alloc), sum(alloc));
for g = 1:nGroups
    fprintf('  %2d) %s: avail=%d, take=%d\n', g, char(groups(g).name), nAvail(g), alloc(g));
end

% ---- OPTIONAL: shuffle group order (blocks contiguous) ----
shuffleGroupOrder = true;
if shuffleGroupOrder
    gOrder = randperm(nGroups);
else
    gOrder = 1:nGroups;
end

% ---- Build final cell list (robust column concatenation) ----
cells = struct('xy',{}, 'row',{}, 'cellindex',{}, 'full',{}, 'PtxTx',{}, 'Ptx',{}, 'Tx',{});
cells = cells(:); % force column

for gg = 1:numel(gOrder)
    g = gOrder(gg);

    takeN = alloc(g);
    if takeN <= 0
        continue;
    end

    idx = randperm(groups(g).nEligible);
    take = idx(1:takeN);

    block = groups(g).cells(take);
    block = block(:); % force column

    cells = [cells; block]; %#ok<AGROW>
end

% ---- MUST DEFINE THESE BEFORE THE WHILE LOOP ----
nUse = numel(cells);
labels = nan(nUse,1); % 1=osc, 0=not, -1=unsure, NaN=unlabeled

fprintf('Final: %d cells queued for labeling.\n', nUse);

if nUse == 0
    warning('No cells selected for labeling after allocation.');
    osc_labels = struct();
    osc_labels.cells = struct([]);
    osc_labels.labels = [];
    osc_labels.chanName = chanName;
    osc_labels.dt_minutes = dt_minutes;
    osc_labels.minPerGroup = p.minpergroup;
    osc_labels.droppedGroups = groupsDropped;
    return;
end


% ---------------- figure layout ----------------
hFig = figure('Name','Oscillation triage (condition-aware sampling)', ...
              'NumberTitle','off','Color','w');

set(hFig,'Units','pixels');
pos = get(hFig,'Position');
pos(3) = 1200;
pos(4) = 450;
set(hFig,'Position',pos);

% fixed x-axis across all cells (in hours)
if Tmax_global < 1
    xMaxHours = dt_hours;
else
    xMaxHours = (Tmax_global-1) * dt_hours;
end

k = 1;
while k>=1 && k<=nUse && ishandle(hFig)
    c = cells(k);
    S = dataloc.d{c.xy};
    tr = S.data.(chanName)(c.row,:);

    Tcell = numel(tr);
    t_hours = (0:Tcell-1) * dt_hours;

    figure(hFig);
    clf(hFig);
    ax = axes('Parent',hFig); hold(ax,'on');

    plot(ax, t_hours, tr, '-', 'LineWidth', 2);

    if Tcell >= 3
        trSmooth = movmean(tr, 3, 'omitnan');
        plot(ax, t_hours, trSmooth, ':', 'LineWidth', 1);
    end

    xlim(ax, [0, xMaxHours]);
    set(ax,'Box','on','FontSize',10);

    ttl = sprintf('Cell %d/%d   xy=%d   cellindex=%d\nPtxTx: %s', ...
        k, nUse, c.xy, c.cellindex, char(c.PtxTx));
    title(ax, ttl, 'Interpreter','none');

    xlabel(ax,'time (h)');
    ylabel(ax, chanName);

    if ~isnan(labels(k))
        ensureTxt = "unlabeled";
        if labels(k)==1, ensureTxt="oscillating"; end
        if labels(k)==0, ensureTxt="not oscillating"; end
        if labels(k)==-1, ensureTxt="unsure"; end
        text(ax, 0.02, 0.93, "current label: " + ensureTxt, 'Units','normalized');
    end

    text(ax, 0.02, 0.05, ...
        'Click in figure then: o=osc, n=not, u=unsure, b=back, q=quit', ...
        'Units','normalized','FontSize',9);

    hold(ax,'off');

    key = '';
    while ishandle(hFig) && isempty(key)
        figure(hFig);
        waitforbuttonpress;
        if ~ishandle(hFig), break; end
        key = lower(get(hFig,'CurrentCharacter'));
        if isempty(key) || ~ismember(key, ['o','n','u','b','q'])
            key = '';
        end
    end
    if ~ishandle(hFig), break; end

    switch key
        case 'o', labels(k) = 1;  k = k+1;
        case 'n', labels(k) = 0;  k = k+1;
        case 'u', labels(k) = -1; k = k+1;
        case 'b', if k>1, k = k-1; end
        case 'q', break;
    end
end

if ishandle(hFig)
    close(hFig);
end

osc_labels = struct();
osc_labels.cells = cells;
osc_labels.labels = labels;
osc_labels.chanName = chanName;
osc_labels.dt_minutes = dt_minutes;
osc_labels.minPerGroup = p.minpergroup;
osc_labels.droppedGroups = groupsDropped;

end

% ======================================================================
% Helper: readable Ptx/Tx strings (same as in ct_make_osc_pulse_df)
% ======================================================================
function [cellLine, preStr, txStr, fullStr] = ct__buildTxStrings(tTxList, tpPerHr)
cellLine = " ";
if ~isempty(tTxList)
    isCell = contains({tTxList.dunit}, 'cells')';
    if any(isCell)
        cellLine = string(tTxList(find(isCell,1,'first')).name);
    end
end

tTxList = tTxList(~contains({tTxList.dunit}, 'cells'));

preParts = strings(0,1);
txParts  = strings(0,1);

for iTx = 1:numel(tTxList)
    tunit = string(tTxList(iTx).tunit);
    if ~ismember(tunit, ["h","d","tp"]), continue; end

    switch tunit
        case "h",  hour = tTxList(iTx).time;
        case "d",  hour = tTxList(iTx).time * 24;
        case "tp", hour = round(tTxList(iTx).time / tpPerHr);
    end

    dose  = tTxList(iTx).dose;
    dunit = string(tTxList(iTx).dunit);
    nm    = string(tTxList(iTx).name);

    if isnan(dose) || strlength(strtrim(dunit))==0
        doseStr = nm;
    else
        doseStr = sprintf('%g%s %s', dose, dunit, nm);
    end
    doseStr = strrep(doseStr, '_',' ');

    if hour < 0
        preParts(end+1,1) = doseStr; %#ok<AGROW>
    else
        txParts(end+1,1)  = doseStr; %#ok<AGROW>
    end
end

preStr  = strjoin(preParts, " + ");
txStr   = strjoin(txParts,  " + ");
fullStr = strjoin([preParts; txParts], " + ");
end
