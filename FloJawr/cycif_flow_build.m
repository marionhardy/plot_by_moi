function T = cycif_flow_build(dataloc, varargin)
%CYCIF_FLOW_BUILD  Assemble a FlowJo-style per-cell event table from dataloc.
%
% T = cycif_flow_build(dataloc, Name, Value, ...)
%
% INPUT
%   dataloc : dataloc struct (uses .d{xy}.data, .d{xy}.cellindex,
%             .platemapd.pmd, and optionally .movieinfo.tsamp).
%
% NAME-VALUE
%   'XYs'             : XY indices to include (default: all non-empty)
%   'SensorSummaries' : subset of {'end','mean','max','auc'}  (default all)
%   'IncludeCoords'   : include XCoord/YCoord/nArea summaries  (default false)
%   'dt'              : minutes per frame for AUC scaling (default: dataloc)
%
% OUTPUT
%   T : table (one row per cell) with columns:
%       cellID, xy, trackID, row, col, Cell, Gene, pTx_Name, pTx_Conc,
%       pTx_Units, pTx_Time, Tx1_Name, Tx2_Name,  followed by every [N×1]
%       field from .data (stains), and for each [N×nT] sensor one column
%       per requested summary: <sensor>_end / _mean / _max / _auc.
%
% Channel / sensor / stain lists are discovered dynamically from .data;
% nothing is hard-coded.

% -- parse inputs ---------------------------------------------------------
p = inputParser;
p.addParameter('XYs', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('SensorSummaries', {'end','mean','max','auc'}, ...
               @(x) iscellstr(x) || isstring(x));                 %#ok<ISCLSTR>
p.addParameter('IncludeCoords', false, @islogical);
p.addParameter('dt', [], @(x) isempty(x) || isnumeric(x));
p.parse(varargin{:});
o = p.Results;
o.SensorSummaries = cellstr(o.SensorSummaries);

if isempty(o.XYs)
    o.XYs = find(~cellfun('isempty', dataloc.d));
end
if isempty(o.dt)
    if isfield(dataloc,'movieinfo') && isfield(dataloc.movieinfo,'tsamp')
        o.dt = dataloc.movieinfo.tsamp;
    else
        o.dt = 1;
    end
end

coordFields = {'XCoord','YCoord','nArea'};

% -- Pass 1: channel inventory across all selected XYs --------------------
stainSet  = {};
sensorSet = {};
for ii = 1:numel(o.XYs)
    xy = o.XYs(ii);
    if isempty(dataloc.d{xy}) || ~isfield(dataloc.d{xy},'data'); continue; end
    D  = dataloc.d{xy}.data;
    fn = fieldnames(D);
    for f = 1:numel(fn)
        v = D.(fn{f});
        if ~isnumeric(v) || isempty(v); continue; end
        if size(v,2) == 1                            % [N×1] → stain
            stainSet{end+1}  = fn{f};                            %#ok<AGROW>
        elseif size(v,2) > 1                         % [N×T] → sensor
            if o.IncludeCoords || ~ismember(fn{f}, coordFields)
                sensorSet{end+1} = fn{f};                        %#ok<AGROW>
            end
        end
    end
end
stainSet  = unique(stainSet);
sensorSet = unique(sensorSet);

% -- Pass 2: per-well row assembly ----------------------------------------
tiles = cell(numel(o.XYs),1);
for ii = 1:numel(o.XYs)
    xy = o.XYs(ii);
    if isempty(dataloc.d{xy}) || ~isfield(dataloc.d{xy},'data'); continue; end
    D = dataloc.d{xy}.data;

    % cell count N
    fn = fieldnames(D);
    N = 0;
    for f = 1:numel(fn)
        if isnumeric(D.(fn{f})) && ~isempty(D.(fn{f}))
            N = size(D.(fn{f}),1); break;
        end
    end
    if N == 0; continue; end

    % trackID
    if isfield(dataloc.d{xy},'cellindex') && ~isempty(dataloc.d{xy}.cellindex)
        trackID = dataloc.d{xy}.cellindex(:);
        if numel(trackID) ~= N; trackID = (1:N).'; end
    else
        trackID = (1:N).';
    end

    t = table();
    t.xy      = repmat(xy, N, 1);
    t.trackID = trackID;

    % platemap metadata
    [wrow, wcol, meta] = local_pm_lookup(dataloc, xy);
    t.row = repmat(wrow, N, 1);
    t.col = repmat(wcol, N, 1);
    mf = fieldnames(meta);
    for f = 1:numel(mf)
        val = meta.(mf{f});
        if isnumeric(val)
            t.(mf{f}) = repmat(val, N, 1);
        else
            t.(mf{f}) = repmat(string(val), N, 1);
        end
    end

    % stains: [N×1] fields, preserve name verbatim
    for f = 1:numel(stainSet)
        nm = stainSet{f};
        if isfield(D,nm) && size(D.(nm),1)==N && size(D.(nm),2)==1
            t.(nm) = D.(nm);
        else
            t.(nm) = nan(N,1);
        end
    end

    % sensors: [N×T] summaries
    for f = 1:numel(sensorSet)
        nm = sensorSet{f};
        haveIt = isfield(D,nm) && size(D.(nm),1)==N && size(D.(nm),2)>1;
        for s = 1:numel(o.SensorSummaries)
            col = [nm '_' o.SensorSummaries{s}];
            if haveIt
                t.(col) = local_summary(D.(nm), o.SensorSummaries{s}, o.dt);
            else
                t.(col) = nan(N,1);
            end
        end
    end

    tiles{ii} = t;
end

T = vertcat(tiles{~cellfun('isempty', tiles)});

% add global cellID as first column
T.cellID = (1:height(T)).';
T = movevars(T, 'cellID', 'Before', 1);

% cast categorical metadata
catCols = intersect(T.Properties.VariableNames, ...
    {'Cell','Gene','pTx_Name','pTx_Units','pTx_Time','Tx1_Name','Tx2_Name','row'});
for c = 1:numel(catCols)
    T.(catCols{c}) = categorical(T.(catCols{c}));
end

end % cycif_flow_build

% -------------------------------------------------------------------------
function [wrow, wcol, meta] = local_pm_lookup(dataloc, xy)
% Map xy → (row letter, col index) and pull per-well metadata fields.
wrow = "?"; wcol = NaN; meta = struct();
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'pmd'); return; end
pm   = dataloc.platemapd.pmd;
pmxy = pm.xy;
rowL = 'ABCDEFGH';
r = 0; c = 0;
for rr = 1:size(pmxy,1)
    for cc = 1:size(pmxy,2)
        v = pmxy{rr,cc};
        if isnumeric(v) && any(v == xy); r = rr; c = cc; break; end
    end
    if r; break; end
end
if r == 0; return; end
wrow = string(rowL(min(r,numel(rowL))));
wcol = c;

getL1 = @(F) local_cell2str(pm.(F){r,c,1});                        % layer-1 cell
if isfield(pm,'Cell');  meta.Cell  = getL1('Cell'); end
if isfield(pm,'Gene');  meta.Gene  = getL1('Gene'); end

if isfield(pm,'pTx') && size(pm.pTx,3) >= 5
    meta.pTx_Name  = local_cell2str(pm.pTx{r,c,1});
    meta.pTx_Conc  = str2double(local_cell2str(pm.pTx{r,c,2}));
    meta.pTx_Units = local_cell2str(pm.pTx{r,c,3});
    meta.pTx_Time  = local_cell2str(pm.pTx{r,c,4});
end
if isfield(pm,'Tx1')
    if size(pm.Tx1,3) >= 1
        meta.Tx1_Name = local_cell2str(pm.Tx1{r,c,1});
    end
    if size(pm.Tx1,3) >= 6
        meta.Tx2_Name = local_cell2str(pm.Tx1{r,c,6});
    end
end
end

function s = local_cell2str(v)
if isempty(v);                             s = "";          return; end
if ischar(v) || isstring(v);               s = string(v);   return; end
if isnumeric(v) && isscalar(v);            s = string(v);   return; end
s = string(v);
end

function y = local_summary(M, kind, dt)
% Compute a per-cell summary from an [N×T] time series, NaN-robust.
switch lower(kind)
    case 'end'
        y = nan(size(M,1),1);
        for i = 1:size(M,1)
            j = find(isfinite(M(i,:)), 1, 'last');
            if ~isempty(j); y(i) = M(i,j); end
        end
    case 'mean'
        y = mean(M, 2, 'omitnan');
    case 'max'
        y = max(M, [], 2, 'omitnan');
    case 'auc'
        y = sum(M, 2, 'omitnan') .* dt;
    otherwise
        error('cycif_flow_build:badSummary','Unknown summary ''%s''.', kind);
end
end
