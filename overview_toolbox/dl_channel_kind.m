function kind = dl_channel_kind(dataloc, channel, varargin)
% DL_CHANNEL_KIND  Classify a data channel as time-series or endpoint.
%
% kind = dl_channel_kind(dataloc, channel)
% kind = dl_channel_kind(dataloc, channel, 'xy', xyIdx)
%
% INPUT  (paired values)
%   dataloc : dataloc struct (v5). Reads dataloc.d{xy}.data.(channel).
%   channel : char/string channel name (e.g. 'HYLIGHT', 'iNOS_nuc').
%   'xy'    : [scalar] which non-empty xy to probe (default: first non-empty).
%
% OUTPUT
%   kind : char, one of
%          'trace'    -> [N x nT], nT > 1  (live-cell time series)
%          'endpoint' -> [N x 1]           (CycIF endpoint intensity, etc.)
%          'missing'  -> channel not present in the probed xy
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): v4's plot_by_MHY assumed every
% d{xy}.data field was an [N x nT] trace and blindly indexed
% (:,firsttp:tracklength). The v5 dataloc mixes [N x 249] live traces with
% [N x 1] CycIF endpoint scalars (ANTIBODY_nuc/_cyt/_cell) in the SAME data
% struct. Passing an endpoint channel to a time-series plot in v4 produced
% an opaque indexing error. This helper lets the v5 plotters detect kind up
% front and either route endpoint channels to endpoint-appropriate plots or
% error clearly. Classification is by column count, matching the schema:
% traces have nT timepoints, endpoints are single-column.
% -------------------------------------------------------------------------

ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'xy', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(ip, varargin{:});
xyIdx = ip.Results.xy;

if ~ischar(channel) && ~isstring(channel)
    error('dl_channel_kind:BadChannel','channel must be a char/string.');
end
channel = char(channel);

% Resolve a probe xy: caller-specified, else first non-empty with data.
if isempty(xyIdx)
    xyIdx = find(cellfun(@(c) ~isempty(c) && isfield(c,'data'), dataloc.d), 1);
    if isempty(xyIdx)
        error('dl_channel_kind:NoData','dataloc.d has no non-empty xy with a data field.');
    end
end

if xyIdx > numel(dataloc.d) || isempty(dataloc.d{xyIdx}) || ~isfield(dataloc.d{xyIdx},'data') ...
        || ~isfield(dataloc.d{xyIdx}.data, channel)
    kind = 'missing';
    return
end

nCols = size(dataloc.d{xyIdx}.data.(channel), 2);
if nCols > 1
    kind = 'trace';
else
    kind = 'endpoint';
end
end
