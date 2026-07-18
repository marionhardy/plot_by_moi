function p = dl_plotparse(plotby, dataloc, extra, varargin)
% DL_PLOTPARSE  Shared value-pair parameter parser for the v5 plot suite.
%
% p = dl_plotparse(plotby, dataloc, extra, varargin)
%
% Every renderer (plot_mean, plot_heatmap, plot_sorted_stacks) parses its
% Name/Value inputs through this so the contract is identical and defined
% ONCE. All inputs are paired values (per lab convention).
%
% INPUT
%   plotby  : 'treatment' | 'celltype' | 'custom' (validated + normalized).
%   dataloc : v5 dataloc (read for movieinfo.tsamp -> p.tktm).
%   extra   : struct of renderer-specific addParameter specs, as
%             extra.<name> = {default, validator}. e.g. for stacks:
%             extra.ntracks = {5, @isscalar}.
%   varargin: the caller's Name/Value pairs.
%
% OUTPUT
%   p : parsed struct with all shared + extra fields, plus derived:
%       p.tktm (frames per hour), p.plotby (normalized), and cell-wrapped
%       channel/subset/exclude/ymn/ymx.
%
% -------------------------------------------------------------------------
% V5 CHANGE: v4 declared ~40 addParameter lines inline; the first v5 port
% copied them into each plot_ file (drift risk). Centralized here. 'save'
% is new (restores SVG export capability). 'plottype' is intentionally
% ABSENT -- the function identity is the plottype now.
% -------------------------------------------------------------------------

ip = inputParser;
ip.KeepUnmatched = true;      % tolerate stray v4 params (e.g. 'plottype')
ip.CaseSensitive = false;

isNNScalar = @(x) isnumeric(x) && isscalar(x) && x >= 0;

% ---- shared params (superset honoring all original plot_by_MHY calls) ----
addParameter(ip, 'channel', {'ERKTR'});
addParameter(ip, 'aftertreatment', 1, isNNScalar);
addParameter(ip, 'closefigs', true, @islogical);   % V5: restore v4 lifecycle
addParameter(ip, 'combinexys', true, @islogical);
addParameter(ip, 'exclude', [], @(x) isempty(x)||isnumeric(x)||iscell(x)||ischar(x));
addParameter(ip, 'facetby', {}, @(x) iscell(x)||isstring(x));
addParameter(ip, 'font_size', 8, isNNScalar);
addParameter(ip, 'groupby', {}, @(x) iscell(x)||isstring(x));
addParameter(ip, 'ifchan', 'CY5_Nuc', @ischar);
addParameter(ip, 'ncells', [], @(x) isempty(x)||isscalar(x));
addParameter(ip, 'nogene', false, @(x) islogical(x)||isnumeric(x));
addParameter(ip, 'plotfromzero', false, @islogical);
addParameter(ip, 'save', true, @islogical);   % V5: default true (save on by default)
addParameter(ip, 'smooth', [], @(x) isempty(x)||isvector(x));
addParameter(ip, 'specsubset', false, @islogical);
addParameter(ip, 'standardizeplots', true, @islogical);   % V5: gates dl_plotstandardize
addParameter(ip, 'alignzero', false, @islogical);          % V5: align yyaxis zero-lines
addParameter(ip, 'subset', [], @(x) isempty(x)||isnumeric(x)||iscell(x)||ischar(x));
addParameter(ip, 'tbeforetx', [], @(x) isempty(x)||isscalar(x));
addParameter(ip, 'tmaxaftertx', [], @(x) isempty(x)||isscalar(x));
addParameter(ip, 'tx_order', 1, isNNScalar);
addParameter(ip, 'ymn', [], @(x) isnumeric(x)||iscell(x));
addParameter(ip, 'ymx', [], @(x) isnumeric(x)||iscell(x));
addParameter(ip, 'zerohrtx', [], @(x) isempty(x)||ischar(x));

% ---- renderer-specific params ----
if ~isempty(extra)
    fn = fieldnames(extra);
    for k = 1:numel(fn)
        spec = extra.(fn{k});
        if numel(spec) >= 2; addParameter(ip, fn{k}, spec{1}, spec{2});
        else;                addParameter(ip, fn{k}, spec{1});
        end
    end
end

parse(ip, varargin{:});
p = ip.Results;

% ---- normalize plotby (custom triggered by groupby/facetby) ----
if ~isempty(p.groupby) || ~isempty(p.facetby); plotby = 'custom'; end
switch lower(plotby)
    case {'treatment','treatments','tx','txs'}; p.plotby = 'treatment';
    case {'cell','celltype','celltypes'};        p.plotby = 'celltype';
    case {'custom','group','grouped'};           p.plotby = 'custom';
    otherwise
        error('dl_plotparse:BadPlotby','plotby must be treatment, celltype, or custom.');
end

% ---- normalize channel to canonical form: {spec1, spec2, ...} where each
%      spec is a cellstr of channels to overlay in one figure. This restores
%      the two-level v4 semantics ('p.channel{k}' is a plot-spec, its
%      elements are overlaid channels) that a flat cell-wrap would collapse.
%      V5 NOTE: the clean-rewrite regression that broke {{'A','B'}} lived
%      here -- a bare {~iscell} wrap left the nested spec un-normalized.
if ischar(p.channel) || isstring(p.channel)
    p.channel = {{char(p.channel)}};                 % 'A'          -> {{'A'}}
elseif iscell(p.channel) && all(cellfun(@(x) ischar(x)||isstring(x), p.channel))
    p.channel = {cellfun(@char, p.channel, 'Un', 0)}; % {'A','B'}    -> {{'A','B'}}
elseif iscell(p.channel)
    p.channel = cellfun(@(s) local_wrap_spec(s), p.channel, 'Un', 0); % {{'A','B'},{'C'}} -> unchanged
end
if ~isempty(p.subset)  && ~iscell(p.subset);  p.subset  = {p.subset};  end
if ~isempty(p.exclude) && ~iscell(p.exclude); p.exclude = {p.exclude}; end
if ~isempty(p.ymn) && ~iscell(p.ymn); p.ymn = {p.ymn}; end
if ~isempty(p.ymx) && ~iscell(p.ymx); p.ymx = {p.ymx}; end

% ---- derive frames-per-hour + convert hour-based windows to frames ----
if isfield(dataloc,'movieinfo') && ~isempty(dataloc.movieinfo.tsamp)
    looptime = dataloc.movieinfo.tsamp;
else
    looptime = 6; warning('dl_plotparse:NoTsamp','No movieinfo.tsamp; assuming 6 min/frame.');
end
p.tktm = 60/looptime;
if ~isempty(p.tmaxaftertx); p.tmaxaftertx = floor(p.tmaxaftertx * p.tktm); end
if ~isempty(p.tbeforetx);   p.tbeforetx   = floor(p.tbeforetx   * p.tktm); end
end

function spec = local_wrap_spec(s)
% Ensure one channel spec is a cellstr (single 'A' -> {'A'}).
if ischar(s) || isstring(s); spec = {char(s)}; else; spec = s; end
end