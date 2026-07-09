function gate = cycif_flow_gate(varargin)
%CYCIF_FLOW_GATE  Interactively draw a gate; return a replayable struct.
%
% gate = cycif_flow_gate(Name, Value, ...)          % draw on current axes
% gate = cycif_flow_gate(T, Name, Value, ...)       % build plot then draw
%
% NAME-VALUE
%   'Name'    : (required) gate name
%   'Shape'   : 'poly' | 'rect' | 'thresh'            default 'poly'
%   'Parent'  : parent gate name (for hierarchies)    default ''
%   'Axes'    : axes handle to draw on                default gca
%   'X'       : X channel name (required if T given; also stored in gate)
%   'Y'       : Y channel name (2D only)              default ''
%   'XTrans'  : 'linear' | 'log10' | 'asinh'          default 'log10'
%   'YTrans'  : same                                  default 'log10'
%   (When T is given, also accepts: 'GroupBy','Sub','Density','MaxPts',...)
%
% OUTPUT
%   gate : struct with fields
%          .name, .shape, .xchan, .ychan, .xtrans, .ytrans,
%          .coords  (stored in TRANSFORMED axis space),
%          .parent, .created
%
% Notes
%   - Gate drawing requires a single-axes figure; for grouped data, either
%     pass Sub to restrict to one group or use Facet='overlay'.
%   - Storing coords in transformed space means downstream code replays the
%     same transform before inpolygon / threshold tests.

% detect a leading table
if ~isempty(varargin) && istable(varargin{1})
    T    = varargin{1};
    args = varargin(2:end);
else
    T    = [];
    args = varargin;
end

p = inputParser; p.KeepUnmatched = true;
p.addParameter('Name','',    @(x) ischar(x) || isstring(x));
p.addParameter('Shape','poly');
p.addParameter('Parent','');
p.addParameter('Axes',[], @(x) isempty(x) || isgraphics(x));
p.addParameter('X','', @(x) ischar(x) || isstring(x));
p.addParameter('Y','', @(x) ischar(x) || isstring(x));
p.addParameter('XTrans','log10');
p.addParameter('YTrans','log10');
p.parse(args{:}); o = p.Results;

assert(~isempty(char(o.Name)), 'cycif_flow_gate: Name is required.');
xchan = char(o.X);
ychan = char(o.Y);

% If a table was supplied, build the plot first (forwarding unmatched args)
if ~isempty(T)
    assert(~isempty(xchan), 'cycif_flow_gate: X is required when T is given.');
    fwd = {'XTrans', o.XTrans, 'YTrans', o.YTrans};
    um  = p.Unmatched; fn = fieldnames(um);
    for k = 1:numel(fn); fwd(end+1:end+2) = {fn{k}, um.(fn{k})}; end
    if ~isempty(ychan)
        ax = cycif_flow_plot(T, xchan, 'Y', ychan, fwd{:});
    else
        ax = cycif_flow_plot(T, xchan, fwd{:});
    end
    assert(numel(ax) == 1, ['cycif_flow_gate: gate drawing requires a ', ...
        'single axes (got %d). Use Sub to restrict to one group, or ', ...
        'Facet=''overlay''.'], numel(ax));
else
    ax = o.Axes; if isempty(ax); ax = gca; end
end

% Draw the ROI in transformed axes coordinates
title(ax, sprintf('Draw gate: %s  (%s)  — finish with double-click', ...
                  char(o.Name), lower(o.Shape)), 'Interpreter','none');

switch lower(o.Shape)
    case 'poly'
        r = drawpolygon(ax,'Color','r','StripeColor','w');
        wait(r);                                % wait for double-click finish
        coords = r.Position;
    case 'rect'
        r = drawrectangle(ax,'Color','r');
        wait(r);
        pp = r.Position;                        % [x y w h]
        coords = [pp(1) pp(2); pp(1)+pp(3) pp(2)+pp(4)];
    case 'thresh'
        r = drawline(ax,'Color','r');
        wait(r);
        coords = r.Position(1,1);               % x-threshold
    otherwise
        error('cycif_flow_gate: unknown shape ''%s''.', o.Shape);
end

title(ax, sprintf('%s', char(o.Name)), 'Interpreter','none');

gate = struct( ...
    'name',    char(o.Name), ...
    'shape',   lower(o.Shape), ...
    'xchan',   xchan, ...
    'ychan',   ychan, ...
    'xtrans',  o.XTrans, ...
    'ytrans',  o.YTrans, ...
    'coords',  coords, ...
    'parent',  char(o.Parent), ...
    'created', datetime('now'));
end
