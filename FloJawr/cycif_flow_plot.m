function [ax, h] = cycif_flow_plot(T, xchan, varargin)
%CYCIF_FLOW_PLOT  Biaxial scatter or 1D histogram of a flow-style table.
%
% [ax, h] = cycif_flow_plot(T, xchan, Name, Value, ...)
%
% INPUT
%   T      : flow table from cycif_flow_build
%   xchan  : name of X channel (column in T)
%
% NAME-VALUE
%   'Y'         : name of Y column (empty → 1D histogram of X)   default ''
%   'GroupBy'   : column name to split by colour or facet         default ''
%   'Facet'     : 'overlay' | 'tile' | 'auto'                     default 'auto'
%   'XTrans'    : 'linear' | 'log10' | 'asinh'                    default 'log10'
%   'YTrans'    : 'linear' | 'log10' | 'asinh'                    default 'log10'
%   'Density'   : colour single-group scatter by density          default true
%   'Sub'       : logical N×1 mask, pre-filter rows               default []
%   'Gate'      : gate struct / struct array to overlay           default []
%   'NBins'     : histogram bins                                  default 100
%   'Alpha'     : marker face alpha                               default 0.25
%   'MaxPts'    : max points plotted per group                    default 2e4
%   'Parent'    : figure handle to draw into                      default []
%
% OUTPUT
%   ax : axes handle(s)
%   h  : graphics handle(s)  (scatter or histogram)
%
% Notes
%   Transforms are applied per axis. Non-positive values are floored at
%   eps under log10. Density colouring only runs for single-group 2D plots.

% -- parse -----------------------------------------------------------------
p = inputParser;
p.addParameter('Y','',       @(x) ischar(x) || isstring(x));
p.addParameter('GroupBy','', @(x) ischar(x) || isstring(x));
p.addParameter('Facet','auto', @(x) any(strcmpi(x,{'auto','overlay','tile'})));
p.addParameter('XTrans','log10');
p.addParameter('YTrans','log10');
p.addParameter('Density',true, @islogical);
p.addParameter('Sub', [],    @(x) isempty(x) || islogical(x));
p.addParameter('Gate', [],   @(x) isempty(x) || isstruct(x));
p.addParameter('NBins',100,  @isnumeric);
p.addParameter('Alpha',0.25, @isnumeric);
p.addParameter('MaxPts',2e4, @isnumeric);
p.addParameter('Parent',[]);
p.parse(varargin{:}); o = p.Results;

ychan  = char(o.Y);
grpcol = char(o.GroupBy);

% subset & validate
if ~isempty(o.Sub); T = T(o.Sub,:); end
assert(height(T) > 0, 'cycif_flow_plot: no rows to plot.');
assert(ismember(xchan, T.Properties.VariableNames), ...
    'cycif_flow_plot: column "%s" not in table.', xchan);

do2D = ~isempty(ychan);
if do2D
    assert(ismember(ychan, T.Properties.VariableNames), ...
        'cycif_flow_plot: column "%s" not in table.', ychan);
end

% transform
xt = local_trans(T.(xchan), o.XTrans);
yt = [];
if do2D; yt = local_trans(T.(ychan), o.YTrans); end

% grouping
if isempty(grpcol)
    g = ones(height(T),1); gL = "all";
else
    assert(ismember(grpcol, T.Properties.VariableNames), ...
        'cycif_flow_plot: GroupBy column "%s" not in table.', grpcol);
    [g, gL] = findgroups(T.(grpcol));
    gL = string(gL);
end
nG = numel(gL);

% facet mode
facet = lower(o.Facet);
if strcmp(facet,'auto')
    if do2D && nG > 1; facet = 'tile'; else; facet = 'overlay'; end
end

% figure
if isempty(o.Parent); fig = gcf; else; fig = ancestor(o.Parent,'figure'); end
clf(fig);
cmap = lines(max(nG,1));

if strcmp(facet,'tile') && nG > 1
    tl = tiledlayout(fig,'flow','TileSpacing','compact','Padding','compact');
    ax = gobjects(nG,1); h = gobjects(nG,1);
    for k = 1:nG
        ax(k) = nexttile(tl);
        m = g == k;
        h(k) = local_draw(ax(k), xt(m), sub_if(yt,m), cmap(k,:), o, do2D, false);
        title(ax(k), gL(k), 'Interpreter','none');
        local_format(ax(k), xchan, ychan, o.XTrans, o.YTrans);
    end
else
    ax = axes(fig); hold(ax,'on');
    h = gobjects(nG,1);
    useDens = o.Density && nG == 1 && do2D;
    for k = 1:nG
        m = g == k;
        h(k) = local_draw(ax, xt(m), sub_if(yt,m), cmap(k,:), o, do2D, useDens);
    end
    if nG > 1
        legend(ax, gL, 'Interpreter','none','Location','best');
    end
    local_format(ax, xchan, ychan, o.XTrans, o.YTrans);
end

% overlay gates
if ~isempty(o.Gate)
    for a = 1:numel(ax)
        for k = 1:numel(o.Gate)
            local_overlay(ax(a), o.Gate(k));
        end
    end
end
end

% -------------------------------------------------------------------------
function y = sub_if(yt, m)
if isempty(yt); y = []; else; y = yt(m); end
end

function xt = local_trans(x, tr)
switch lower(tr)
    case 'linear'; xt = x;
    case 'log10';  xt = log10(max(x, eps));
    case 'asinh';  xt = asinh(x/150);                 % biexp-like, cofactor 150
    otherwise;     error('cycif_flow:badTrans','Unknown transform ''%s''.', tr);
end
end

function h = local_draw(ax, xt, yt, col, o, do2D, useDens)
% Subsample if large
n = numel(xt);
if n > o.MaxPts
    idx = randperm(n, o.MaxPts);
    xt  = xt(idx);
    if do2D; yt = yt(idx); end
end
if do2D
    if useDens && n > 200
        c = local_density(xt, yt);
        h = scatter(ax, xt, yt, 6, c, 'filled', 'MarkerFaceAlpha', o.Alpha);
        colormap(ax, parula);
    else
        h = scatter(ax, xt, yt, 6, col, 'filled', 'MarkerFaceAlpha', o.Alpha);
    end
else
    xok = xt(isfinite(xt));
    if isempty(xok)
        h = gobjects(1); return;
    end
    edges = linspace(min(xok), max(xok), o.NBins+1);
    h = histogram(ax, xt, edges, 'Normalization','probability', ...
                  'FaceColor',col, 'EdgeColor','none', 'FaceAlpha',0.5);
end
end

function c = local_density(xt, yt)
% 2D histogram density, smoothed, mapped back to each point.
c  = zeros(size(xt));
ok = isfinite(xt) & isfinite(yt);
if nnz(ok) < 50; return; end
n  = 80;
xe = linspace(min(xt(ok)), max(xt(ok)), n+1);
ye = linspace(min(yt(ok)), max(yt(ok)), n+1);
ix = discretize(xt(ok), xe);
iy = discretize(yt(ok), ye);
vv = ~isnan(ix) & ~isnan(iy);
H  = accumarray([iy(vv) ix(vv)], 1, [n n]);
if exist('imgaussfilt','file') == 2
    H = imgaussfilt(H, 1.5);
end
okIdx = find(ok);
okIdx = okIdx(vv);
c(okIdx) = H(sub2ind([n n], iy(vv), ix(vv)));
end

function local_format(ax, xchan, ychan, xt, yt)
grid(ax,'on'); box(ax,'on');
xlabel(ax, sprintf('%s  [%s]', xchan, xt), 'Interpreter','none');
if ~isempty(ychan)
    ylabel(ax, sprintf('%s  [%s]', ychan, yt), 'Interpreter','none');
else
    ylabel(ax, 'probability');
end
end

function local_overlay(ax, G)
% Draw a gate boundary on an axes (coords are in TRANSFORMED space).
hold(ax,'on');
switch lower(G.shape)
    case 'rect'
        c = G.coords;          % 2×2: [x1 y1; x2 y2]
        plot(ax, [c(1,1) c(2,1) c(2,1) c(1,1) c(1,1)], ...
                 [c(1,2) c(1,2) c(2,2) c(2,2) c(1,2)], ...
                 'r-', 'LineWidth', 1.25);
    case 'poly'
        V = G.coords;
        plot(ax, [V(:,1); V(1,1)], [V(:,2); V(1,2)], ...
                 'r-', 'LineWidth', 1.25);
    case 'thresh'
        v  = G.coords(1);
        yl = ylim(ax);
        plot(ax, [v v], yl, 'r--', 'LineWidth', 1.25);
end
end
