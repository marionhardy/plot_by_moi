function [P, counts, h] = plot_stacked_composition(T, varargin)
% plot_stacked_composition  Stacked proportion bars from any table.
%
%   One bar per 'y' category, each normalized to 100%, split into segments
%   colored by 'x'. Works on any table with two categorical-izable columns.
%
%   Reads as: "composition of X within each Y".
%     y = the bars        (e.g. 'cluster'  -> one bar per cluster)
%     x = the fill/stacks (e.g. 'tx2_label'-> each bar split by treatment)
%
% Usage:
%   plot_stacked_composition(Tf, 'x','tx2_label', 'y','cluster');
%   plot_stacked_composition(Tf, 'x','cluster',   'y','tx2_label');   % flipped
%   [P,counts] = plot_stacked_composition(Tf, 'x','tx2_label','y','cluster', ...
%                                            'norm','count');
%
% REQUIRED NAME-VALUE
%   'x'   char   column whose categories become the stacked segments (fill)
%   'y'   char   column whose categories become the bars
%
% OPTIONAL NAME-VALUE
%   'norm'       'prop' (default, each bar=100%) | 'count' (raw cell counts)
%   'order'      bar order: 'sorted' (default) | 'abundance' | 'none'
%   'cmap'       Nx3 colors for the x (segment) categories
%   'font_size'  double                                          (default 9)
%   'ax'         axes handle to draw into                        (default new fig)
%
% OUTPUT
%   P       [nBars x nSeg] proportion (or count if norm='count') matrix
%   counts  [nBars x nSeg] raw counts
%   h       bar object handles

ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'x', '', @(s) ischar(s)||isstring(s));
addParameter(ip, 'y', '', @(s) ischar(s)||isstring(s));
addParameter(ip, 'norm', 'prop', @(s) any(strcmpi(s,{'prop','count'})));
addParameter(ip, 'order', 'sorted', @(s) any(strcmpi(s,{'sorted','abundance','none'})));
addParameter(ip, 'cmap', [], @(x) isempty(x)||(isnumeric(x)&&size(x,2)==3));
addParameter(ip, 'font_size', 9, @(x) isnumeric(x)&&isscalar(x));
addParameter(ip, 'ax', [], @(x) isempty(x)||isgraphics(x));
parse(ip, varargin{:}); p = ip.Results;

xVar = char(p.x);  yVar = char(p.y);
assert(~isempty(xVar) && ~isempty(yVar), ...
    'Both ''x'' (segments) and ''y'' (bars) must be given.');
vn = T.Properties.VariableNames;
assert(ismember(xVar, vn), 'Table has no column ''%s''.', xVar);
assert(ismember(yVar, vn), 'Table has no column ''%s''.', yVar);

% categorical axes
xc = categorical(string(T.(xVar)));   % segments
yc = categorical(string(T.(yVar)));   % bars
segLevels = categories(xc);
barLevels = categories(yc);

% count matrix [bar x segment]
counts = zeros(numel(barLevels), numel(segLevels));
for i = 1:numel(barLevels)
    inBar = yc == barLevels{i};
    for j = 1:numel(segLevels)
        counts(i,j) = sum(inBar & xc == segLevels{j});
    end
end

% bar ordering
switch lower(p.order)
    case 'abundance'
        [~, ord] = sort(sum(counts,2), 'descend');
    case 'sorted'
        num = str2double(barLevels);
        if all(~isnan(num)), [~,ord] = sort(num); else, [~,ord] = sort(string(barLevels)); end
    otherwise
        ord = (1:numel(barLevels))';
end
barLevels = barLevels(ord);
counts    = counts(ord,:);

% normalization
if strcmpi(p.norm,'prop')
    rowTot = sum(counts,2); rowTot(rowTot==0) = 1;
    P = counts ./ rowTot;
    Yplot = P*100; yLab = sprintf('%% of %s', strrep(yVar,'_',' '));
    yTop  = 100;
else
    P = counts;
    Yplot = counts; yLab = 'cell count';
    yTop  = [];
end

% colors for segments
if isempty(p.cmap)
    nS = numel(segLevels);
    if nS <= 7, cmap = lines(nS); else, cmap = parula(nS); end
else
    cmap = p.cmap;
end

% draw
if isempty(p.ax)
    figure('Color','w','Units','inches', ...
        'Position',[1 1 max(5,0.6*numel(barLevels)+2.5) 4]);
    ax = axes;
else
    ax = p.ax;
end

h = bar(ax, Yplot, 'stacked', 'EdgeColor','none');
for j = 1:numel(h), h(j).FaceColor = cmap(min(j,size(cmap,1)),:); end

set(ax, 'XTick', 1:numel(barLevels), 'XTickLabel', barLevels);
xlabel(ax, strrep(yVar,'_',' '));
ylabel(ax, yLab);
if ~isempty(yTop), ylim(ax, [0 yTop]); end
xlim(ax, [0.5 numel(barLevels)+0.5]);
legend(ax, h, strrep(segLevels,'_',' '), 'Location','eastoutside', 'Interpreter','none');
title(ax, sprintf('%s composition per %s', strrep(xVar,'_',' '), strrep(yVar,'_',' ')), ...
    'Interpreter','none');
set(ax, 'FontSize', p.font_size, 'TickLabelInterpreter','none', 'Box','off');
end