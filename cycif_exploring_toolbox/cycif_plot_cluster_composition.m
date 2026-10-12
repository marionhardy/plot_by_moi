function [P, counts, h] = cycif_plot_cluster_composition(membership, varargin)
% cycif_plot_cluster_composition  Stacked bars of a label's proportion per cluster.
%
%   One bar per cluster, each normalized to 100%, stacked by the categories of
%   a chosen label column (default 'tx2_label'). Returns the proportion matrix
%   so the numbers are traceable, not just drawn.
%
% Usage:
%   [memb, Data] = cycif_cluster_membership(Tf, 'Data', Data);
%   cycif_plot_cluster_composition(memb);                       % tx2_label
%   cycif_plot_cluster_composition(memb, 'label_var','tx1_label');
%   [P,counts] = cycif_plot_cluster_composition(memb, 'order','abundance');
%
% INPUT
%   membership -- table with a 'cluster' column and the label column.
%
% NAME-VALUE
%   'label_var'  char   column to break each cluster down by   (default 'tx2_label')
%   'order'      char   cluster x-order: 'sorted' (default) | 'abundance' (most
%                       cells first) | 'none' (as encountered)
%   'cmap'       Nx3    colors for label categories            (default lines/parula)
%   'font_size'  double                                        (default 9)
%   'ax'         handle plot into existing axes                (default new figure)
%
% OUTPUT
%   P       -- [nClusters x nLabels] proportion matrix (rows sum to 1)
%   counts  -- [nClusters x nLabels] raw cell counts
%   h       -- handle to the bar objects

ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'label_var', 'tx2_label', @(x) ischar(x)||isstring(x));
addParameter(ip, 'order',     'sorted',    @(x) any(strcmpi(x,{'sorted','abundance','none'})));
addParameter(ip, 'cmap',      [],          @(x) isempty(x)||(isnumeric(x)&&size(x,2)==3));
addParameter(ip, 'font_size', 9,           @(x) isnumeric(x)&&isscalar(x));
addParameter(ip, 'ax',        [],          @(x) isempty(x)||isgraphics(x));
parse(ip, varargin{:}); p = ip.Results;
labelVar = char(p.label_var);

% --- validate columns ---
assert(ismember('cluster', membership.Properties.VariableNames), ...
    'membership has no ''cluster'' column.');
assert(ismember(labelVar, membership.Properties.VariableNames), ...
    'membership has no ''%s'' column.', labelVar);

% --- categorical axes ---
clu = categorical(string(membership.cluster));
lab = categorical(string(membership.(labelVar)));

cluLevels = categories(clu);          % cluster order (as-encountered via categories)
labLevels = categories(lab);          % label categories = stack segments

% --- count matrix [cluster x label] ---
counts = zeros(numel(cluLevels), numel(labLevels));
for i = 1:numel(cluLevels)
    inClu = clu == cluLevels{i};
    for j = 1:numel(labLevels)
        counts(i,j) = sum(inClu & lab == labLevels{j});
    end
end

% --- cluster ordering ---
switch lower(p.order)
    case 'abundance'
        [~, ord] = sort(sum(counts,2), 'descend');
    case 'sorted'
        % numeric-aware sort if clusters are numbers, else lexical
        num = str2double(cluLevels);
        if all(~isnan(num)), [~,ord] = sort(num); else, [~,ord] = sort(string(cluLevels)); end
    otherwise
        ord = (1:numel(cluLevels))';
end
cluLevels = cluLevels(ord);
counts    = counts(ord,:);

% --- proportions within each cluster (rows sum to 1) ---
rowTot = sum(counts,2);
rowTot(rowTot==0) = 1;                 % guard empty clusters
P = counts ./ rowTot;

% --- colors ---
if isempty(p.cmap)
    nL = numel(labLevels);
    if nL <= 7, cmap = lines(nL); else, cmap = parula(nL); end
else
    cmap = p.cmap;
end

% --- plot ---
if isempty(p.ax)
    figure('Color','w','Units','inches','Position',[1 1 max(5,0.6*numel(cluLevels)+2) 4]);
    ax = axes;
else
    ax = p.ax;
end

h = bar(ax, P*100, 'stacked', 'EdgeColor','none');
for j = 1:numel(h), h(j).FaceColor = cmap(min(j,size(cmap,1)),:); end

set(ax, 'XTick', 1:numel(cluLevels), 'XTickLabel', cluLevels);
xlabel(ax, 'Cluster');
ylabel(ax, '% of cluster');
ylim(ax, [0 100]);
xlim(ax, [0.5 numel(cluLevels)+0.5]);
legend(ax, h, strrep(labLevels,'_',' '), 'Location','eastoutside', 'Interpreter','none');
title(ax, sprintf('%s composition per cluster', strrep(labelVar,'_',' ')), 'Interpreter','none');
set(ax, 'FontSize', p.font_size, 'TickLabelInterpreter','none', 'Box','off');
end