function fig = cycif_plot_scatter(T, markers, group_by, opts)
% cycif_plot_scatter  Pairwise marker scatter matrix, colored by group.
%
% INPUT
%   T         — table from cycif_build_table
%   markers   — {1×M} cell, pass [] to auto-detect
%   group_by  — column in T to color by          (default 'tx1_label')
%   opts      — struct:
%                 .max_pts  max cells to plot (random subsample; default 2000)
%                 .alpha    marker alpha                       (default 0.3)
%                 .mksz     marker size                        (default 4)

if nargin < 3 || isempty(group_by), group_by = 'tx1_label'; end
if nargin < 4, opts = struct(); end
if isempty(markers), markers = i_auto_markers(T); end

max_pts = i_getopt(opts, 'max_pts', 2000);
alpha   = i_getopt(opts, 'alpha',   0.3);
mksz    = i_getopt(opts, 'mksz',    4);

grp_levels = unique(T.(group_by), 'stable');
nGrp = numel(grp_levels);
cmap = lines(nGrp);

% Subsample for speed
if height(T) > max_pts
    idx = randperm(height(T), max_pts);
    T   = T(idx,:);
end

grp_idx = arrayfun(@(i) find(grp_levels == T.(group_by)(i)), (1:height(T))');
Xmat    = table2array(T(:, markers));
nMk     = numel(markers);

fig = figure('Color','w','Name', ['Pairwise scatter | ' group_by], ...
             'Position',[80 80 130*nMk 120*nMk]);

for r = 1:nMk
    for c = 1:nMk
        ax = subplot(nMk, nMk, (r-1)*nMk + c);
        hold(ax,'on');

        if r == c
            % Diagonal: per-group KDE/histogram
            for gi = 1:nGrp
                v = Xmat(grp_idx==gi, r);
                v = v(~isnan(v));
                if numel(v) > 5
                    [f,x] = ksdensity(v);
                    plot(ax, x, f, 'Color', cmap(gi,:), 'LineWidth', 1.2);
                end
            end
        else
            scatter(ax, Xmat(:,c), Xmat(:,r), mksz, cmap(grp_idx,:), ...
                    'filled', 'MarkerFaceAlpha', alpha);
        end

        % Axis labels on outer edges only
        if r == nMk, xlabel(ax, markers{c}, 'Interpreter','none','FontSize',7); end
        if c == 1,   ylabel(ax, markers{r}, 'Interpreter','none','FontSize',7); end
        ax.XTickLabel = {}; ax.YTickLabel = {};
        ax.FontSize = 7;
        hold(ax,'off');
    end
end

% Legend
lh = gobjects(nGrp,1);
for gi = 1:nGrp
    lh(gi) = scatter(NaN,NaN, 20, cmap(gi,:), 'filled', ...
                     'DisplayName', string(grp_levels{gi}));
end
legend(lh, 'Location','eastoutside', 'Interpreter','none','FontSize',8);
sgtitle(fig, ['Pairwise scatter | group: ' group_by], 'FontWeight','bold');
end