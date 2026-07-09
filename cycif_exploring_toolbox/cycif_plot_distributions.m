function fig = cycif_plot_distributions(T, varargin)
% cycif_plot_distributions  Per-marker overlaid KDE histograms, grouped and
%                            optionally faceted by treatment variables.
%
% INPUT
%   T  — table from cycif_build_table
%
% NAME-VALUE OPTIONS
%   'markers'    {1xM} cell  marker names to plot; [] = auto-detect
%   'group_by'   string      T column to color/group by    (default 'tx1_label')
%   'facet_by'   string      T column to facet into rows   (default '' = no facet)
%   'n_bins'     scalar      histogram bins                (default 40)
%   'norm'       string      'probability'|'pdf'|'count'   (default 'probability')
%   'show_hist'  logical     bar histogram under KDE       (default false)
%   'alpha'      scalar      KDE fill alpha                (default 0.05)
%   'lw'         scalar      KDE line width                (default 1.5)
%   'colormap'   string      'bright'|'parula'|any MATLAB map (default 'bright')
%   'use_log2'   logical     log2(max(x,1)) before plot    (default true)
%   'xlabel'     string      x-axis label override         (default auto)
%   'fig_pos'    [x y w h]   figure position               (default auto)
%   'save_path'  string      full path to save SVG; '' = no save
%
% OUTPUT
%   fig  — figure handle
%
% USAGE
%   cycif_plot_distributions(T);
%   cycif_plot_distributions(T, 'group_by','tx1_label', 'facet_by','ptx_pool_label');
%   cycif_plot_distributions(T, 'use_log2',false, 'colormap','parula');

p = inputParser();
addParameter(p, 'markers',   []);
addParameter(p, 'group_by',  'tx1_label');
addParameter(p, 'facet_by',  '');
addParameter(p, 'n_bins',    40);
addParameter(p, 'norm',      'probability');
addParameter(p, 'show_hist', false);
addParameter(p, 'alpha',     0.05);
addParameter(p, 'lw',        1.5);
addParameter(p, 'colormap',  'bright');
addParameter(p, 'use_log2',  true);
addParameter(p, 'xlabel',    '');
addParameter(p, 'fig_pos',   []);
addParameter(p, 'save_path', '');
parse(p, varargin{:});

markers   = p.Results.markers;
group_by  = p.Results.group_by;
facet_by  = p.Results.facet_by;
n_bins    = p.Results.n_bins;
norm      = p.Results.norm;
show_hist = p.Results.show_hist;
alpha     = p.Results.alpha;
lw        = p.Results.lw;
palette   = p.Results.colormap;
use_log2  = p.Results.use_log2;
fig_pos   = p.Results.fig_pos;

if ~isempty(p.Results.xlabel)
    xlabel_str = p.Results.xlabel;
elseif use_log2
    xlabel_str = 'log_2 intensity';
else
    xlabel_str = 'raw intensity (a.u.)';
end

if isempty(markers), markers = i_auto_markers(T); end
nMk = numel(markers);

grp_levels = unique(T.(group_by), 'stable');
nGrp       = numel(grp_levels);
cmap       = i_colormap(nGrp, palette);

use_facet = ~isempty(facet_by);
if use_facet
    facet_levels = unique(T.(facet_by), 'stable');
else
    facet_levels = {''};
end
nFacet = numel(facet_levels);
nCols  = nMk;
nRows  = nFacet;

% --- Pre-compute shared x-range per marker column across all facets -------
% Each column (marker) gets the same x-axis limits regardless of facet row,
% so distributions are directly comparable across glucose conditions.
col_xlim = zeros(nMk, 2);
for mi = 1:nMk
    all_vals = T.(markers{mi});
    all_vals = all_vals(~isnan(all_vals));
    if use_log2, all_vals(~isnan(all_vals)) = log2(max(all_vals(~isnan(all_vals)), 1)); end
    if isempty(all_vals)
        col_xlim(mi,:) = [0 1];
    else
        pad = 0.05 * range(all_vals);
        col_xlim(mi,:) = [min(all_vals)-pad, max(all_vals)+pad];
    end
end

if isempty(fig_pos)
    fig_pos = [40 40 min(300*nCols, 2400) 220*nRows + 60];
end

fig = figure('Color','w', 'Position', fig_pos, ...
             'Name', sprintf('Distributions | group=%s facet=%s', group_by, facet_by));

for fi = 1:nFacet
    if use_facet
        Tf          = T(T.(facet_by) == facet_levels(fi), :);
        facet_title = string(facet_levels(fi));
    else
        Tf          = T;
        facet_title = '';
    end

    for mi = 1:nMk
        ax = subplot(nRows, nCols, (fi-1)*nCols + mi);
        hold(ax, 'on');

        all_vals = Tf.(markers{mi});
        all_vals = all_vals(~isnan(all_vals));
        if use_log2, all_vals(~isnan(all_vals)) = log2(max(all_vals(~isnan(all_vals)), 1)); end
        if isempty(all_vals)
            title(ax, markers{mi}, 'Interpreter','none', 'FontSize', 8);
            xlim(ax, col_xlim(mi,:));
            hold(ax,'off'); continue;
        end

        % Use shared x-range for bin edges so all facets bin identically
        edges = linspace(col_xlim(mi,1), col_xlim(mi,2), n_bins + 1);

        for gi = 1:nGrp
            v = Tf.(markers{mi})(Tf.(group_by) == grp_levels(gi));
            v = v(~isnan(v));
            if use_log2, v(~isnan(v)) = log2(max(v(~isnan(v)), 1)); end
            if numel(v) < 3, continue; end

            if show_hist
                histogram(ax, v, edges, 'Normalization', norm, ...
                          'FaceColor', cmap(gi,:), 'FaceAlpha', 0.05, ...
                          'EdgeColor', 'none');
            end

            bw = 1.06 * std(v) * numel(v)^(-1/5);
            if bw == 0, bw = 0.1; end
            [f, xi] = ksdensity(v, 'Bandwidth', bw);

            if strcmp(norm, 'count')
                bin_w = edges(2) - edges(1);
                f = f * numel(v) * bin_w;
            elseif strcmp(norm, 'probability')
                bin_w = edges(2) - edges(1);
                f = f * bin_w;
            end

            fill(ax, [xi, fliplr(xi)], [f, zeros(size(f))], ...
                 cmap(gi,:), 'FaceAlpha', alpha, 'EdgeColor', 'none');
            plot(ax, xi, f, 'Color', cmap(gi,:), 'LineWidth', lw);
        end

        % Apply shared x-limits for this marker column
        xlim(ax, col_xlim(mi,:));

        title(ax, markers{mi}, 'Interpreter','none', 'FontSize', 9);
        if mi == 1
            if use_facet
                ylabel(ax, [char(facet_title) newline norm], 'FontSize', 8, 'Interpreter','none');
            else
                ylabel(ax, norm, 'FontSize', 8);
            end
        else
            ylabel(ax, '');
        end
        xlabel(ax, xlabel_str, 'FontSize', 7, 'Interpreter','none');
        ax.FontSize = 8;
        grid(ax, 'on'); box(ax, 'off');
        hold(ax, 'off');
    end
end

% --- Shared legend on dedicated invisible axes ----------------------------
lh = gobjects(nGrp, 1);
for gi = 1:nGrp
    lh(gi) = patch(NaN, NaN, cmap(gi,:), 'FaceAlpha', max(alpha, 0.4), ...
                   'EdgeColor', cmap(gi,:), 'LineWidth', lw, ...
                   'DisplayName', char(grp_levels(gi)));
end
ax_leg = axes(fig, 'Position',[0.92 0 0.08 1], 'Visible','off');
leg = legend(ax_leg, lh, 'Interpreter','none', 'FontSize', 8);
leg.Box      = 'off';
leg.Units    = 'normalized';
leg.Position = [0.92, 0.5 - leg.Position(4)/2, leg.Position(3), leg.Position(4)];

sgtitle(fig, sprintf('CyCIF distributions  |  group: %s  |  facet: %s', ...
        group_by, facet_by), 'FontWeight','bold', 'FontSize', 11, 'Interpreter','none');

% --- Save as SVG ----------------------------------------------------------
if ~isempty(p.Results.save_path)
    fpath = p.Results.save_path;
    if ~strcmpi(fpath(end-3:end), '.svg'), fpath = [fpath '.svg']; end

    % Fix paper size to match the screen figure exactly — prevents SVG crop
    fig.Units            = 'inches';
    fig_sz               = fig.Position(3:4);   % [width, height] in inches
    fig.PaperUnits       = 'inches';
    fig.PaperSize        = fig_sz;
    fig.PaperPosition    = [0, 0, fig_sz];
    fig.PaperPositionMode = 'manual';

    print(fig, fpath, '-dsvg', '-painters');
    fprintf('Saved SVG: %s\n', fpath);
end
end


% =========================================================================
function markers = i_auto_markers(T)
    cols   = T.Properties.VariableNames;
    is_num = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
    is_id  = ismember(cols, {'xy','cellid','ptx_is_spike'});
    markers = cols(is_num & ~is_id);
end

function cmap = i_colormap(n, palette)
% [n x 3] colormap for categorical groups.
% Palettes: 'bright' | 'parula' | 'lines' | any MATLAB colormap name.
if nargin < 2 || isempty(palette), palette = 'bright'; end

switch lower(palette)
    case 'bright'
        base = [ ...
            0.259, 0.263, 0.282;  % 1  charcoal
            1.000, 0.435, 0.761;  % 2  hot pink
            0.647, 0.000, 0.380;  % 3  deep magenta
            0.463, 0.337, 0.875;  % 4  medium purple
            0.059, 0.420, 0.996;  % 5  vivid blue
            0.776, 0.855, 0.992;  % 6  soft lavender-blue
            0.898, 0.702, 0.847;  % 7  soft pink
            0.875, 0.918, 0.996;  % 8  soft periwinkle
            0.647, 0.647, 0.647;  % 9  mid grey
            0.106, 0.106, 0.106]; % 10 near-black
        if n <= size(base,1)
            cmap = base(1:n,:);
        else
            cmap = lines(n);
            warning('i_colormap: n=%d exceeds palette size, using lines().', n);
        end

    case 'pal1'
        base = [ ...
            0.259, 0.263, 0.282;  % 1 charcoal
            1.000, 0.435, 0.761;  % 2 hot pink
            0.647, 0.000, 0.380;  % 3 deep magenta
            0.463, 0.337, 0.875;  % 4 medium purple
            0.106, 0.106, 0.106;  % 5 near-black
            0.776, 0.855, 0.992]; % 6 soft lavender-blue
        if n <= size(base,1)
            cmap = base(1:n,:);
        else
            cmap = lines(n);
            warning('i_colormap: n=%d exceeds palette size, using lines().', n);
        end

    case 'pal2'
        base = [ ...
            0.945, 0.431, 0.745;  % 1 hot pink
            0.898, 0.702, 0.847;  % 2 soft pink
            0.059, 0.420, 0.996;  % 3 vivid blue
            0.875, 0.918, 0.996;  % 4 soft periwinkle
            0.710, 0.620, 0.800;  % 5 muted lavender
            0.118, 0.122, 0.098;  % 6 near-black
            0.647, 0.647, 0.647]; % 7 mid grey
        if n <= size(base,1)
            cmap = base(1:n,:);
        else
            cmap = lines(n);
            warning('i_colormap: n=%d exceeds palette size, using lines().', n);
        end

    case 'parula'
        raw = parula(max(n,2));
        idx = round(linspace(1, size(raw,1), n));
        cmap = raw(idx,:);

    case 'lines'
        cmap = lines(n);

    otherwise
        try
            raw  = feval(lower(palette), max(n,2));
            idx  = round(linspace(1, size(raw,1), n));
            cmap = raw(idx,:);
        catch
            warning('i_colormap: unknown palette ''%s'', using lines().', palette);
            cmap = lines(n);
        end
end
end