function dl_plotstandardize(figH, kind, channel, p)
% DL_PLOTSTANDARDIZE  Harmonize axis/color limits across tiles of ONE figure.
%
% dl_plotstandardize(figH, kind, channel, p)
%
% Post-draw pass: after a renderer has drawn every tile of a figure, this
% walks that figure's axes and sets a common CLim (heatmap) or y-limit
% (mean / sorted stacks) so tiles are visually comparable. Scope is
% per-figure (one facet), matching v4.
%
% INPUT (paired at the renderer call site)
%   figH    : figure handle to standardize.
%   kind    : 'heatmap' | 'mean' | 'sorted_stacks'
%   channel : cellstr of channels in this figure (1 = single, 2 = yyaxis).
%   p       : param struct; uses p.ymn, p.ymx (cells), p.standardizeplots.
%
% OUTPUT
%   none (mutates axes of figH in place).
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): extracts v4's ~200-line inline
% subplot_standardizer into a discrete shared function so the three
% renderers stay flat. Limit math is faithful to v4 (plot_by_MHY.m
% L1104-1204): heatmap -> common CLim from p.ymn/p.ymx or data range;
% sorted_stacks/mean -> common yyaxis limits per channel. NANCOLOR is NOT
% applied this pass (deferred by decision); MATLAB default NaN rendering.
% CLim for heatmaps now lives HERE, not in a per-tile clim() call inside
% the renderer (the per-tile approach was not the v4 mechanism and
% interacted badly with the colormap).
% -------------------------------------------------------------------------

if ~p.standardizeplots; return; end
if isempty(figH) || ~isgraphics(figH); return; end

ax = findall(figH, 'Type', 'Axes');
if isempty(ax); return; end

switch kind
    case 'heatmap'
        % Common CLim across all heatmap tiles.
        if ~isempty(p.ymn) && ~isempty(p.ymx)
            cvals = [p.ymn{1}, p.ymx{1}];
        else
            cmins = arrayfun(@(a) a.CLim(1), ax);
            cmaxs = arrayfun(@(a) a.CLim(2), ax);
            % v4 nudges the auto range slightly inward
            cvals = [min(cmins)*1.05, max(cmaxs)*0.95];
        end
        set(ax, 'CLim', cvals);

    case {'mean','sorted_stacks'}
        % Common y-limits; per-channel when 2-channel yyaxis overlay.
        Y = get(ax, 'YAxis');
        if iscell(Y); Y = [Y{:}]; end   % ax x nAxes rulers
        if numel(channel) > 1
            % Y is [2 x nAxes]: row1 = left, row2 = right
            if ~isempty(p.ymn) && ~isempty(p.ymx)
                lLim = [p.ymn{1}, p.ymx{1}];
                rLim = [p.ymn{2}, p.ymx{2}];
            else
                lRange = [Y(1,:).Limits]; rRange = [Y(2,:).Limits];
                lLim = [min(lRange(:)), max(lRange(:))];
                rLim = [min(rRange(:)), max(rRange(:))];
            end
            % V5: align zero-lines so red/blue crossings are meaningful.
            % Rescale the RIGHT axis to match the LEFT axis's zero fraction
            % (fraction of the range below zero), preserving its span.
            if p.alignzero
                fLeft = -lLim(1) / (lLim(2) - lLim(1));   % 0..1 below-zero frac
                rSpan = rLim(2) - rLim(1);
                rLim  = [-fLeft*rSpan, (1-fLeft)*rSpan];
            end
            set(Y(1,:), 'Limits', lLim);
            set(Y(2,:), 'Limits', rLim);
        else
            if ~isempty(p.ymn) && ~isempty(p.ymx)
                set(Y(:), 'Limits', [p.ymn{1}, p.ymx{1}]);
            else
                Range = [Y(:).Limits];
                set(Y(:), 'Limits', [min(Range(:)), max(Range(:))]);
            end
        end
end
end
