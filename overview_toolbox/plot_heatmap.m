function figs = plot_heatmap(plotby, dataloc, varargin)
% PLOT_HEATMAP  Single-cell heatmaps, one tile per group, selectable sort.
%
% figs = plot_heatmap(plotby, dataloc, Name, Value, ...)
%
% Opens ONE figure per facet level; each tile is a cells x time heatmap
% for one group. Single channel (heatmaps don't overlay).
%
% INPUT (Name/Value; see dl_plotparse for the shared set)
%   'channel'  : single channel, e.g. {'HYLIGHT'}
%   'sortmode' : 'none' | 'nan' | 'mean' | 'delta'   [renderer-specific]
%                none  = track order            (v4 'heatmap')
%                nan   = ascending NaN count     (v4 'sorted heatmap')
%                mean  = descending mean signal  (v4 'mean sorted heatmap')
%                delta = descending (max-min)    (v4 'delta heatmap', fixed)
%   'ymn','ymx': color-axis limits (clim), e.g. {-0.2}/{0.3}
%   'ncells'   : cap rows per tile.
%
% OUTPUT
%   figs : array of figure handles (one per facet).
%
% -------------------------------------------------------------------------
% V5 (DatalocHandler_v4 -> v5): clean rewrite consolidating v4's FOUR
% heatmap cases into one 'sortmode' param. BUGFIX carried: v4's 'delta
% heatmap' computed a range but never sorted by it (dead code) -> 'delta'
% here actually sorts by (max-min). Prep -> dl_plotprep; plumbing ->
% dl_plotrender (figure-per-facet + SVG); params -> dl_plotparse. Channel-
% kind guarded (crit #5). Identical ct_trackvis_MHY 'heatmap' call as v4.
% -------------------------------------------------------------------------

extra.sortmode = {'nan', @(x) any(validatestring(x,{'none','nan','mean','delta'}))};
p = dl_plotparse(plotby, dataloc, extra, varargin{:});

% p.channel is a list of plot-specs (each a cellstr). Heatmaps don't overlay,
% so take the FIRST channel of each spec; one figure-set per channel.
heatChans = cellfun(@(spec) spec{1}, p.channel, 'Un', 0);

prep = dl_plotprep(p.plotby, dataloc, p);
R    = dl_plotrender();
dl   = prep.dataloc;
figs = gobjects(0);

for iCh = 1:numel(heatChans)
    chan = heatChans{iCh};
    if strcmp(dl_channel_kind(dataloc, chan), 'endpoint')
        error('plot_heatmap:EndpointChannel', ...
            'Channel "%s" is an endpoint ([N x 1]) measurement, not a time series.', chan);
    end

    for s = 1:numel(prep.facetby)
        figH = R.newFigure(sprintf('plot_heatmap | %s | %s', chan, char(prep.facetby{s})));
        figs(end+1) = figH; %#ok<AGROW>
        T = tiledlayout(figH,'flow','TileSpacing','compact','Padding','compact');
        tileAx = gobjects(0); tileGrp = {};   % collect tiles for per-condition SVG

        for g = 1:numel(prep.groupby)
            xys = intersect(prep.idx.(prep.facetby{s}), prep.idx.(prep.groupby{g}));
            for xy = xys(:)'
                if xy > numel(dl.d) || isempty(dl.d{xy}) || ~isfield(dl.d{xy}.data, chan); continue; end
                if ~prep.goodxy(xy); continue; end

                tw = R.timeWindow(dl, xy, chan, prep.linetp{cellfun(@(v)any(v==xy),prep.xymat)}, p);
                D  = dl.d{xy}.data.(chan)(:, tw.firsttp:tw.tracklength);

                switch p.sortmode                                   % row order
                    case 'none';  I = 1:size(D,1);
                    case 'nan';   [~,I] = sort(sum(isnan(D),2));
                    case 'mean';  [~,I] = sort(mean(D,2,'omitnan'),'descend');
                    case 'delta'; [~,I] = sort(max(D,[],2,'omitnan')-min(D,[],2,'omitnan'),'descend');
                end
                D = D(I,:);
                if ~isempty(p.ncells); D = D(1:min(p.ncells,size(D,1)),:); end

                AH = nexttile(T);
                rules = R.formatRules(chan, 'heatmap');   % V5: parula colormap
                ct_trackvis_MHY(AH, 'heatmap', D, 'cmap', rules.CmapColor, 'nolabel', true);
                set(AH,'YDir','reverse'); axis(AH,'tight');
                R.timeAxis(AH, tw, p);
                box(AH,'on'); set(AH,'Layer','top','TickDir','out');  % frame; no styleAxis (keep image tight)
                title(AH, R.safeTitle(prep.legname{g}), 'Interpreter','none');
                tileAx(end+1)  = AH;                 %#ok<AGROW>
                tileGrp{end+1} = prep.groupby{g};    %#ok<AGROW>
            end
        end

        title(T, sprintf('%s | %s | heatmap (%s)', strrep(dl.file.base,'_',' '), ...
            strrep(chan,'_',' '), p.sortmode), 'Interpreter','none');
        fontsize(figH, p.font_size, 'points'); fontname(figH,'Calibri');
        % V5: CLim harmonized across tiles here (not per-tile), then per-tile save
        dl_plotstandardize(figH, 'heatmap', {chan}, p);
        for it = 1:numel(tileAx)
            R.exportTileSVG(tileAx(it), dl, struct('plottype',['heatmap_' p.sortmode], ...
                'channel',chan, 'facet',char(prep.facetby{s}), 'group',tileGrp{it}), p);
        end
        if p.closefigs && p.save; close(figH); end   % V5: guard on save (v4 closed unconditionally but assumed saving)
    end
end
end
