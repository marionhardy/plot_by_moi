function figs = plot_mean(plotby, dataloc, varargin)
% PLOT_MEAN  Per-condition mean traces (with quantile bands).
%
% figs = plot_mean(plotby, dataloc, Name, Value, ...)
%
% Opens ONE figure per facet level (cell type when plotby='treatment',
% treatment when plotby='celltype'), tiling one mean trace per group.
% Overlays 1-2 channels on left/right y-axes.
%
% INPUT (all Name/Value pairs; superset shared across the v5 suite)
%   plotby      : 'treatment' | 'celltype' | 'custom'
%   dataloc     : v5 dataloc (legacy tolerated; see dl_check_platemap_compat)
%   'channel'   : {'HYLIGHT'} or overlay {{'HYLIGHT','RAMPKAR2'}}
%   'combinexys': pool replicate xys (default true)
%   'ymn','ymx' : per-channel y-limits, cell (e.g. {-0.2,0} / {0.15,0.2})
%   'qtiles'    : draw quantile bands (default true)      [renderer-specific]
%   'overlapmeans': overlay all groups' means in one tile (default true)  "
%   'subset','exclude','specsubset','nogene','smooth','ncells',
%   'aftertreatment','tmaxaftertx','tbeforetx','plotfromzero','zerohrtx',
%   'tx_order','font_size','save' : see dl_plotparse.
%
% OUTPUT
%   figs : array of figure handles (one per facet), for caller-side save.
%
% -------------------------------------------------------------------------
% V5 (DatalocHandler_v4 -> v5): clean rewrite. Replaces v4's 'mean' case
% inside a 560-line switch. Prep -> dl_plotprep (crit #4/#6); plumbing ->
% dl_plotrender (figure-per-facet fix + SVG export); params -> dl_plotparse
% ('plottype' gone; the function IS the plottype). Channel-kind guarded via
% dl_channel_kind (crit #5). Byte-for-byte identical ct_trackvis_MHY 'mean'
% call as v4, so trace + quantile-band output matches.
% -------------------------------------------------------------------------

extra.qtiles       = {true,  @islogical};
extra.overlapmeans = {true,  @islogical};
p = dl_plotparse(plotby, dataloc, extra, varargin{:});

prep = dl_plotprep(p.plotby, dataloc, p);
R    = dl_plotrender();
dl   = prep.dataloc;
figs = gobjects(0);

% p.channel is a list of plot-specs; each spec is a cellstr of channels to
% overlay (left/right yyaxis). One figure-set per spec, per facet.
for iSpec = 1:numel(p.channel)
    chans = p.channel{iSpec};                 % cellstr, 1-2 channels
    nCh   = numel(chans);

    % endpoint channels ([N x 1] CycIF) can't go on a time-series plot (crit #5)
    for c = 1:nCh
        if strcmp(dl_channel_kind(dataloc, chans{c}), 'endpoint')
            error('plot_mean:EndpointChannel', ...
                'Channel "%s" is an endpoint ([N x 1]) measurement, not a time series.', chans{c});
        end
    end

    for s = 1:numel(prep.facetby)
        figH = R.newFigure(sprintf('plot_mean | %s | %s', strjoin(chans,'+'), char(prep.facetby{s})));
        figs(end+1) = figH; %#ok<AGROW>
        T = tiledlayout(figH,'flow','TileSpacing','compact','Padding','compact');
        tileAx = gobjects(0); tileGrp = {};   % collect tiles for per-condition SVG

        for g = 1:numel(prep.groupby)
            xys = intersect(prep.idx.(prep.facetby{s}), prep.idx.(prep.groupby{g}));
            % V5 FIX (overlapmeans semantics): v4 resets firstxy PER GROUP so
            % each treatment gets its OWN tile; only replicate XYs *within* a
            % group overlay. The prior v5 rewrite kept a single figure-scoped
            % flag, collapsing every treatment into one axis. Restores v4
            % (plot_by_MHY.m ~L547-558).
            firstxy = 0;
            for xy = xys(:)'
                if xy > numel(dl.d) || isempty(dl.d{xy}) || ~isfield(dl.d{xy}.data, chans{1}); continue; end
                if ~prep.goodxy(xy); continue; end

                if p.overlapmeans
                    firstxy = firstxy + 1;
                    if firstxy == 1; AH = nexttile(T); hold(AH,'on'); end
                else
                    AH = nexttile(T); hold(AH,'on');
                end

                tw = R.timeWindow(dl, xy, chans{1}, prep.linetp{cellfun(@(v)any(v==xy),prep.xymat)}, p);

                for c = 1:nCh
                    if nCh > 1; if c==1; yyaxis(AH,'left'); else; yyaxis(AH,'right'); end; end
                    if iscell(p.ymn); ymn = p.ymn{c}; ymx = p.ymx{c}; else; ymn = p.ymn; ymx = p.ymx; end
                    D = dl.d{xy}.data.(chans{c})(:, tw.firsttp:tw.tracklength);
                    if ~isempty(p.smooth); D = movmean(D, p.smooth, 2); end
                    rules = R.formatRules(chans{c});
                    ct_trackvis_MHY(AH, 'mean', D, 'cmap', rules.CmapColor, ...
                        'ymn', ymn, 'ymx', ymx, 'PLOTQTILES', p.qtiles);
                    set(AH,'YColor','k');
                    ylabel(AH, rules.Ylabel, 'Color', rules.CmapColor);
                end

                % decorate axis + title once per tile: overlap -> first xy of
                % the group; non-overlap -> every xy is its own tile.
                if ~p.overlapmeans || firstxy == 1
                    R.txLines(AH, tw);
                    R.timeAxis(AH, tw, p);
                    xlim(AH, [tw.firsttp-1, tw.tracklength]);
                    title(AH, R.safeTitle(prep.legname{g}), 'Interpreter','none');
                    R.styleAxis(AH);
                    tileAx(end+1)  = AH;                 %#ok<AGROW> collect for export
                    tileGrp{end+1} = prep.groupby{g};    %#ok<AGROW>
                end
            end
        end

        title(T, sprintf('%s | %s | mean', strrep(dl.file.base,'_',' '), ...
            strjoin(strrep(chans,'_',' '),', ')), 'Interpreter','none');
        fontsize(figH, p.font_size, 'points'); fontname(figH,'Calibri');
        % V5: harmonize y-limits across tiles (per-channel if yyaxis overlay)
        dl_plotstandardize(figH, 'mean', chans, p);
        % V5: one SVG per condition (tile), named base_plottype_channel_facet_group
        for it = 1:numel(tileAx)
            R.exportTileSVG(tileAx(it), dl, struct('plottype','mean', ...
                'channel',strjoin(chans,'-'), 'facet',char(prep.facetby{s}), ...
                'group',tileGrp{it}), p);
        end
        if p.closefigs && p.save; close(figH); end   % V5: guard on save (v4 closed unconditionally but assumed saving)
    end
end
end
