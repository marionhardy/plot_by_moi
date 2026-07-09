function figs = plot_sorted_stacks(plotby, dataloc, varargin)
% PLOT_SORTED_STACKS  Individual cell traces, NaN-completeness sorted.
%
% figs = plot_sorted_stacks(plotby, dataloc, Name, Value, ...)
%
% Opens ONE figure per facet level; within it, each group contributes a
% block of single-cell tiles (one trace per tile), cells sorted by
% ascending NaN count (most-complete tracks first). Reproduces v4's
% 'sorted stacks v2'. Overlays 1-2 channels on left/right y-axes.
%
% Per the v5 spec, each tile title is ONLY "Cell <origID>" -- the long
% condition string that overflowed tiles in v4 is gone (condition context
% lives in the figure/layout title instead).
%
% INPUT (Name/Value; see dl_plotparse for the shared set)
%   'channel'  : {'HYLIGHT'} or overlay {{'HYLIGHT','RAMPKAR2'}}
%   'ntracks'  : cells per group block (default 5)   [renderer-specific]
%   'nstacks'  : stack multiplier (default 1)         "
%   'ymn','ymx': per-channel y-limits (cell).
%
% OUTPUT
%   figs : array of figure handles (one per facet).
%
% -------------------------------------------------------------------------
% V5 (DatalocHandler_v4 -> v5): clean rewrite of v4's 'sorted stacks v2'.
% Title reduced to "Cell <id>" (original request). Prep -> dl_plotprep;
% plumbing -> dl_plotrender (figure-per-facet fix + SVG); params ->
% dl_plotparse. Channel-kind guarded (crit #5). Original cell ID recovered
% via dl.d{xy}.source (populated by dl_combine_xys) so titles/ct_cell_video
% stay traceable after replicate combining.
% -------------------------------------------------------------------------

extra.ntracks = {5, @(x) isnumeric(x)&&isscalar(x)};
extra.nstacks = {1, @(x) isnumeric(x)&&isscalar(x)};
p = dl_plotparse(plotby, dataloc, extra, varargin{:});

prep = dl_plotprep(p.plotby, dataloc, p);
R    = dl_plotrender();
dl   = prep.dataloc;
figs = gobjects(0);

% p.channel: list of plot-specs; each spec a cellstr of overlaid channels.
for iSpec = 1:numel(p.channel)
    chans = p.channel{iSpec};
    nCh   = numel(chans);

    for c = 1:nCh
        if strcmp(dl_channel_kind(dataloc, chans{c}), 'endpoint')
            error('plot_sorted_stacks:EndpointChannel', ...
                'Channel "%s" is an endpoint ([N x 1]) measurement, not a time series.', chans{c});
        end
    end

    for s = 1:numel(prep.facetby)
        figH = R.newFigure(sprintf('plot_sorted_stacks | %s | %s', strjoin(chans,'+'), char(prep.facetby{s})));
        figs(end+1) = figH; %#ok<AGROW>
        T = tiledlayout(figH,'flow','TileSpacing','compact','Padding','compact');

        for g = 1:numel(prep.groupby)
            xys = intersect(prep.idx.(prep.facetby{s}), prep.idx.(prep.groupby{g}));
            for xy = xys(:)'
                if xy > numel(dl.d) || isempty(dl.d{xy}) || ~isfield(dl.d{xy}.data, chans{1}); continue; end
                if ~prep.goodxy(xy); continue; end

                tw = R.timeWindow(dl, xy, chans{1}, prep.linetp{cellfun(@(v)any(v==xy),prep.xymat)}, p);

                % how many cells this block: cap to availability
                nCells = size(dl.d{xy}.data.(chans{1}),1);
                if ~isempty(p.ncells) && p.ncells <= nCells; nCells = p.ncells; end
                nWant = p.ntracks * p.nstacks;
                if nCells < nWant; nWant = floor(nCells/p.nstacks)*p.nstacks; end

                % sort by NaN-completeness over the plotted window (v4 logic)
                [~,I] = sort(sum(isnan(dl.d{xy}.data.(chans{1})(:, tw.firsttp:tw.tracklength)),2));

                for k = 1:nWant
                    AH = nexttile(T); hold(AH,'on');
                    for c = 1:nCh
                        if nCh > 1; if c==1; yyaxis(AH,'left'); else; yyaxis(AH,'right'); end; end
                        D = dl.d{xy}.data.(chans{c})(:, tw.firsttp:tw.tracklength);
                        if ~isempty(p.smooth); D = movmean(D, p.smooth, 2); end
                        rules = R.formatRules(chans{c});
                        plot(AH, D(I(k),:), 'Color', rules.CmapColor, 'LineStyle','-');
                        set(AH,'YColor','k');
                        ylabel(AH, rules.Ylabel, 'Color', rules.CmapColor);
                        if iscell(p.ymn) && numel(p.ymn)>=c; ylim(AH,[p.ymn{c}, p.ymx{c}]); end
                    end

                    % title = "Cell <origID>" only (recover ID through .source)
                    if isfield(dl.d{xy},'source') && ~isempty(dl.d{xy}.source) && I(k) <= size(dl.d{xy}.source,1)
                        cid = dl.d{xy}.source(I(k),2);
                    else
                        cid = I(k);
                    end
                    title(AH, sprintf('Cell %d', cid), 'FontWeight','normal');

                    R.txLines(AH, tw);
                    R.timeAxis(AH, tw, p);
                    xlim(AH, [tw.firsttp-1, tw.tracklength]);
                    R.styleAxis(AH);
                    hold(AH,'off');
                end
            end
        end

        title(T, sprintf('%s | %s | sorted stacks | %s', strrep(dl.file.base,'_',' '), ...
            strjoin(strrep(chans,'_',' '),', '), strrep(char(prep.facetby{s}),'_',' ')), ...
            'Interpreter','none');
        fontsize(figH, p.font_size, 'points'); fontname(figH,'Calibri');
        % V5: harmonize yyaxis limits across cell tiles (per-channel)
        dl_plotstandardize(figH, 'sorted_stacks', chans, p);
        % V5: sorted-stacks is a block of many single-cell tiles per facet;
        % saved as ONE file per facet (per-cell splitting would be excessive).
        % Name: base_sortedstacks_channel_facet (no timestamp; re-runs overwrite).
        R.exportSVG(figH, dl, sprintf('sortedstacks_%s_%s', ...
            strjoin(chans,'-'), char(prep.facetby{s})), p);
        if p.closefigs && p.save; close(figH); end   % V5: guard on save (v4 closed unconditionally but assumed saving)
    end
end
end
