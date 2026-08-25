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
        fprintf('plot_sorted_stacks: facet %d/%d (%s)...\n', ...
            s, numel(prep.facetby), char(prep.facetby{s}));

        % --- Pre-count the tiles this facet will draw, so the FIGURE can be
        % sized to a fixed per-tile dimension (tiles stay the same size
        % whether you ask for 10, 50, or 100 cells; the figure grows). ---
        nTiles = 0;
        contribXY = [];      % which wells actually contribute tiles
        contribN  = [];      % how many tiles each contributes
        for g = 1:numel(prep.groupby)
            xys = intersect(prep.idx.(prep.facetby{s}), prep.idx.(prep.groupby{g}));
            for xy = xys(:)'
                if xy > numel(dl.d) || isempty(dl.d{xy}) || ~isfield(dl.d{xy}.data, chans{1}); continue; end
                if ~prep.goodxy(xy); continue; end
                nAvail = size(dl.d{xy}.data.(chans{1}),1);
                nWantP = p.ntracks * p.nstacks;
                if ~isempty(p.ncells); nWantP = min(nWantP, p.ncells); end
                nWantP = min(nWantP, nAvail);
                % per-replicate breakdown (mirrors the stratified draw below)
                if isfield(dl.d{xy},'source') && ~isempty(dl.d{xy}.source)
                    srcXY0 = dl.d{xy}.source(:,1);
                else
                    srcXY0 = repmat(xy, nAvail, 1);
                end
                reps0  = unique(srcXY0,'stable');
                perRep0 = floor(nWantP / numel(reps0));
                for r = 1:numel(reps0)
                    contribXY(end+1) = reps0(r);                                  %#ok<AGROW>
                    contribN(end+1)  = min(perRep0, sum(srcXY0==reps0(r)));       %#ok<AGROW>
                end
                nTiles = nTiles + nWantP;
            end
        end
        if nTiles == 0; continue; end

        nCols = ceil(sqrt(nTiles));            % near-square grid
        nRows = ceil(nTiles / nCols);
        tileW = 2.2;  tileH = 1.6;             % inches per tile (fixed size)
        figW  = min(nCols*tileW, 60);          % cap so it never exceeds a sane max
        figH_in = min(nRows*tileH + 0.6, 40);  % +0.6 for the layout title

        figH = R.newFigure(sprintf('plot_sorted_stacks | %s | %s', strjoin(chans,'+'), char(prep.facetby{s})));
        set(figH, 'Units','inches', 'Position',[0.5 0.5 figW figH_in]);
        figs(end+1) = figH; %#ok<AGROW>
        T = tiledlayout(figH, nRows, nCols, 'TileSpacing','compact','Padding','compact');

        for g = 1:numel(prep.groupby)
            xys = intersect(prep.idx.(prep.facetby{s}), prep.idx.(prep.groupby{g}));
            for xy = xys(:)'
                if xy > numel(dl.d) || isempty(dl.d{xy}) || ~isfield(dl.d{xy}.data, chans{1}); continue; end
                if ~prep.goodxy(xy); continue; end

                tw = R.timeWindow(dl, xy, chans{1}, prep.linetp{cellfun(@(v)any(v==xy),prep.xymat)}, p);

                % --- Replicate-stratified sampling -----------------------
                % combinexys pools technical replicates into one pseudo-xy,
                % so a single global NaN-ranking would draw most/all tiles
                % from whichever well happened to segment best. Instead:
                % split ntracks equally across the contributing wells, rank
                % by NaN-completeness WITHIN each, and keep replicate blocks
                % contiguous so boundaries are visible.
                Dfull  = dl.d{xy}.data.(chans{1})(:, tw.firsttp:tw.tracklength);
                nanCnt = sum(isnan(Dfull), 2);

                if isfield(dl.d{xy},'source') && ~isempty(dl.d{xy}.source)
                    srcXY = dl.d{xy}.source(:,1);          % original well per cell
                else
                    srcXY = repmat(xy, size(Dfull,1), 1);  % not pooled
                end
                reps   = unique(srcXY, 'stable');
                nReps  = numel(reps);
                nWant  = p.ntracks * p.nstacks;
                if ~isempty(p.ncells); nWant = min(nWant, p.ncells); end
                perRep = floor(nWant / nReps);             % equal split

                I = [];
                for r = 1:nReps
                    idxR   = find(srcXY == reps(r));
                    [~,ord] = sort(nanCnt(idxR));          % rank within replicate
                    take   = min(perRep, numel(idxR));
                    I      = [I; idxR(ord(1:take))];       %#ok<AGROW> blocks stay contiguous
                end
                % top up from the best remaining tracks if equal split left a
                % shortfall (e.g. a replicate had fewer cells than perRep)
                if numel(I) < nWant
                    rest = setdiff((1:size(Dfull,1))', I, 'stable');
                    [~,ord] = sort(nanCnt(rest));
                    I = [I; rest(ord(1:min(nWant-numel(I), numel(rest))))];
                end
                nWant = numel(I);
                % ---------------------------------------------------------

                for k = 1:nWant
                    AH = nexttile(T); hold(AH,'on');
                    xf = tw.firsttp:tw.tracklength;   % absolute frame numbers (x-data)
                    for c = 1:nCh
                        if nCh > 1; if c==1; yyaxis(AH,'left'); else; yyaxis(AH,'right'); end; end
                        D = dl.d{xy}.data.(chans{c})(:, tw.firsttp:tw.tracklength);
                        if ~isempty(p.smooth); D = movmean(D, p.smooth, 2); end
                        rules = R.formatRules(chans{c});
                        plot(AH, xf, D(I(k),:), 'Color', rules.CmapColor, 'LineStyle','-');
                        set(AH,'YColor','k');
                        if iscell(p.ymn) && numel(p.ymn)>=c; ylim(AH,[p.ymn{c}, p.ymx{c}]); end
                    end

                    % title = "Cell <origID> (XY n)" -- XY makes each tile
                    % self-identifying, so a condition backed by 3 replicate
                    % wells is visibly distinguishable from one backed by 1.
                    if isfield(dl.d{xy},'source') && ~isempty(dl.d{xy}.source) && I(k) <= size(dl.d{xy}.source,1)
                        cid    = dl.d{xy}.source(I(k),2);   % original cell id
                        origXY = dl.d{xy}.source(I(k),1);   % original well (pre-combine)
                    else
                        cid    = I(k);
                        origXY = xy;
                    end
                    title(AH, sprintf('Cell %d (XY%d)', cid, origXY), 'FontWeight','normal');

                    R.txLines(AH, tw);
                    R.timeAxis(AH, tw, p);
                    xlim(AH, [tw.firsttp, tw.tracklength]);   % exact crop, no blank pad
                    R.styleTile(AH, nCh);      % declutter: thin ticks, drop repeated y-labels
                    hold(AH,'off');
                end
            end
        end

        % figure title carries the constant context (cell type, channels) so
        % tiles don't have to repeat it; color-codes the channels since the
        % per-tile y-labels were removed for decluttering.
        if nCh > 1
            chanTag = sprintf('%s (L, blue) vs %s (R, red)', chans{1}, chans{2});
        else
            chanTag = chans{1};
        end
        title(T, sprintf('%s  |  %s  |  sorted stacks  |  %s', ...
            strrep(dl.file.base,'_',' '), strrep(chanTag,'_',' '), ...
            strrep(char(prep.facetby{s}),'_',' ')), 'Interpreter','none');

        % Well composition: ntracks is PER WELL, so a condition backed by 3
        % replicate wells contributes 3x the tiles. State it explicitly so
        % conditions with different replicate counts aren't silently compared.
        wellTxt = strjoin(arrayfun(@(a,b) sprintf('XY%d(%d)', a, b), ...
            contribXY, contribN, 'Un', 0), ', ');
        subtitle(T, sprintf('%d replicate well(s): %s  |  %d traces total  |  ntracks=%d split equally, NaN-ranked within each well', ...
            numel(contribXY), wellTxt, nTiles, p.ntracks), ...
            'Interpreter','none', 'FontSize', max(p.font_size-1,6));

        fontsize(figH, p.font_size, 'points'); fontname(figH,'Calibri');
        % V5: harmonize yyaxis limits across cell tiles (per-channel)
        dl_plotstandardize(figH, 'sorted_stacks', chans, p);
        drawnow;   % flush this facet figure to screen so the user sees it build
        % V5: sorted-stacks is a block of many single-cell tiles per facet;
        % saved as ONE file per facet (per-cell splitting would be excessive).
        % Name: base_sortedstacks_channel_facet (no timestamp; re-runs overwrite).
        % V5: encode the manual time window (minutes) in the filename so
        % different window "forks" of the same condition save to distinct files.
        winTag = '';
        if ~isempty(p.tstart_min) || ~isempty(p.tend_min)
            lo = ''; hi = '';
            if ~isempty(p.tstart_min); lo = num2str(p.tstart_min); end
            if ~isempty(p.tend_min);   hi = num2str(p.tend_min);   end
            winTag = sprintf('_t%s-%smin', lo, hi);
        end
        R.exportSVG(figH, dl, sprintf('sortedstacks_%s_%s%s', ...
            strjoin(chans,'-'), char(prep.facetby{s}), winTag), p);
        if p.closefigs && p.save; close(figH); end   % V5: guard on save (v4 closed unconditionally but assumed saving)
    end
end
end