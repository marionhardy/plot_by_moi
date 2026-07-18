function R = dl_plotrender()
% DL_PLOTRENDER  Shared rendering helpers for the v5 plot suite.
%
% R = dl_plotrender() returns a struct of function handles so the three
% renderers share ONE implementation of the cross-cutting plumbing (figure
% creation, time-axis ticks, tx-lines, channel colors, per-xy time window,
% SVG export). Call e.g. R = dl_plotrender(); R.newFigure(name).
%
% HANDLES
%   R.newFigure(name)                       -> figure handle (one per facet)
%   R.formatRules(channel)                  -> struct array, .Ylabel/.CmapColor
%   R.timeWindow(dataloc,xy,chan1,linetp,p) -> struct .firsttp/.tracklength/
%                                              .thisTX/.txs
%   R.timeAxis(AH, tw, p)                    -> set x-ticks in hours (findZero)
%   R.txLines(AH, tw)                        -> draw treatment xlines
%   R.exportSVG(figH, dataloc, tag, p)       -> save .svg to fig dir
%   R.safeTitle(txt, maxChars)               -> truncated, interpreter-safe
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): v4 scattered this plumbing through
% main_plotting_func, subplot_standardizer, and a local findZero. The first
% v5 port duplicated findZero/format-rules into each plot_ file. Here it is
% ONE shared module. newFigure fixes the figure-per-facet regression by
% design (every facet gets its own figure explicitly). exportSVG restores
% v4's printstyle='svg' capability that the port had dropped. timeAxis is
% v4's findZero verbatim in behaviour, just parameterized on a tw struct.
% -------------------------------------------------------------------------

R.newFigure   = @newFigure;
R.formatRules = @formatRules;
R.timeWindow  = @timeWindow;
R.timeAxis    = @timeAxis;
R.txLines     = @txLines;
R.exportSVG   = @exportSVG;
R.safeTitle   = @safeTitle;
R.styleAxis   = @styleAxis;
R.styleTile   = @styleTile;
end


function styleTile(AH, nCh)
% Decluttered styling for dense tile grids (sorted stacks). Keeps per-tile
% axes (light structural change) but removes the visual noise that made the
% 50-tile figure unreadable: repeated channel y-labels dropped, ticks thin
% and outward, tick LABELS thinned so numbers don't collide at small size.
% Traces keep the full pixel budget. nCh = 1 or 2 (yyaxis).
box(AH, 'on');
set(AH, 'Layer','top', 'TickDir','out', 'LineWidth',0.4, ...
        'XMinorTick','off', 'YMinorTick','off', 'TickLength',[0.03 0.03]);

if nCh > 1
    for side = {'left','right'}
        yyaxis(AH, side{1});
        ylabel(AH, '');
        thinYTicks(AH);
    end
    yyaxis(AH, 'left');           % leave active side deterministic
else
    ylabel(AH, '');
    thinYTicks(AH);
end
end


function thinYTicks(AH)
% Keep only min / 0 / max y-ticks so the axis reads without a dense stack.
yl = ylim(AH);
if yl(1) < 0 && yl(2) > 0
    yticks(AH, unique([yl(1) 0 yl(2)]));
else
    yticks(AH, unique([yl(1) yl(2)]));
end
end


function styleAxis(AH)
% Clean, modern tile framing. yyaxis can leave the top spine missing; box on
% closes it. Thin axis lines, ticks pointing out, no minor ticks.
box(AH, 'on');
set(AH, 'Layer', 'top', 'TickDir', 'out', 'LineWidth', 0.5, 'XMinorTick','off','YMinorTick','off');
end


function figH = newFigure(name)
% One figure per facet level (fixes the v4->v5 figure-per-facet regression).
figH = figure('Name', name, 'NumberTitle', 'off', 'Color', 'w');
end


function tw = timeWindow(dataloc, xy, chan1, linetp_xy, p)
% Resolve the plotted time window for one xy: treatment tp, first/last tp.
% Mirrors v4's firsttp/tracklength/thisTX logic exactly.
% linetp_xy may be a numeric vector (v5 prep: ALL sub-treatment times) or a
% legacy char/cell string; coerce to a numeric row of positive timepoints.
if iscell(linetp_xy); linetp_xy = [linetp_xy{:}]; end
if ischar(linetp_xy) || isstring(linetp_xy)
    treat = str2double(linetp_xy);
else
    treat = double(linetp_xy);
end
u = unique(treat(:)'); u = u(~isnan(u)); u = u(u>0);
txs = u;

thisTX = [];
if ~isempty(txs)
    if numel(txs) < p.aftertreatment; thisTX = txs(end);
    else; thisTX = txs(p.aftertreatment); end
end

tMax = size(dataloc.d{xy}.data.(chan1), 2);
if isempty(p.tmaxaftertx); tracklength = tMax; else; tracklength = p.tmaxaftertx + thisTX; end
if tracklength > tMax; tracklength = tMax; end

firsttp = 1;
if p.plotfromzero
    if ~isempty(txs); firsttp = thisTX; end
elseif ~isempty(p.tbeforetx)
    firsttp = thisTX - p.tbeforetx;
    if firsttp > 0
        if ~isempty(txs); txs = txs - firsttp; end
    else
        firsttp = 1;
    end
end
if ~isempty(txs); txs = txs(txs < tracklength); end

tw = struct('firsttp',firsttp, 'tracklength',tracklength, 'thisTX',thisTX, 'txs',txs);
end


function txLines(AH, tw)
% Fine dotted vertical lines at every treatment timepoint.
if ~isempty(tw.txs)
    xline(AH, tw.txs, ':', 'LineWidth', 0.5, 'Color', [0.35 0.35 0.35], 'Alpha', 0.9);
end
end


function timeAxis(AH, tw, p)
% Time axis in HOURS, zeroed on the treatment. Chooses a label step that
% yields ~6-8 ticks on clean hour boundaries so labels never collide.
% (v4 used a fixed 2h step which turned long movies into a black block.)
spanHrs = (tw.tracklength - tw.firsttp) / p.tktm;

% pick an hour step from a "nice" ladder targeting <= ~8 ticks
ladder = [1 2 5 10 15 20 25 50];
stepHr = ladder(find(spanHrs./ladder <= 8, 1));
if isempty(stepHr); stepHr = ceil(spanHrs/8); end
step = p.tktm * stepHr;                       % back to frames

if isempty(p.zerohrtx); zeroTx = tw.thisTX;
else
    if numel(tw.txs) < p.zerohrtx; zeroTx = tw.txs(end); else; zeroTx = tw.txs(p.zerohrtx); end
end

if ~isempty(zeroTx) && ~p.plotfromzero
    back = zeroTx:-step:tw.firsttp; back = sort(back(2:end));
    xt   = [back, zeroTx:step:tw.tracklength];
    xl   = round((xt - zeroTx)/p.tktm);
elseif ~isempty(zeroTx) && p.plotfromzero
    xt = tw.firsttp:step:tw.tracklength;
    xl = round((xt - tw.firsttp)/p.tktm);
else
    xt = tw.firsttp:step:tw.tracklength;
    xl = round(xt/p.tktm);
end
xticks(AH, xt); xticklabels(AH, string(xl)); xtickangle(AH, 0);
AH.XAxis.TickLength = [0.02 0.025];            % modest tick marks
xlabel(AH, '');
end


function PlotRules = formatRules(channels, plottype)
% Channel -> (Ylabel, CmapColor). Trimmed from v4 PlotFormatRules.
% V5 FIX: 'plottype' arg restores v4's heatmap override (L1285) that the
% first v5 rewrite dropped -- heatmaps force the parula COLORMAP instead of
% a per-sensor hex color. Passing a hex string as a colormap was the
% flat-blue-heatmap root cause.
if nargin < 2; plottype = ''; end
List = {...
 'RAMPKAR',  'AMPK Activity (a.u.)',      '#C30F0E', '^(RAMPKAR2)$'; ...
 'Hylight',  'FBP level (a.u.)',          '#1976D2', '(n)?(HYLIGHT)'; ...
 'Perceval', 'ATP:ADP ratio',             '#168039', '(n)?(PERCEVAL)'; ...
 'AMPK',     'AMPK Activity (a.u.)',      '#04BFBF', '^(AMPKAR)$'; ...
 'ERK',      'ERK Activity (a.u.)',       '#C30F0E', '(n)?E(R)?K(TR|AR)'; ...
 'JNK',      'JNK Activity (a.u.)',       '#D17300', '(n)?JNKTR'; ...
 'AKT',      'AKT Activity (a.u.)',       '#168039', '(n)?AKT(2)?KTR'; ...
 'P65KTR',   'NFkB Activity (a.u.)',      '#168039', '(n)?P65(KTR)?' ...
 };
if ischar(channels); channels = {channels}; end
hits = cellfun(@(x)find(~cellfun(@isempty, regexpi(x, List(:,4))),1), channels, 'Un', 0);
for ii = 1:numel(channels)
    if isempty(hits{ii})
        PlotRules(ii).Ylabel = [channels{ii} ' (a.u.)']; %#ok<AGROW>
        PlotRules(ii).CmapColor = 'red';
    else
        PlotRules(ii).Ylabel    = List{hits{ii},2}; %#ok<AGROW>
        PlotRules(ii).CmapColor = List{hits{ii},3};
    end
    % V5 FIX (root cause of flat-blue heatmap): heatmaps use the parula
    % COLORMAP, not the per-sensor hex line color (v4 PlotFormatRules L1285).
    if contains(plottype, 'heatmap'); PlotRules(ii).CmapColor = 'parula'; end
end
end


function txt = safeTitle(txt, maxChars)
% Wrap dense condition strings onto 2 lines at the first '+' (e.g. splits
% "Glucose17mM+IL6+LPS+10Apuro1X" -> "Glucose17mM" / "IL6+LPS+10Apuro1X"),
% then ellipsis-truncate the tail if still too long. Prevents the single
% dense line that stretched across each tile.
if nargin < 2; maxChars = 22; end
if isstring(txt); txt = char(txt); end
plus = strfind(txt, '+');
if ~isempty(plus)
    head = txt(1:plus(1)-1);
    tail = txt(plus(1)+1:end);
    if numel(tail) > maxChars; tail = [tail(1:maxChars-1) char(8230)]; end
    txt = {head, tail};                       % 2-line title (cellstr)
elseif numel(txt) > maxChars
    txt = [txt(1:maxChars-1) char(8230)];
end
end


function exportSVG(figH, dataloc, tag, p)
% Save the WHOLE figure as one vector SVG per facet.
% Uses print(...,'-dsvg','-painters') rather than exportgraphics/saveas:
% print operates on the FIGURE (not a tile axes), so it sidesteps the
% "multiple coordinate systems" / container errors that exportgraphics
% raises on a yyaxis tile. '-painters' forces true vector output (not a
% flattened raster), which is required for Illustrator/Inkscape editing.
if ~p.save; return; end
if isfield(dataloc,'fold') && isfield(dataloc.fold,'fig') && ~isempty(dataloc.fold.fig)
    figdir = dataloc.fold.fig;
else
    figdir = pwd;
    warning('dl_plotrender:NoFigDir','No dataloc.fold.fig; saving SVG to pwd.');
end
if ~exist(figdir,'dir'); mkdir(figdir); end

% V5 NAMING: use the experiment DATE (leading yyyy-mm-dd of the folder/exp
% name) as the prefix, not the full base string. The date is the stable,
% human-meaningful identifier; the long base produced unwieldy names, and
% makeValidName mangled the dashes (2026-06-25 -> x2026_06_25). We keep the
% dash-form date, sanitize only the tag, and join with '__'.
dateStr = local_expdate(dataloc);
tagStr  = matlab.lang.makeValidName(tag);      % tag already carries facet/chan
fname   = sprintf('%s__%s', dateStr, tagStr);

set(figH, 'Renderer', 'painters');            % force vector, not opengl raster
fontname(figH, 'Calibri');
try
    print(figH, fullfile(figdir, [fname '.svg']), '-dsvg', '-painters');
catch ME
    warning('dl_plotrender:SVGExportFailed', ...
        'SVG export failed for %s (%s); skipping.', fname, ME.message);
end
end


function d = local_expdate(dataloc)
% Extract the experiment date (yyyy-mm-dd) from the start of the folder/exp
% name. Falls back through exp.name -> file.base -> 'nodate'.
src = '';
if isfield(dataloc,'exp') && isfield(dataloc.exp,'name') && ~isempty(dataloc.exp.name)
    src = dataloc.exp.name;
elseif isfield(dataloc,'file') && isfield(dataloc.file,'base') && ~isempty(dataloc.file.base)
    src = dataloc.file.base;
end
tok = regexp(char(src), '^(\d{4}-\d{2}-\d{2})', 'tokens', 'once');
if ~isempty(tok); d = tok{1}; else; d = 'nodate'; end
end