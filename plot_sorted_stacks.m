function plot_sorted_stacks(dataloc, varargin)
% plot_sorted_stacks  Stacked single-cell traces with optional multi-channel overlay.
%
% Each facet level gets one figure. Within the figure, rows are cells; columns
% are groups. For multi-channel calls, channels are overlaid per row using
% yyaxis left/right. Cell matching across channels uses cellindex (inner join:
% only cells present AND non-all-NaN in ALL channels are kept). Rows are
% sorted by NaN count in channel{1} (fewest NaNs first). Row 1 of each column
% is a mean +/- shading panel; rows 2..ntracks+1 are individual cells.
%
% USAGE
%   plot_sorted_stacks(dataloc, ...
%       'channel',     {'HYLIGHT','RAMPKAR2'}, ...
%       'facetby',     {'Gluc'}, ...
%       'groupby',     {'Oligo','Vehicle'}, ...
%       'ntracks',     10, ...
%       'tmaxaftertx', 6, ...
%       'tbeforetx',   1)
%
% INPUTS
%   dataloc       v5 dataloc struct (.d, .platemapd.pmd, .movieinfo)
%
% PARAMETERS
%   channel       cell of channel name strings        (default: {'HYLIGHT'})
%   facetby       cell of pmd token strings           (default: {})
%   groupby       cell of pmd token strings           (default: {})
%   ntracks       int, max cells per group            (default: 10)
%   combinexys    logical, pool XYs within a group    (default: true)
%   tmaxaftertx   double hours, crop after tx         (default: [])
%   tbeforetx     double hours, crop before tx        (default: [])
%   smooth        int, movmean window frames          (default: [])
%   ylims         cell of [mn mx], one per channel    (default: auto)
%   errortype     'std'|'sem'|'iqr'                  (default: 'std')
%   font_size     double pt                           (default: 8)
%   exclude       cell of exclusion token strings     (default: {})
%   savedir       char, output folder                 (default: pwd/figures)
%   closefigs     logical                             (default: true)
%   track_alpha   double 0-1, raw trace opacity       (default: 0.20)

%% ---- Input parsing -------------------------------------------------------
ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'channel',     {'HYLIGHT'}, @(x) ischar(x)||isstring(x)||iscell(x));
addParameter(ip, 'facetby',     {},          @(x) iscell(x)||isstring(x)||ischar(x));
addParameter(ip, 'groupby',     {},          @(x) iscell(x)||isstring(x)||ischar(x));
addParameter(ip, 'ntracks',     10,          @(x) isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'combinexys',  true,        @islogical);
addParameter(ip, 'tmaxaftertx', [],          @(x) isempty(x)||isscalar(x));
addParameter(ip, 'tbeforetx',   [],          @(x) isempty(x)||isscalar(x));
addParameter(ip, 'smooth',      [],          @(x) isempty(x)||(isnumeric(x)&&isscalar(x)));
addParameter(ip, 'ylims',       {},          @(x) iscell(x)||isempty(x));
addParameter(ip, 'errortype',   'std',       @(x) any(strcmpi(x,{'std','sem','iqr'})));
addParameter(ip, 'font_size',   8,           @(x) isnumeric(x)&&isscalar(x));
addParameter(ip, 'exclude',     {},          @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'savedir',     '',          @(x) ischar(x)||isstring(x));
addParameter(ip, 'closefigs',   true,        @islogical);
addParameter(ip, 'track_alpha', 0.20,        @(x) isnumeric(x)&&isscalar(x));
parse(ip, varargin{:}); p = ip.Results;

% Normalize to row cell arrays
if ~iscell(p.channel); p.channel = {p.channel}; end
p.channel = cellfun(@char, p.channel(:)', 'UniformOutput', false);
p.facetby = cellstr(p.facetby); if iscolumn(p.facetby), p.facetby = p.facetby'; end
p.groupby = cellstr(p.groupby); if iscolumn(p.groupby), p.groupby = p.groupby'; end
p.exclude = cellstr(p.exclude); if iscolumn(p.exclude), p.exclude = p.exclude'; end

nChan  = numel(p.channel);
isDual = nChan > 1;

% Per-channel y-limits, padded to nChan
ylims = p.ylims;
if ~iscell(ylims), ylims = {}; end
while numel(ylims) < nChan, ylims{end+1} = []; end

%% ---- Timing --------------------------------------------------------------
if isfield(dataloc,'movieinfo') && ~isempty(dataloc.movieinfo.tsamp)
    looptime = dataloc.movieinfo.tsamp;
else
    looptime = 6;
    warning('plot_sorted_stacks: movieinfo.tsamp missing, assuming 6 min/frame.');
end
dt_hr = looptime / 60;

%% ---- Channel colors (left=blue, right=red) ------------------------------
chanColors = {[0.15 0.45 0.80]; [0.80 0.15 0.15]};
if nChan > 2, chanColors = num2cell(lines(nChan), 2); end

%% ---- Platemap levels -----------------------------------------------------
pmd = dataloc.platemapd.pmd;
facetFields = pss_resolve_fields(pmd, p.facetby, {'pTx','Tx1','Tx2'});
groupFields = pss_resolve_fields(pmd, p.groupby, {'Tx1','Tx2','pTx'});

facetLevels = pss_find_levels(pmd, p.facetby, facetFields, dt_hr);
if isempty(facetLevels)
    facetLevels = pss_find_levels(pmd, {''}, facetFields, dt_hr);
end
groupLevels = pss_find_levels(pmd, p.groupby, groupFields, dt_hr);
if isempty(groupLevels)
    groupLevels = pss_find_levels(pmd, {''}, groupFields, dt_hr);
end

% Exclusions
if ~isempty(p.exclude)
    exL = lower(p.exclude);
    facetLevels = facetLevels(~arrayfun(@(x) any(contains(lower(x.title),exL)), facetLevels));
    groupLevels = groupLevels(~arrayfun(@(x) any(contains(lower(x.title),exL)), groupLevels));
end

if isempty(facetLevels)
    warning('plot_sorted_stacks: no facet levels found after filtering.'); return;
end

%% ---- Global Tx anchor (hours) -------------------------------------------
allTimes = [[facetLevels.time_hr], [groupLevels.time_hr]];
allTimes = allTimes(isfinite(allTimes) & allTimes > 0);
tx_hr    = min([allTimes, NaN]);

grpCols = lines(max(1, numel(groupLevels)));

%% ---- Save directory ------------------------------------------------------
saveDir = char(p.savedir);
if isempty(saveDir), saveDir = fullfile(pwd, 'figures'); end
if ~exist(saveDir,'dir'), mkdir(saveDir); end

expName = '';
if isfield(dataloc,'file') && isfield(dataloc.file,'base')
    expName = strrep(dataloc.file.base,'_',' ');
end

%% =========================================================================
%% One figure per facet level
%% =========================================================================
for f = 1:numel(facetLevels)
    facetXY    = facetLevels(f).xys(:)';
    facetTitle = strrep(facetLevels(f).title,'_',' ');

    activeGroups = pss_build_active_groups(facetXY, groupLevels, grpCols, facetTitle);

    % grpData{g}{ch} = [nCells x T] inner-joined, sorted, subsampled
    % tAxes{g}       = [1 x T] hours (0 = tx)
    % grpCellIDs{g}  = [nCells x 1] cellindex values
    [grpData, tAxes, grpCellIDs] = pss_gather_all_groups( ...
        dataloc, activeGroups, p.channel, tx_hr, dt_hr, p);

    if ~any(cellfun(@(g) ~isempty(g{1}), grpData)), continue; end

    nGrp  = numel(activeGroups);
    nRows = p.ntracks + 1;   % row 1 = mean panel
    figW  = min(24, max(5, 3.5*nGrp));
    figH  = min(22, max(4, 0.9*nRows + 1.5));
    fig   = figure('Color','w','Units','inches','Position',[1 1 figW figH]);
    tl    = tiledlayout(nRows, nGrp, 'TileSpacing','compact', 'Padding','compact');
    sgtitle(tl, sprintf('%s  |  %s  |  %s', ...
        expName, strjoin(p.channel,' + '), facetTitle), 'Interpreter','none');

    %% ---- Row 1: mean +/- panel per group --------------------------------
    for g = 1:nGrp
        ax = nexttile(tl, g);
        hold(ax,'on'); box(ax,'on');
        tAxis = tAxes{g};
        if isempty(tAxis)
            title(ax, strrep(activeGroups(g).title,'_',' '), 'Interpreter','none');
            continue;
        end
        pss_plot_mean_panel(ax, grpData{g}, tAxis, p.channel, ...
            chanColors, ylims, p.errortype, isDual);
        if isfinite(tx_hr)
            xline(ax, 0, '--k', 'LineWidth',0.8, 'HandleVisibility','off');
        end
        title(ax, strrep(activeGroups(g).title,'_',' '), ...
            'Interpreter','none', 'FontSize', p.font_size+1);
        set(ax,'XTickLabel',{});
    end

    %% ---- Rows 2..nRows: one cell per row --------------------------------
    for iRow = 1:p.ntracks
        for g = 1:nGrp
            ax    = nexttile(tl, iRow*nGrp + g);
            tAxis = tAxes{g};
            nCells = size(grpData{g}{1}, 1);

            if isempty(tAxis) || iRow > nCells
                axis(ax,'off'); continue;
            end

            hold(ax,'on'); box(ax,'on');
            pss_plot_stack_row(ax, grpData{g}, iRow, tAxis, ...
                p.channel, chanColors, ylims, isDual);

            if isfinite(tx_hr)
                xline(ax, 0, '--k', 'LineWidth',0.6, 'HandleVisibility','off');
            end

            % Cell ID label
            if ~isempty(grpCellIDs{g}) && iRow <= numel(grpCellIDs{g})
                ylabel(ax, sprintf('c%d', grpCellIDs{g}(iRow)), ...
                    'FontSize', max(6,p.font_size-1), ...
                    'Rotation',0, 'HorizontalAlignment','right');
            end

            isLastRow = (iRow == p.ntracks) || (iRow == nCells);
            if isLastRow
                xlabel(ax, 'Time (hours)', 'FontSize', p.font_size);
            else
                set(ax,'XTickLabel',{});
            end
            set(ax,'YTick',[]);
        end
    end

    fontsize(fig, p.font_size, 'points');

    chanStr  = regexprep(strjoin(p.channel,'_'),'[^a-zA-Z0-9_]','_');
    facetStr = regexprep(facetLevels(f).title,'[^a-zA-Z0-9]','_');
    fname    = regexprep(sprintf('sorted_stacks_%s_%s',chanStr,facetStr),'_+','_');
    set(fig,'Renderer','painters');
    print(fig, fullfile(saveDir, fname), '-dsvg');
    fprintf('plot_sorted_stacks: saved %s\n', fullfile(saveDir,[fname '.svg']));
    if p.closefigs, close(fig); end
end
end % main

%% =========================================================================
%% PLOTTING HELPERS
%% =========================================================================

function pss_plot_mean_panel(ax, chanMats, tAxis, channels, ...
        chanColors, ylims, errortype, isDual)
    nChan = numel(channels);
    for ch = 1:nChan
        mat = chanMats{ch};
        if isempty(mat), continue; end
        col = chanColors{min(ch,numel(chanColors))};
        if isDual
            if ch==1, yyaxis(ax,'left');  set(ax,'YColor',col);
            else,     yyaxis(ax,'right'); set(ax,'YColor',col);
            end
        end
        mu       = mean(mat,1,'omitnan');
        [lo, hi] = pss_error_bounds(mat, errortype);
        xP = [tAxis, fliplr(tAxis)]; yP = [hi, fliplr(lo)];
        bad = isnan(xP)|isnan(yP);
        if ~all(bad)
            xP(bad)=[]; yP(bad)=[];
            patch(ax, xP, yP, col, 'FaceAlpha',0.20, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        plot(ax, tAxis, mu, '-', 'Color',col, 'LineWidth',1.8, ...
            'DisplayName', strrep(channels{ch},'_',' '));
        if ~isempty(ylims{ch}), ylim(ax, ylims{ch}); end
    end
    if isDual, yyaxis(ax,'left'); end
    legend(ax,'Location','best','FontSize',6);
end

function pss_plot_stack_row(ax, chanMats, iRow, tAxis, channels, ...
        chanColors, ylims, isDual)
    nChan = numel(channels);
    for ch = 1:nChan
        mat = chanMats{ch};
        if isempty(mat)||iRow>size(mat,1), continue; end
        col = chanColors{min(ch,numel(chanColors))};
        if isDual
            if ch==1, yyaxis(ax,'left');  set(ax,'YColor',col);
            else,     yyaxis(ax,'right'); set(ax,'YColor',col);
            end
        end
        plot(ax, tAxis, mat(iRow,:), '-', 'Color',col, ...
            'LineWidth',0.8, 'HandleVisibility','off');
        if ~isempty(ylims{ch}), ylim(ax, ylims{ch}); end
    end
    if isDual, yyaxis(ax,'left'); end
end

%% =========================================================================
%% DATA HELPERS
%% =========================================================================

function [grpData, tAxes, grpCellIDs] = pss_gather_all_groups( ...
        dataloc, activeGroups, channels, tx_hr, dt_hr, p)
    nG = numel(activeGroups);
    grpData=cell(1,nG); tAxes=cell(1,nG); grpCellIDs=cell(1,nG);
    for g = 1:nG
        xys = activeGroups(g).xys(:)';
        if ~p.combinexys, xys = xys(1); end

        rawMats=cell(1,numel(channels)); rawIDs=cell(1,numel(channels));
        tAxisRaw=[];
        for ch = 1:numel(channels)
            [mat,ids,tAx] = pss_gather_channel(dataloc, xys, channels{ch}, ...
                tx_hr, dt_hr, p.tmaxaftertx, p.tbeforetx);
            rawMats{ch}=mat; rawIDs{ch}=ids;
            if ch==1 && ~isempty(tAx), tAxisRaw=tAx; end
        end

        if isempty(rawMats{1})
            grpData{g}=repmat({[]},1,numel(channels));
            tAxes{g}=[]; grpCellIDs{g}=[]; continue;
        end

        % Inner join: intersect cellids across channels
        commonIDs = rawIDs{1};
        for ch = 2:numel(channels)
            if ~isempty(rawIDs{ch}), commonIDs=intersect(commonIDs,rawIDs{ch});
            else, commonIDs=[]; end
        end

        if isempty(commonIDs)
            grpData{g}=repmat({[]},1,numel(channels));
            tAxes{g}=[]; grpCellIDs{g}=[]; continue;
        end

        % Align rows to commonIDs order
        aligned=cell(1,numel(channels));
        for ch=1:numel(channels)
            [~,loc]=ismember(commonIDs, rawIDs{ch});
            aligned{ch}=rawMats{ch}(loc,:);
        end

        % Drop cells all-NaN in ANY channel
        keepMask=true(numel(commonIDs),1);
        for ch=1:numel(channels)
            keepMask=keepMask & ~all(isnan(aligned{ch}),2);
        end
        commonIDs=commonIDs(keepMask);
        for ch=1:numel(channels), aligned{ch}=aligned{ch}(keepMask,:); end

        if isempty(commonIDs)
            grpData{g}=repmat({[]},1,numel(channels));
            tAxes{g}=[]; grpCellIDs{g}=[]; continue;
        end

        % Sort by NaN count in ch1 (fewest first)
        [~,sIdx]=sort(sum(isnan(aligned{1}),2),'ascend');
        commonIDs=commonIDs(sIdx);
        for ch=1:numel(channels), aligned{ch}=aligned{ch}(sIdx,:); end

        % Subsample to ntracks
        nKeep=min(p.ntracks, size(aligned{1},1));
        commonIDs=commonIDs(1:nKeep);
        for ch=1:numel(channels), aligned{ch}=aligned{ch}(1:nKeep,:); end

        % Smooth
        if ~isempty(p.smooth)
            for ch=1:numel(channels)
                aligned{ch}=movmean(aligned{ch},p.smooth,2,'omitnan');
            end
        end

        grpData{g}=aligned; tAxes{g}=tAxisRaw; grpCellIDs{g}=commonIDs;
    end
end

function [mat, cellids, tAxis] = pss_gather_channel( ...
        dataloc, xys, chan, tx_hr, dt_hr, tmaxaftertx_hr, tbeforetx_hr)
    mat=[]; cellids=[]; tAxis=[];
    for xy = xys
        if xy<1||xy>numel(dataloc.d), continue; end
        S=dataloc.d{xy};
        if isempty(S)||~isstruct(S)||~isfield(S,'data'), continue; end
        if ~isfield(S.data,chan), continue; end
        if ~isfield(S,'cellindex')||isempty(S.cellindex), continue; end
        M=S.data.(chan); ids=S.cellindex(:);
        if isempty(M)||~isnumeric(M), continue; end
        if size(M,1)~=numel(ids)
            warning('plot_sorted_stacks: XY %d cellindex/data size mismatch for %s, skipping.',xy,chan);
            continue;
        end
        T1=size(mat,2); T2=size(M,2);
        if ~isempty(mat)
            if T1<T2, mat=[mat, NaN(size(mat,1),T2-T1)]; end %#ok<AGROW>
            if T2<T1, M  =[M,   NaN(size(M,1),  T1-T2)]; end
        end
        mat=[mat; M]; cellids=[cellids; ids]; %#ok<AGROW>
    end
    if isempty(mat), return; end
    [mat, tAxis] = pss_crop(mat, tx_hr, dt_hr, tmaxaftertx_hr, tbeforetx_hr);
end

%% =========================================================================
%% PLATEMAP HELPERS
%% =========================================================================

function groups = pss_build_active_groups(facetXY, groupLevels, grpCols, facetTitle)
    groups=struct('title',{},'xys',{},'col',{});
    for g=1:numel(groupLevels)
        xyg=intersect(facetXY, groupLevels(g).xys(:)');
        if ~isempty(xyg)
            groups(end+1).title=groupLevels(g).title; %#ok<AGROW>
            groups(end).xys=xyg; groups(end).col=grpCols(g,:);
        end
    end
    if isempty(groups)
        groups(1).title=facetTitle; groups(1).xys=facetXY; groups(1).col=[0.2 0.2 0.8];
    end
end

function fields = pss_resolve_fields(pmd, requested, defaults)
    fields={};
    for i=1:numel(requested)
        if ~isempty(requested{i})&&isfield(pmd,requested{i}), fields{end+1}=requested{i}; end %#ok<AGROW>
    end
    if isempty(fields)
        for i=1:numel(defaults)
            if isfield(pmd,defaults{i}), fields{end+1}=defaults{i}; end %#ok<AGROW>
        end
    end
end

function levels = pss_find_levels(pmd, tokens, fields, dt_hr)
    tokens=cellstr(tokens); [R,C]=size(pmd.xy); gMap=containers.Map();
    for t=1:numel(tokens)
        tok=tokens{t}; if isempty(tok), continue; end

        % '+' means AND: 'MK8722+Oligo' keeps only wells that received BOTH.
        subToks = strtrim(strsplit(tok,'+'));
        isCombo = numel(subToks) > 1;

        if isCombo
            % Collect XY sets per sub-token (across ALL fields), then intersect
            subXYsets = cell(1,numel(subToks));
            allTimes  = [];
            for st = 1:numel(subToks)
                subXYsets{st} = [];
                for fi=1:numel(fields)
                    fname=fields{fi}; if ~isfield(pmd,fname), continue; end
                    A=pmd.(fname);
                    for r=1:R; for c=1:C
                        nm=pss_safe_str(A,r,c,1);
                        if isempty(nm)||~contains(lower(nm),lower(subToks{st})), continue; end
                        xy_here=pmd.xy{r,c}; if isempty(xy_here), continue; end
                        subXYsets{st}=[subXYsets{st}, xy_here(:)'];
                        if any(strcmp(fname,{'Tx1','Tx2','Tx3','Tx4'}))
                            t_raw=pss_safe_num(A,r,c,4);
                            t_unit=lower(strtrim(pss_safe_str(A,r,c,5)));
                            if isfinite(t_raw)
                                switch t_unit
                                    case {'h','hr','hrs','hour','hours'}, allTimes(end+1)=t_raw;
                                    otherwise, allTimes(end+1)=t_raw*dt_hr;
                                end
                            end
                        end
                    end; end
                end
            end
            % AND-intersection
            comboXYs = subXYsets{1};
            for st=2:numel(subToks), comboXYs=intersect(comboXYs,subXYsets{st}); end
            if isempty(comboXYs), continue; end
            key    = strjoin(sort(lower(subToks)),'+');
            dispNm = strjoin(subToks,' + ');
            time_hr = min([allTimes, NaN]);
            if isKey(gMap,key)
                e=gMap(key); e.xys=[e.xys,comboXYs];
                if isfinite(time_hr), e.times(end+1)=time_hr; end
                gMap(key)=e;
            else
                e.title=dispNm; e.xys=comboXYs; e.dose=NaN; e.times=[];
                if isfinite(time_hr), e.times=time_hr; end
                gMap(key)=e;
            end

        else
            % Standard single-token search
            for fi=1:numel(fields)
                fname=fields{fi}; if ~isfield(pmd,fname), continue; end
                A=pmd.(fname);
                for r=1:R; for c=1:C
                    nm=pss_safe_str(A,r,c,1);
                    if isempty(nm)||~contains(lower(nm),lower(tok)), continue; end
                    dose=pss_safe_num(A,r,c,2); unit=pss_safe_str(A,r,c,3);
                    time_hr=NaN;
                    if any(strcmp(fname,{'Tx1','Tx2','Tx3','Tx4'}))
                        t_raw=pss_safe_num(A,r,c,4);
                        t_unit=lower(strtrim(pss_safe_str(A,r,c,5)));
                        if isfinite(t_raw)
                            switch t_unit
                                case {'h','hr','hrs','hour','hours'}, time_hr=t_raw;
                                otherwise, time_hr=t_raw*dt_hr;
                            end
                        end
                    end
                    if isnan(dose), key=sprintf('%s_NaN_%s',nm,unit); dispNm=nm;
                    else, key=sprintf('%s_%.6f_%s',nm,dose,unit);
                          dispNm=sprintf('%s %.4g %s',nm,dose,unit);
                    end
                    xy_here=pmd.xy{r,c}; if isempty(xy_here), continue; end
                    if isKey(gMap,key)
                        e=gMap(key); e.xys=[e.xys,xy_here(:)'];
                        if isfinite(time_hr), e.times(end+1)=time_hr; end
                        gMap(key)=e;
                    else
                        e.title=dispNm; e.xys=xy_here(:)'; e.dose=dose; e.times=[];
                        if isfinite(time_hr), e.times=time_hr; end
                        gMap(key)=e;
                    end
                end; end
            end
        end % isCombo
    end
    keys=gMap.keys;
    if isempty(keys), levels=struct('title',{},'xys',{},'dose',{},'time_hr',{}); return; end
    levels=repmat(struct('title','','xys',[],'dose',NaN,'time_hr',NaN),1,numel(keys));
    for k=1:numel(keys)
        e=gMap(keys{k}); levels(k).title=e.title; levels(k).xys=unique(e.xys);
        levels(k).dose=e.dose; levels(k).time_hr=min([e.times,NaN]);
    end
    doses=[levels.dose]; isNaN=isnan(doses);
    [~,I]=sort(doses(~isNaN)); order=[find(~isNaN(I)),find(isNaN)]; %#ok<FNDSB>
    if numel(order)==numel(levels), levels=levels(order); end
end

function [mat,tAxis] = pss_crop(mat, tx_hr, dt_hr, tmaxaftertx_hr, tbeforetx_hr)
    T=size(mat,2); txF=1;
    if isfinite(tx_hr), txF=max(1,min(T,round(tx_hr/dt_hr)+1)); end
    t1=1; if ~isempty(tbeforetx_hr)&&isfinite(tbeforetx_hr)
              t1=max(1,txF-round(tbeforetx_hr/dt_hr)); end
    t2=T; if ~isempty(tmaxaftertx_hr)&&isfinite(tmaxaftertx_hr)
              t2=min(T,txF+round(tmaxaftertx_hr/dt_hr)); end
    mat=mat(:,t1:t2); tAxis=((t1:t2)-txF)*dt_hr;
end

function [lo,hi] = pss_error_bounds(mat, errortype)
    mu=mean(mat,1,'omitnan');
    switch lower(errortype)
        case 'sem', n=sum(~isnan(mat),1); s=std(mat,0,1,'omitnan')./sqrt(max(n,1));
                    lo=mu-s; hi=mu+s;
        case 'iqr', pcts=prctile(mat,[25 75],1); lo=pcts(1,:); hi=pcts(2,:);
        otherwise,  s=std(mat,0,1,'omitnan'); lo=mu-s; hi=mu+s;
    end
end

function s = pss_safe_str(A,r,c,k)
    s=''; if k>size(A,3), return; end
    try; v=A{r,c,k};
        if ischar(v),   s=strtrim(v);       return; end
        if isstring(v), s=strtrim(char(v));  return; end
    catch; end
end

function x = pss_safe_num(A,r,c,k)
    x=NaN; if k>size(A,3), return; end
    try; v=A{r,c,k};
        if isnumeric(v)&&isscalar(v), x=double(v); return; end
        if ischar(v)||isstring(v)
            vv=str2double(strtrim(char(v))); if isfinite(vv), x=vv; end
        end
    catch; end
end