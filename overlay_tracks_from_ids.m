function overlay_tracks_from_ids(dataloc, GroupStruct, varargin)
% Overlay single-cell tracks for Osc vs Non with optional per-XY treatment filtering.
% X-axis is HOURS, with 0 at the selected treatment time parsed from metadata.
% GroupStruct must have fields: .Oscillatory [ids], .NonOscillatory [ids]

% ---------- args ----------
ip = inputParser; ip.CaseSensitive=false;
addParameter(ip,'channel','HYLIGHT');
addParameter(ip,'group','osc',@(s) any(strcmpi(s,{'osc','non','both'})));
addParameter(ip,'ntracks',20,@(x) isnumeric(x)&&isscalar(x)&&x>0);  % total when group='both', per-group otherwise
addParameter(ip,'smooth',[],@(x) isempty(x)||isscalar(x));
addParameter(ip,'tile',true,@islogical);
addParameter(ip,'tmax_hours',[],@(x) isempty(x)||isscalar(x));      % crop window length (hours from movie start)
addParameter(ip,'tx_filter',{},@(x) iscell(x) || isstring(x) || ischar(x));
addParameter(ip,'tx_mode','any',@(s) any(strcmpi(s,{'any','all','none'})));
addParameter(ip,'show_std',true,@islogical);     % mean ± SD shading
addParameter(ip,'std_alpha',0.18,@isscalar);
addParameter(ip,'invert_layers',true,@islogical);% plot Non first, Osc last
addParameter(ip,'show_mean',true,@islogical);    % draw mean line + legend handles
addParameter(ip,'full_only',false,@islogical);   % do NOT enforce full tracks by default
addParameter(ip,'balance_groups',true,@islogical);% when 'both', sample ntracks/2 per group
addParameter(ip,'frame_minutes',[],@(x) isempty(x)||isscalar(x));   % minutes per frame (e.g., 6)
addParameter(ip,'anchor_tx',[]); % (unused now; anchor is from metadata)
parse(ip,varargin{:}); p=ip.Results; chan=char(p.channel);

% colors (Non behind/red, Osc on top/blue)
colOsc = [0.15 0.45 0.80];
colNon = [0.75 0.20 0.25];

% minutes per frame (priority: user arg > dataloc.meta > default 6)
if ~isempty(p.frame_minutes)
    mPerFrame = p.frame_minutes;
elseif isfield(dataloc,'movieinfo') && isfield(dataloc.movieinfo,'tsamp') && ~isempty(dataloc.movieinfo.tsamp)
    mPerFrame = dataloc.movieinfo.tsamp;
else
    mPerFrame = 6; % fallback
end
framesPerHour = 60/mPerFrame;

% experiment name for figure title
expName = '';
if isfield(dataloc,'file') && isfield(dataloc.file,'base') && ~isempty(dataloc.file.base)
    expName = strrep(dataloc.file.base,'_',' ');
end

% ---------- XYs with data ----------
goodXY = find(cellfun(@(x) ~isempty(x) && isfield(x,'cellindex') && ~isempty(x.cellindex), dataloc.d));
xyList = goodXY(:);

% ---------- optional Tx filtering (per-XY) ----------
if ~isempty(p.tx_filter)
    want = lower(string(p.tx_filter));
    keep = false(size(xyList));
    for i = 1:numel(xyList)
        tag = lower(local_xy_tx_compact(dataloc, xyList(i), {})); % uses Tx1..Tx4 only
        has = arrayfun(@(w) contains(tag, w), want);
        switch lower(p.tx_mode)
            case 'all',  keep(i) = all(has);
            case 'any',  keep(i) = any(has);
            case 'none', keep(i) = ~any(has);
        end
    end
    xyList = xyList(keep);
end

% ---------- plotting ----------
if p.tile
    fig = figure('Color','w','Position',[100 100 1600 900]);
    tl = tiledlayout(fig,'flow','TileSpacing','compact','Padding','compact');
    if ~isempty(expName), sgtitle(tl, expName); end
else
    fig = [];
end

for xy = xyList(:)'
    D = dataloc.d{xy};
    if ~isfield(D,'data') || ~isfield(D.data,chan), continue; end
    X = D.data.(chan);
    T = size(X,2);
    if ~isempty(p.tmax_hours), T = min(T, floor(p.tmax_hours * framesPerHour)); end
    if T < 2, continue; end
    cellIDs_xy = D.cellindex(:);

    % ---- treatment anchor (frame) from metadata + title string WITH DOSES
    [t0_frame, txLabelShort] = local_pick_anchor_frame_and_label(dataloc, xy, p.tx_filter);
    if isempty(t0_frame) || ~isfinite(t0_frame), t0_frame = 1; end
    tH = ((1:T) - t0_frame) / framesPerHour; % hours, t=0 at first matching Tx

    % group order controls draw stack
    switch lower(p.group)
        case 'osc'
            groups = {{'Oscillatory',    GroupStruct.Oscillatory(:),    colOsc}};
        case 'non'
            groups = {{'NonOscillatory', GroupStruct.NonOscillatory(:), colNon}};
        case 'both'
            if p.invert_layers
                groups = {{'NonOscillatory', GroupStruct.NonOscillatory(:), colNon}; ...
                          {'Oscillatory',    GroupStruct.Oscillatory(:),    colOsc}};
            else
                groups = {{'Oscillatory',    GroupStruct.Oscillatory(:),    colOsc}; ...
                          {'NonOscillatory', GroupStruct.NonOscillatory(:), colNon}};
            end
    end

    if p.tile
        ax = nexttile(tl);
    else
        figure('Color','w','Position',[100 100 1600 900]); ax = gca;
    end
    hold(ax,'on'); plotted=false; hMeans = gobjects(0);

    % balanced sampling logic
    if strcmpi(p.group,'both') && p.balance_groups
        nPer = max(1, floor(p.ntracks/2));
    else
        nPer = p.ntracks;
    end

    for g = 1:numel(groups)
        gname = groups{g}{1}; IDs = groups{g}{2}; col = groups{g}{3};
        if isempty(IDs), continue; end

        rowsAll = find(ismember(cellIDs_xy, IDs));
        if isempty(rowsAll), continue; end

        rowsKeep = rowsAll;
        if p.full_only
            rowsKeep = rowsAll(~any(isnan(X(rowsAll,1:T)),2));
            if isempty(rowsKeep), continue; end
        end

        if numel(rowsKeep) > nPer, rows = rowsKeep(randperm(numel(rowsKeep), nPer));
        else,                     rows = rowsKeep;
        end

        dat = X(rows, 1:T);
        if ~isempty(p.smooth), dat = movmean(dat, p.smooth, 2, 'omitnan'); end

        h = plot(ax, tH, dat', 'LineWidth', 0.6, 'Color', col, 'HandleVisibility','off');
        local_setalpha(ax, h, 0.75, col);

        m = mean(dat,1,'omitnan');
        if p.show_std && size(dat,1) >= 2
            s = std(dat,0,1,'omitnan');
            patch(ax, [tH, fliplr(tH)], [m - s, fliplr(m + s)], col, ...
                  'FaceAlpha', p.std_alpha, 'EdgeColor','none', 'HandleVisibility','off');
        end
        if p.show_mean
            hMeans(end+1) = plot(ax, tH, m, '-', 'LineWidth', 1.2, 'Color', col, 'DisplayName', gname); %#ok<AGROW>
        end
        plotted = true;
    end

    if ~plotted, delete(ax); continue; end

    xline(ax, 0, '--', 'Color', [0 0 0 0.6], 'HandleVisibility','off');
    title(ax, ['XY ' char(string(xy)) ' | ' char(txLabelShort)]);
    xlabel(ax,'hours'); ylabel(ax, strrep(chan,'_',' '));

    L = gobjects(0,1);
    if p.show_mean
        L = hMeans(isgraphics(hMeans));
    else
        for g = 1:numel(groups)
            gname = groups{g}{1}; col = groups{g}{3};
            L(end+1) = plot(ax, NaN, NaN, '-', 'LineWidth', 1.2, ...
                            'Color', col, 'DisplayName', gname); %#ok<AGROW>
        end
    end
    if any(isgraphics(L)), legend(ax, L, 'Location','best'); end
end
end

% =================== helpers ===================

function [t0_frame, titleLabel] = local_pick_anchor_frame_and_label(dataloc, xy, tx_filter)
% Parse Tx details (5-plane blocks). Choose earliest matching frame; build label (name + dose).
tx = local_tx_details_for_xy(dataloc, xy); % struct array: .name .dose .frame
parts = arrayfun(@(t) strtrim(strjoin({t.name, t.dose}, ' ')), tx, 'UniformOutput', false);
titleLabel = strtrim(regexprep(strjoin(unique(parts(~cellfun(@isempty,parts))), ' + '),'\s+',' '));

t0_frame = [];
if isempty(tx), return; end
if ~isempty(tx_filter)
    want = lower(string(tx_filter));
    sel  = arrayfun(@(t) any(contains(lower(string(t.name)), want)), tx);
    frSel = [tx(sel).frame]; frSel = frSel(isfinite(frSel) & frSel>0);
    if ~isempty(frSel), t0_frame = min(frSel); return; end
end
frAll = [tx.frame]; frAll = frAll(isfinite(frAll) & frAll>0);
if ~isempty(frAll), t0_frame = min(frAll); end
end

function tx = local_tx_details_for_xy(dataloc, xy)
% Return struct array .name .dose .frame for this XY using 5-plane blocks:
% per TxX: [name(1), dose(2), name2(3), time(4), spare(5)] × numTreatments
tx = struct('name',{},'dose',{},'frame',{});
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'pmd'), return; end
pmd = dataloc.platemapd.pmd;
[r,c] = local_find_rc_for_xy(pmd, xy); if isnan(r), return; end

txFields = {'Tx1','Tx2','Tx3','Tx4'}; % <-- ONLY Tx fields (no pTx)
for s = 1:numel(txFields)
    f = txFields{s}; if ~isfield(pmd,f) || isempty(pmd.(f)), continue; end
    A = pmd.(f); kmax = size(A,3); if kmax==0, continue; end
    nblk = floor(kmax/5);
    for b = 1:nblk
        base = 5*(b-1);
        nmParts = strings(1,0);
        for pp = 1:min(3,kmax-base)
            nmParts(end+1) = string(local_pull(A, r, c, base+pp)); %#ok<AGROW>
        end
        nm = strtrim(strjoin(nmParts(nmParts~=""), " "));
        dose = ''; if base+2 <= kmax, dose = strtrim(char(string(local_pull(A, r, c, base+2)))); end
        frm  = NaN;
        if base+4 <= kmax
            v = A(r,c,base+4); if iscell(v), v=v{1}; end
            if isnumeric(v) && isscalar(v) && isfinite(v) && v>0, frm = double(v); end
        end
        if (~isempty(nm) || ~isempty(dose) || (isfinite(frm)&&frm>0))
            tx(end+1) = struct('name',char(nm),'dose',dose,'frame',frm); %#ok<AGROW>
        end
    end
end

if ~isempty(tx)
    fr = [tx.frame];
    [~,I] = sortrows([~isfinite(fr(:)), fr(:)], [1 2]);
    tx = tx(I);
end
end

function label = local_xy_tx_compact(dataloc, xy, tx_filter)
% Build compact title like "Oligo 3 uM + MK8722 10 uM" from Tx1..Tx4 only.
label = '';
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'pmd'), return; end
pmd = dataloc.platemapd.pmd;
[r,c] = local_find_rc_for_xy(pmd, xy); if isnan(r), return; end

txFields = {'Tx1','Tx2','Tx3','Tx4'};
parts = {};
isUnit = @(s) ~isempty(regexpi(s,'^(uM|µM|nM|mM|pM|%|ug/mL|µg/mL|ng/mL|mg/mL)$','once'));
hasLetters = @(s) any(isletter(s));
trimstr = @(s) strtrim(char(string(s)));

for s = 1:numel(txFields)
    f = txFields{s};
    if ~isfield(pmd,f) || isempty(pmd.(f)), continue; end
    A = pmd.(f); kmax = size(A,3); if kmax==0, continue; end
    nblk = max(1, floor(kmax/5));           % handle 2-plane and 5-plane layouts
    for b = 1:nblk
        base = 5*(b-1);
        toks = cell(1,3);
        for pp = 1:min(3,kmax-base)
            toks{pp} = trimstr(local_pull(A,r,c,base+pp));
        end
        toks = toks(~cellfun(@isempty,toks));

        if isempty(toks), continue; end

        % 1) NAME: first token with letters that is NOT a pure unit
        name = '';
        for t = 1:numel(toks)
            if hasLetters(toks{t}) && ~isUnit(toks{t})
                name = toks{t}; break;
            end
        end
        if isempty(name), continue; end

        % 2) DOSE:
        dose = '';
        % 2a) Look for combined "number+unit" in any token
        gotCombined = false;
        for t = 1:numel(toks)
            m = regexp(toks{t},'^\s*(\d+(\.\d+)?)\s*([a-zA-Zµ/%].*)$','once','tokens');
            if ~isempty(m)
                dose = [m{1}{1} ' ' strtrim(m{1}{3})];
                gotCombined = true; break;
            end
        end
        % 2b) Else pair last numeric with a unit token
        if ~gotCombined
            numTok = '';
            unitTok = '';
            for t = 1:numel(toks)
                if ~isempty(regexp(toks{t},'^\s*\d+(\.\d+)?\s*$', 'once'))
                    numTok = toks{t};
                elseif isUnit(toks{t})
                    unitTok = toks{t};
                end
            end
            if ~isempty(numTok) && ~isempty(unitTok)
                dose = [numTok ' ' unitTok];
            elseif ~isempty(numTok)
                dose = numTok; % no unit available
            end
        end

        if ~isempty(dose), parts{end+1} = [name ' ' dose]; %#ok<AGROW>
        else,              parts{end+1} = name;            %#ok<AGROW>
        end
    end
end

% Optional filter (e.g., {'oligo','MK8722'})
if ~isempty(tx_filter)
    want = lower(string(tx_filter));
    keep = cellfun(@(p) any(contains(lower(p), want)), parts);
    parts = parts(keep);
end

% Deduplicate and join
if ~isempty(parts)
    parts = unique(parts,'stable');
    label = regexprep(strjoin(parts,' + '), '\s+', ' ');
end
end

function [r,c] = local_find_rc_for_xy(pmd, xy)
% Find row/col in platemap for this XY (robust to different xy storage)
r = NaN; c = NaN;
try
    [H,W,~] = size(pmd.Cell);
catch
    return;
end
for rr=1:H
    for cc=1:W
        xys = [];
        try
            xys = pmd.xy{rr,cc};
        catch
            k = sub2ind(size(pmd.Cell), rr, cc);
            if isfield(pmd,'xy') && k<=numel(pmd.xy), xys = pmd.xy{k}; end
        end
        if isempty(xys), continue; end
        if any(isfinite(xys) & (round(double(xys)) == xy)), r=rr; c=cc; return; end
    end
end
end

function val = local_pull(A, r, c, plane)
% Safe pull from 3D platemap array; returns '' for empty; numerics -> char
val = '';
if isempty(A) || r<1 || c<1 || r>size(A,1) || c>size(A,2) || plane>size(A,3), return; end
v = A(r,c,plane);
if iscell(v), v = v{1}; end
if isempty(v), return; end
if isstring(v), v = char(v); end
if ischar(v), val = strtrim(v); return; end
if isnumeric(v)
    vv = v(:); vv = vv(isfinite(vv));
    if isempty(vv), val=''; return; end
    if abs(vv(1)) > 1e6, val = sprintf('%.3g', vv(1));
    else,               val = sprintf('%g',   vv(1));
    end
end
end

function local_setalpha(ax, hLines, a, col)
% Try RGBA on line Color; otherwise blend with axes color as fallback.
try
    for i=1:numel(hLines)
        if isgraphics(hLines(i)), hLines(i).Color = [col, a]; end
    end
catch
    bg = get(ax,'Color'); if ischar(bg) || (isstring(bg) && strlength(bg)>0), bg = [1 1 1]; end
    blend = @(c) (1-a)*bg + a*c;
    for i=1:numel(hLines)
        if isgraphics(hLines(i)), hLines(i).Color = blend(col); end
    end
end
end
