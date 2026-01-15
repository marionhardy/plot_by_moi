function osc_labels = label_oscillators_all(dataloc, chanName, maxCells, xyInclude, minTrackPts)
% label_oscillators_all  Manual triage across selected XYs for one channel.
%
% osc_labels = label_oscillators_all(dataloc, chanName, maxCells, xyInclude)
%
% Inputs:
%   dataloc   : your main dataloc struct with .d cell array
%   chanName  : string, field in S.data (e.g. 'HYLIGHT', 'Ch1', etc.)
%   maxCells  : (optional) max number of cells to label (e.g. 100–150)
%   xyInclude : (optional) list OR logical mask of xy indices to include
%
% Behaviour:
%   - Restricts to XYs in xyInclude (or all XYs if empty).
%   - For each cell, computes effective track length = # finite values
%     in S.data.(chanName)(row,:).
%   - Sorts cells by track length (descending) and shows the longest ones.
%   - Uses a fixed x-axis (time in hours) across all cells.
%   - Controls in the figure:
%       o = oscillating
%       n = not oscillating
%       u = unsure
%       b = back
%       q = quit
%
% NOTE: edit dt_hours below to match your imaging interval.

    if nargin < 3 || isempty(maxCells)
        maxCells = inf;
    end

    % ==== EDIT THIS LINE TO MATCH YOUR DATA ====
    % hours between frames (e.g. 10 min = 10/60 h)
    dt_hours = 10/60;  % <-- change this if needed
    % ===========================================

    nXY = numel(dataloc.d);

    % ---- build XY inclusion mask ----
    if nargin < 4 || isempty(xyInclude)
        xyMask = true(1, nXY);
    else
        if islogical(xyInclude)
            xyMask = xyInclude(:).';
            if numel(xyMask) ~= nXY
                error('xyInclude logical mask must have length numel(dataloc.d).');
            end
        else
            % assume numeric indices
            xyMask = false(1, nXY);
            xyInclude = xyInclude(:).';
            xyInclude = xyInclude(xyInclude>=1 & xyInclude<=nXY);
            xyMask(xyInclude) = true;
        end
    end

    % ---- build global list of cells from allowed XYs + track lengths ----
    cells   = struct('xy',{}, 'row',{}, 'cellindex',{});
    lengths = [];    % effective track length (finite points) per cell
    Tmax    = 0;     % global max number of timepoints across all candidate cells

    for xy = 1:nXY
        if ~xyMask(xy), continue; end
        S = dataloc.d{xy};
        if isempty(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        if ~isfield(S,'cellindex') || isempty(S.cellindex)
            continue;
        end

        X = S.data.(chanName);   % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end
        [nCellsXY, Txy] = size(X);
        Tmax = max(Tmax, Txy);   % update global max length

        for i = 1:nCellsXY
            tr = X(i,:);
            effLen = sum(isfinite(tr));      % number of finite timepoints
            if effLen < minTrackPts
                continue;                    % skip truncated tracks
            end
        
            cells(end+1).xy        = xy;          
            cells(end).row         = i;
            cells(end).cellindex   = S.cellindex(i);
            lengths(end+1,1)       = effLen;
        end
    end

    n = numel(cells);
    if n == 0
        warning('No cells found for channel %s in the selected XYs.', chanName);
        osc_labels = struct('cells',cells,'labels',[],'chanName',chanName);
        return;
    end

    % ---- sort cells by track length (longest first) ----
    [~, order] = sort(lengths, 'descend');
    cells   = cells(order);
    lengths = lengths(order); %#ok<NASGU>

    % ---- restrict to the top maxCells ----
    nUse   = min(n, maxCells);
    cells  = cells(1:nUse);
    labels = nan(nUse,1);  % 1=osc, 0=not, -1=unsure

    % ---- create figure, wider and white background ----
    hFig = figure('Name','Oscillation triage', ...
                  'NumberTitle','off', ...
                  'Color','w');

    set(hFig, 'Units','pixels');
    pos = get(hFig,'Position');   % [left bottom width height]
    pos(3) = 1000;                % width
    pos(4) = 400;                 % height
    set(hFig,'Position',pos);

    % precompute global x-axis limit in hours
    if Tmax < 1
        xMaxHours = dt_hours;   % avoid zero range
    else
        xMaxHours = (Tmax-1) * dt_hours;
    end

    k = 1;
    while k>=1 && k<=nUse && ishandle(hFig)
        c = cells(k);
        S = dataloc.d{c.xy};
        tr = S.data.(chanName)(c.row,:);
        Tcell   = numel(tr);
        t_hours = (0:Tcell-1) * dt_hours;

        figure(hFig);  % ensure this figure has focus
        clf(hFig);
        ax = axes('Parent',hFig);
        hold(ax,'on');

        % main raw trace: thicker line
        plot(ax, t_hours, tr, '-', 'LineWidth', 2);

        % optional smoothed trace as dotted overlay
        if Tcell >= 3
            trSmooth = movmean(tr, 3);
            plot(ax, t_hours, trSmooth, ':', 'LineWidth', 1);
        end

        % fixed x-axis across all cells (in hours)
        xlim(ax, [0, xMaxHours]);

        hold(ax,'off');
        set(ax, 'Box','on', 'FontSize',10);

        title(ax, sprintf('Cell %d/%d    xy=%d  cellindex=%d', ...
              k, nUse, c.xy, c.cellindex));
        xlabel(ax,'time (h)');
        ylabel(ax,chanName);

        if ~isnan(labels(k))
            txt = 'current label: ';
            switch labels(k)
                case 1, labtxt = 'oscillating';
                case 0, labtxt = 'not oscillating';
                case -1, labtxt = 'unsure';
                otherwise, labtxt = 'unlabeled';
            end
            text(ax, 0.02, 0.95, [txt labtxt], 'Units','normalized');
        end

        % instructions overlay
        text(ax, 0.02, 0.05, 'click in this figure, then: o=osc, n=not, u=unsure, b=back, q=quit', ...
             'Units','normalized','FontSize',9);

        % ---- wait for a *valid* key from this figure ----
        key = '';
        while ishandle(hFig) && isempty(key)
            % make sure the figure is active; then wait
            figure(hFig);
            waitforbuttonpress;
            if ~ishandle(hFig), break; end
            key = get(hFig,'CurrentCharacter');  % single char
            if isempty(key), key = ''; continue; end
            key = lower(key);
            if ~ismember(key, ['o','n','u','b','q'])
                key = '';  % ignore other keys (enter, arrows, etc.)
            end
        end
        if ~ishandle(hFig), break; end

        switch key
            case 'o'
                labels(k) = 1; k = k+1;
            case 'n'
                labels(k) = 0; k = k+1;
            case 'u'
                labels(k) = -1; k = k+1;
            case 'b'
                if k > 1, k = k-1; end
            case 'q'
                break;
            otherwise
                % shouldn't happen, but just stay on same cell
        end
    end

    if ishandle(hFig)
        close(hFig);
    end

    osc_labels = struct('cells',cells, ...
                        'labels',labels, ...
                        'chanName',chanName);
end
