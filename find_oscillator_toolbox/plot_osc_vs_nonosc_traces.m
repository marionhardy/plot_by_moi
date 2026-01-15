function plot_osc_vs_nonosc_traces(dataloc, chanName, dt, varargin)
% plot_osc_vs_nonosc_traces  Overlay traces by oscillation class (0/1/2),
% aligned to Tx1 (t=0) and cropped to a fixed time window.
%
% plot_osc_vs_nonosc_traces(dataloc, chanName, dt, ...)
%
% Inputs:
%   dataloc   : struct with .d{xy}.data.(chanName) and classification fields
%   chanName  : channel name (e.g. 'HYLIGHT')
%   dt        : sampling interval in minutes
%
% Name-value options:
%   'maxNon'     : max # class-0 cells to plot (default 50)
%   'maxMed'     : max # class-1 cells to plot (default 50)
%   'maxHigh'    : max # class-2 cells to plot (default 50)
%   'maxCells'   : convenience: set all three maxima to this value
%                  (per-class cap; e.g. 50 -> 50 non, 50 med, 50 high)
%   'xyInclude'  : vector of xy indices OR logical mask to include (default: all)
%   'useBinaryIfNoClass' : if oscClass missing, use isOsc (default true)
%   'binaryOscAs' : if using binary isOsc, map true->2 (high) or true->1 (med) (default 2)
%
%   Alignment + cropping:
%   'tPre'       : hours BEFORE Tx1 to include (default 4)
%   'tPost'      : hours AFTER  Tx1 to include (default 12)
%   'txField'    : which platemap field to use for alignment (default 'Tx1')
%   'txSlice'    : which slice within pmd.(txField) stores time (default 4)
%
% Behavior:
%   - For each XY, find Tx1 frame from dataloc.platemapd.pmd.(txField){r,c,txSlice}
%     using pmd.xy to map xy -> (row,col).
%   - Crop each cell trace to [Tx1 - tPre, Tx1 + tPost] and shift time so Tx1 = 0 h.
%   - Keep longest tracks first per class when subsampling.
%   - Plot 3 subplots: class0 (non), class1 (medium), class2 (high/perfect).
%
% Notes:
%   - If you only have S.osc.(chanName).isOsc (binary), you won't truly get
%     class 1 vs 2 separation unless you have oscClass saved. We provide a
%     mapping option for convenience.

    ip = inputParser; ip.CaseSensitive = false;
    addParameter(ip,'maxNon',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxMed',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxHigh',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxCells',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(ip,'xyInclude',[],@(x)isnumeric(x) || islogical(x));
    addParameter(ip,'useBinaryIfNoClass',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'binaryOscAs',uint8(2),@(x)isnumeric(x)&&isscalar(x)&&any(x==[1 2]));

    % Tx alignment + cropping
    addParameter(ip,'tPre',4,@(x)isnumeric(x)&&isscalar(x)&&x>=0);     % hours before Tx1
    addParameter(ip,'tPost',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % hours after Tx1
    addParameter(ip,'txField','Tx1',@(s)ischar(s) || isstring(s));    % pmd field
    addParameter(ip,'txSlice',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1); % slice storing time/tp

    parse(ip,varargin{:});
    P = ip.Results;

    % ---- convenience: maxCells means "per class" ----
    if ~isempty(P.maxCells)
        P.maxNon  = P.maxCells;
        P.maxMed  = P.maxCells;
        P.maxHigh = P.maxCells;
    end

    dt_hours = dt/60;

    % ---- build XY inclusion mask ----
    nXY = numel(dataloc.d);
    if isempty(P.xyInclude)
        xyMask = true(1,nXY);
    else
        if islogical(P.xyInclude)
            xyMask = P.xyInclude(:).';
            if numel(xyMask) ~= nXY
                error('xyInclude logical mask must have length numel(dataloc.d).');
            end
        else
            xyMask = false(1,nXY);
            idx = P.xyInclude(:).';
            idx = idx(idx>=1 & idx<=nXY);
            xyMask(idx) = true;
        end
    end

    % ---- prepare mapping xy -> (row,col) using pmd.xy ----
    hasPMD = isfield(dataloc,'platemapd') && isfield(dataloc.platemapd,'pmd') ...
          && isfield(dataloc.platemapd.pmd,'xy');
    if hasPMD
        pmd = dataloc.platemapd.pmd;
        xyGrid = pmd.xy; % cell [8x12], entries like numeric xy
    else
        pmd = [];
        xyGrid = [];
    end

    % ---- collect traces and lengths ----
    tr0 = {}; tr1 = {}; tr2 = {};
    len0 = []; len1 = []; len2 = [];
    Tmax = 0;  % max length after cropping

    for xy = 1:nXY
        if ~xyMask(xy), continue; end
        S = dataloc.d{xy};
        if isempty(S) || ~isstruct(S) || ~isfield(S,'data') || ~isfield(S.data,chanName)
            continue;
        end
        if ~isfield(S,'osc') || ~isfield(S.osc,chanName)
            continue;
        end

        X = S.data.(chanName); % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end

        % Prefer oscClass if present
        hasClass = isfield(S.osc.(chanName),'oscClass');
        hasBin   = isfield(S.osc.(chanName),'isOsc');

        if hasClass
            cls = S.osc.(chanName).oscClass(:);
            cls = uint8(cls);
        elseif P.useBinaryIfNoClass && hasBin
            isOsc = S.osc.(chanName).isOsc(:);
            cls = zeros(size(isOsc), 'uint8');
            cls(isOsc) = uint8(P.binaryOscAs); % map binary osc to class 2 (default) or class 1
        else
            continue;
        end

        [nCellsXY,Txy] = size(X);
        if numel(cls) ~= nCellsXY
            warning('XY %d: class label length mismatch; skipping.', xy);
            continue;
        end

        % ---- find Tx frame for this XY from pmd ----
        TxFrame = 1; % fallback
        if hasPMD && isfield(pmd, char(P.txField))
            try
                % locate row/col for this xy
                [r,c] = find(cellfun(@(v) isequal(v,xy), xyGrid));
                if ~isempty(r)
                    txArr = pmd.(char(P.txField));
                    if ndims(txArr) >= 3
                        tVal = txArr{r(1),c(1),P.txSlice};
                    else
                        tVal = txArr{r(1),c(1)};
                    end
                    if ischar(tVal) || isstring(tVal)
                        tVal = str2double(tVal);
                    end
                    if isfinite(tVal)
                        TxFrame = tVal;
                    end
                end
            catch
                % keep fallback TxFrame=1
            end
        end

        if isempty(TxFrame) || ~isfinite(TxFrame) || TxFrame < 1
            TxFrame = 1;
        end

        % ---- define crop window indices based on TxFrame ----
        % Tx time in hours relative to movie frame 1:
        t0_hr = (TxFrame - 1) * dt_hours;

        % desired window in hours relative to Tx: [-tPre, +tPost]
        tStart_hr = -P.tPre;
        tEnd_hr   =  P.tPost;

        % convert desired absolute hours to indices
        iStart = max(1, round((t0_hr + tStart_hr)/dt_hours) + 1);
        iEnd   = min(Txy, round((t0_hr + tEnd_hr)/dt_hours) + 1);

        if iEnd <= iStart
            continue;
        end

        % update Tmax using CROPPED length
        Tmax = max(Tmax, iEnd - iStart + 1);

        % ---- collect per-cell traces (CROPPED) ----
        for i = 1:nCellsXY
            tr = X(i, iStart:iEnd);
            if all(~isfinite(tr)), continue; end

            effLen = sum(isfinite(tr));
            if effLen < 5, continue; end

            switch cls(i)
                case 0
                    tr0{end+1} = tr;        %#ok<AGROW>
                    len0(end+1) = effLen;   %#ok<AGROW>
                case 1
                    tr1{end+1} = tr;        %#ok<AGROW>
                    len1(end+1) = effLen;   %#ok<AGROW>
                case 2
                    tr2{end+1} = tr;        %#ok<AGROW>
                    len2(end+1) = effLen;   %#ok<AGROW>
                otherwise
                    % ignore
            end
        end
    end

    n0 = numel(tr0);
    n1 = numel(tr1);
    n2 = numel(tr2);

    if n0==0 && n1==0 && n2==0
        warning('No traces with oscillation labels found for %s.', chanName);
        return;
    end

    fprintf('Found class0(non)=%d, class1(med)=%d, class2(high)=%d traces for %s.\n', ...
        n0, n1, n2, chanName);

    % ---- subsample by longest length if too many ----
    if n0 > P.maxNon
        [~, ord] = sort(len0, 'descend');
        keep = ord(1:P.maxNon);
        tr0 = tr0(keep); len0 = len0(keep);
        n0 = numel(tr0);
    end
    if n1 > P.maxMed
        [~, ord] = sort(len1, 'descend');
        keep = ord(1:P.maxMed);
        tr1 = tr1(keep); len1 = len1(keep);
        n1 = numel(tr1);
    end
    if n2 > P.maxHigh
        [~, ord] = sort(len2, 'descend');
        keep = ord(1:P.maxHigh);
        tr2 = tr2(keep); len2 = len2(keep);
        n2 = numel(tr2);
    end

    % ---- helper to build a padded matrix (NaN for missing timepoints) ----
    function M = buildMatrix(trCell, TmaxLocal)
        M = NaN(numel(trCell), TmaxLocal);
        for ii = 1:numel(trCell)
            v = trCell{ii};
            L = min(numel(v), TmaxLocal);
            M(ii,1:L) = v(1:L);
        end
    end

    M0 = []; M1 = []; M2 = [];
    if n0 > 0, M0 = buildMatrix(tr0, Tmax); end
    if n1 > 0, M1 = buildMatrix(tr1, Tmax); end
    if n2 > 0, M2 = buildMatrix(tr2, Tmax); end

    % ---- time axis: Tx is 0 ----
    th = linspace(-P.tPre, P.tPost, Tmax);

    % ---- plotting ----
    figure('Name',sprintf('Osc classes: %s (Tx-aligned)', chanName), ...
           'NumberTitle','off','Color','w');

    % Class 0: non
    subplot(3,1,1); hold on;
    if ~isempty(M0)
        for i = 1:size(M0,1)
            plot(th, M0(i,:), 'Color',[0.7 0.7 1]); % light blue
        end
        plot(th, nanmean(M0,1), 'b-', 'LineWidth', 2);
    end
    xline(0,'k-','LineWidth',1);
    hold off; xlim([th(1) th(end)]);
    ylabel(chanName);
    title(sprintf('Class 0: Non-osc (n = %d)', n0));
    box on;

    % Class 1: medium
    subplot(3,1,2); hold on;
    if ~isempty(M1)
        for i = 1:size(M1,1)
            plot(th, M1(i,:), 'Color',[0.8 0.8 0.8]); % light gray
        end
        plot(th, nanmean(M1,1), 'k-', 'LineWidth', 2);
    end
    xline(0,'k-','LineWidth',1);
    hold off; xlim([th(1) th(end)]);
    ylabel(chanName);
    title(sprintf('Class 1: Medium (n = %d)', n1));
    box on;

    % Class 2: high
    subplot(3,1,3); hold on;
    if ~isempty(M2)
        for i = 1:size(M2,1)
            plot(th, M2(i,:), 'Color',[1 0.7 0.7]); % light red
        end
        plot(th, nanmean(M2,1), 'r-', 'LineWidth', 2);
    end
    xline(0,'k-','LineWidth',1);
    hold off; xlim([th(1) th(end)]);
    xlabel(sprintf('time (h) relative to %s', char(P.txField)));
    ylabel(chanName);
    title(sprintf('Class 2: Strong (n = %d)', n2));
    box on;
end
