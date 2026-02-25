function plot_osc_matrix(dataloc, chanName, dt, varargin)
% plot_osc_matrix  Bulk view of traces by oscillation class (0/1/2),
% with Tx-alignment + cropping, plus peak raster.
%
% DEFAULT BEHAVIOR (updated):
%   - Prioritize FULL-window traces (i.e., traces that cover the entire requested
%     crop window [tPre, tPost] around Tx, or the full XY trace if not aligning).
%   - If not enough full traces to hit maxNon/maxMed/maxHigh, backfill with the
%     LONGEST-AVAILABLE partial traces.
%   - Still requires minimum finite content and a finite-fraction criterion.
%
% Output:
%   Top:    heatmap of ROW-SCALED traces
%   Bottom: peak raster
%
% Class meaning (recommended):
%   0 = non-oscillatory
%   1 = weak/medium oscillatory
%   2 = strong/perfect oscillatory

    ip = inputParser; ip.CaseSensitive = false;

    addParameter(ip,'xyInclude',[],@(x)isnumeric(x) || islogical(x));
    addParameter(ip,'classChan',chanName,@(s)ischar(s)||isstring(s));

    % per-class caps
    addParameter(ip,'maxNon',100,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxMed',100,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxHigh',100,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'maxCells',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));

    % class source
    addParameter(ip,'useBinaryIfNoClass',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'binaryOscAs',uint8(2),@(x)isnumeric(x)&&isscalar(x)&&any(x==[1 2]));

    % peak detection
    addParameter(ip,'Pmin',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);   % hours
    addParameter(ip,'Pmax',6,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'smoothHrs',[],@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'promFactor',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>=0);

    % alignment/cropping
    addParameter(ip,'alignToTx',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'tPre',4,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(ip,'tPost',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip,'txField','Tx1',@(s)ischar(s) || isstring(s));
    addParameter(ip,'txSlice',4,@(x)isnumeric(x)&&isscalar(x)&&x>=1);

    % row scaling
    addParameter(ip,'rowscale','minmax',@(s) any(strcmpi(s,{'none','minmax','zscore'})));
    addParameter(ip,'clipZ',3,@(x)isnumeric(x)&&isscalar(x)&&x>0);

    parse(ip,varargin{:});
    P = ip.Results;

    dt_hours = dt/60;

    if ~isempty(P.maxCells)
        P.maxNon  = P.maxCells;
        P.maxMed  = P.maxCells;
        P.maxHigh = P.maxCells;
    end

    if isempty(P.smoothHrs)
        P.smoothHrs = P.Pmin/4;
    end

    classChan = string(P.classChan);

    % ---- defaults for "full-first then longest partial" filtering ----
    minPts = 5;
    minFiniteFrac = 0.90;   % require >= 90% finite points within the extracted window

    % ---- XY mask ----
    nXY = numel(dataloc.d);
    if isempty(P.xyInclude)
        xyMask = true(1,nXY);
    else
        if islogical(P.xyInclude)
            xyMask = P.xyInclude(:).';
            if numel(xyMask) ~= nXY
                error('xyInclude logical mask must match numel(dataloc.d).');
            end
        else
            xyMask = false(1,nXY);
            idx = P.xyInclude(:).';
            idx = idx(idx>=1 & idx<=nXY);
            xyMask(idx) = true;
        end
    end

    % ---- platemap mapping xy -> (row,col) for tx alignment ----
    hasPMD = isfield(dataloc,'platemapd') && isfield(dataloc.platemapd,'pmd') ...
          && isfield(dataloc.platemapd.pmd,'xy');
    if hasPMD
        pmd = dataloc.platemapd.pmd;
        xyGrid = pmd.xy;
    else
        pmd = [];
        xyGrid = [];
    end

    % ---- collect traces by class: FULL-window first, then LONGEST partial ----
    tr0_full = {}; tr1_full = {}; tr2_full = {};
    tr0_part = {}; tr1_part = {}; tr2_part = {};
    len0_full = []; len1_full = []; len2_full = [];
    len0_part = []; len1_part = []; len2_part = [];

    TmaxWanted = 0;   % desired window length (max over XYs of wantLen)
    TmaxSeen   = 0;   % max actual extracted length seen (for partial-only scenarios)

    for xy = 1:nXY
        if ~xyMask(xy), continue; end

        S = dataloc.d{xy};
        if isempty(S) || ~isstruct(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        if ~isfield(S,'osc') || ~isfield(S.osc, classChan)
            continue;
        end

        X = S.data.(chanName); % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end
        [nCellsXY, Txy] = size(X);

        % ---- class vector ----
        hasClass = isfield(S.osc.(classChan),'oscClass');
        hasBin   = isfield(S.osc.(classChan),'isOsc');

        if hasClass
            cls = uint8(S.osc.(classChan).oscClass(:));
        elseif P.useBinaryIfNoClass && hasBin
            isOsc = logical(S.osc.(chanName).isOsc(:));
            cls = zeros(size(isOsc), 'uint8');
            cls(isOsc) = uint8(P.binaryOscAs);
        else
            continue;
        end

        if numel(cls) ~= nCellsXY
            warning('XY %d: class label length mismatch; skipping.', xy);
            continue;
        end

        % ---- desired crop window for this XY ----
        if P.alignToTx
            TxFrame = 1; % fallback
            if hasPMD && isfield(pmd, char(P.txField))
                try
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
                end
            end
            if isempty(TxFrame) || ~isfinite(TxFrame) || TxFrame < 1
                TxFrame = 1;
            end

            t0_hr    = (TxFrame - 1) * dt_hours;
            iStart   = max(1, round((t0_hr - P.tPre)/dt_hours) + 1);
            iEndWant = round((t0_hr + P.tPost)/dt_hours) + 1;
        else
            iStart   = 1;
            iEndWant = Txy;
        end

        if iEndWant <= iStart
            continue;
        end

        wantLen = iEndWant - iStart + 1;
        TmaxWanted = max(TmaxWanted, wantLen);

        % If the XY can't even reach iStart, skip
        if Txy < iStart
            continue;
        end

        % ---- per-cell extraction with full-first logic ----
        isFullXY = (Txy >= iEndWant);
        iEndAvail = min(Txy, iEndWant);
        availLenXY = iEndAvail - iStart + 1;
        TmaxSeen = max(TmaxSeen, availLenXY);

        for i = 1:nCellsXY
            tr = X(i, iStart:iEndAvail);

            ok = isfinite(tr);
            if nnz(ok) < minPts
                continue;
            end
            if nnz(ok) < minFiniteFrac * numel(tr)
                continue;
            end

            switch cls(i)
                case 0
                    if isFullXY
                        tr0_full{end+1} = tr; len0_full(end+1) = numel(tr); %#ok<AGROW>
                    else
                        tr0_part{end+1} = tr; len0_part(end+1) = numel(tr); %#ok<AGROW>
                    end
                case 1
                    if isFullXY
                        tr1_full{end+1} = tr; len1_full(end+1) = numel(tr); %#ok<AGROW>
                    else
                        tr1_part{end+1} = tr; len1_part(end+1) = numel(tr); %#ok<AGROW>
                    end
                case 2
                    if isFullXY
                        tr2_full{end+1} = tr; len2_full(end+1) = numel(tr); %#ok<AGROW>
                    else
                        tr2_part{end+1} = tr; len2_part(end+1) = numel(tr); %#ok<AGROW>
                    end
            end
        end
    end

    % ---- choose per-class: full first, then longest partial to fill cap ----
    [tr0, len0] = local_take_pref_full(tr0_full, len0_full, tr0_part, len0_part, P.maxNon);
    [tr1, len1] = local_take_pref_full(tr1_full, len1_full, tr1_part, len1_part, P.maxMed);
    [tr2, len2] = local_take_pref_full(tr2_full, len2_full, tr2_part, len2_part, P.maxHigh);

    n0 = numel(tr0); n1 = numel(tr1); n2 = numel(tr2);
    if n0==0 && n1==0 && n2==0
        warning('No traces found for %s.', chanName);
        return;
    end

    % ---- decide final Tmax for display ----
    % If we got any full-window traces at all, use TmaxWanted.
    anyFull = ~isempty(tr0_full) || ~isempty(tr1_full) || ~isempty(tr2_full);
    if anyFull
        Tmax = TmaxWanted;
    else
        Tmax = max(1, TmaxSeen);
    end

    % ---- concatenate (class0 -> class1 -> class2) ----
    traces = [tr0(:); tr1(:); tr2(:)];
    clsAll = uint8([zeros(n0,1); ones(n1,1); 2*ones(n2,1)]);
    nCells = numel(traces);

    % boundaries for separators (before clustering)
    b01 = n0 + 0.5;
    b12 = n0 + n1 + 0.5;

    % ---- build RAW trace matrix [cells x Tmax] ----
    Mraw = NaN(nCells, Tmax);
    for i = 1:nCells
        v = traces{i};
        L = min(numel(v), Tmax);
        Mraw(i,1:L) = v(1:L);
    end

    % ---- ROW SCALING for display ----
    Mdisp = Mraw;
    switch lower(P.rowscale)
        case 'none'
            % keep raw
        case 'minmax'
            for i = 1:nCells
                row = Mraw(i,:);
                ok = isfinite(row);
                if nnz(ok) < 2
                    Mdisp(i,:) = NaN; continue;
                end
                mn = min(row(ok)); mx = max(row(ok));
                if mx <= mn
                    Mdisp(i,:) = NaN; continue;
                end
                row2 = row;
                row2(ok) = (row(ok) - mn) ./ (mx - mn);
                Mdisp(i,:) = row2;
            end
        case 'zscore'
            for i = 1:nCells
                row = Mraw(i,:);
                ok = isfinite(row);
                if nnz(ok) < 3
                    Mdisp(i,:) = NaN; continue;
                end
                mu = mean(row(ok));
                sd = std(row(ok));
                if sd == 0 || ~isfinite(sd)
                    Mdisp(i,:) = NaN; continue;
                end
                row2 = row;
                row2(ok) = (row(ok) - mu) ./ sd;
                row2(ok) = max(-P.clipZ, min(P.clipZ, row2(ok)));
                Mdisp(i,:) = row2;
            end
    end

    % ---- peak raster computed from RAW (not scaled) ----
    Mpk = zeros(nCells, Tmax);
    minDistSamples = max(1, round(P.Pmin / dt_hours * 0.5));
    smoothSamples  = max(1, round(P.smoothHrs / dt_hours));

    for i = 1:nCells
        v = Mraw(i,:);
        if all(~isfinite(v)), continue; end
        v2 = v; v2(isnan(v2)) = 0;

        if smoothSamples > 1
            vSmooth = movmean(v2, smoothSamples);
        else
            vSmooth = v2;
        end

        sigma = std(vSmooth);
        if sigma == 0 || ~isfinite(sigma), continue; end
        prom = P.promFactor * sigma;

        [~,locs] = findpeaks(vSmooth, ...
            'MinPeakDistance', minDistSamples, ...
            'MinPeakProminence', prom);

        locs = locs(locs>=1 & locs<=Tmax);
        Mpk(i,locs) = 1;
    end

    % ---- cluster WITHIN each class using Mdisp (scaled) ----
    idx0 = find(clsAll==0);
    idx1 = find(clsAll==1);
    idx2 = find(clsAll==2);

    ord0 = idx0(cluster_order_rows(Mdisp(idx0,:)));
    ord1 = idx1(cluster_order_rows(Mdisp(idx1,:)));
    ord2 = idx2(cluster_order_rows(Mdisp(idx2,:)));

    order = [ord0(:); ord1(:); ord2(:)];

    % apply ordering everywhere
    Mraw   = Mraw(order,:);
    Mdisp  = Mdisp(order,:);
    Mpk    = Mpk(order,:);
    clsAll = clsAll(order);

    % boundaries remain counts-based (same class sizes)
    n0 = numel(ord0);
    n1 = numel(ord1);
    n2 = numel(ord2);
    b01 = n0 + 0.5;
    b12 = n0 + n1 + 0.5;

    % ---- time axis ----
    if P.alignToTx
        % Use requested window if any full exists; otherwise use actual Tmax
        if anyFull
            th = linspace(-P.tPre, P.tPost, Tmax);
        else
            th = linspace(-P.tPre, -P.tPre + (Tmax-1)*dt_hours, Tmax);
        end
        xlab = sprintf('time (h) relative to %s', char(P.txField));
    else
        th = (0:Tmax-1)*dt_hours;
        xlab = 'time (h)';
    end

    % ---- plot ----
    figure('Name',sprintf('Osc matrix (3-class): %s', chanName), ...
           'NumberTitle','off','Color','w');

    subplot(2,1,1);
    imagesc(th, 1:nCells, Mdisp);
    set(gca,'YDir','reverse');
    colormap(gca, parula);
    colorbar;
    ylabel('cells (0 bottom \rightarrow 2 top)');
    xlim([th(1) th(end)]);
    hold on;
    if n0>0 && (n1>0 || n2>0), plot([th(1) th(end)], [b01 b01], 'k-', 'LineWidth', 1.2); end
    if n1>0 && n2>0,           plot([th(1) th(end)], [b12 b12], 'k-', 'LineWidth', 1.2); end
    if P.alignToTx, xline(0,'k-','LineWidth',1); end
    hold off; box on;
    title(sprintf('%s (0=non, 1=weak, 2=strong)  rowscale=%s', chanName, lower(P.rowscale)));

    subplot(2,1,2);
    imagesc(th, 1:nCells, Mpk);
    set(gca,'YDir','normal');
    colormap(gca, gray);
    colorbar;
    xlabel(xlab);
    ylabel('cells (same order)');
    title('Peak raster');
    xlim([th(1) th(end)]);
    hold on;
    if n0>0 && (n1>0 || n2>0), plot([th(1) th(end)], [b01 b01], 'r-', 'LineWidth', 1.2); end
    if n1>0 && n2>0,           plot([th(1) th(end)], [b12 b12], 'r-', 'LineWidth', 1.2); end
    if P.alignToTx, xline(0,'r-','LineWidth',1); end
    hold off; box on;
end


function [trOut, lenOut] = local_take_pref_full(trF, lenF, trP, lenP, maxN)
% Keep up to maxN. Prefer full-window traces first. If not enough, fill with
% longest partial traces.

    trOut = {};
    lenOut = [];

    % start with full traces
    if ~isempty(trF)
        trOut = trF(:);
        lenOut = lenF(:);
    end

    % if too many full traces, trim (they're generally equal length, but be safe)
    if numel(trOut) > maxN
        [~, ord] = sort(lenOut, 'descend');
        ord = ord(1:maxN);
        trOut = trOut(ord);
        lenOut = lenOut(ord);
        return;
    end

    % fill remaining with longest partial
    need = maxN - numel(trOut);
    if need <= 0 || isempty(trP)
        return;
    end

    [~, ordP] = sort(lenP(:), 'descend');
    ordP = ordP(1:min(need, numel(ordP)));

    trOut  = [trOut;  trP(ordP(:))];
    lenOut = [lenOut; lenP(ordP(:))];
end


function ordLocal = cluster_order_rows(Min)
% Min is [n x T], should be row-scaled already.
% Returns permutation indices 1..n that makes similar rows adjacent.

    n = size(Min,1);
    if n <= 2
        ordLocal = 1:n;
        return;
    end

    X = Min;
    X(~isfinite(X)) = 0;

    D = pdist(X, 'correlation');
    Z = linkage(D, 'average');

    try
        ordLocal = optimalleaforder(Z, D);
    catch
        % fallback if optimalleaforder not available
        [~, ordLocal] = sort(X(:,1));
    end
end
