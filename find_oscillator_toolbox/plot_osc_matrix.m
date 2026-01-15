function plot_osc_matrix(dataloc, chanName, dt, varargin)
% plot_osc_matrix  Bulk view of traces by oscillation class (0/1/2),
% with Tx-alignment + cropping, plus peak raster.
%
% NEW: per-row scaling (min-max or zscore) so each cell is visible.
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

    % ---- collect traces by class (cropped if requested) ----
    tr0 = {}; tr1 = {}; tr2 = {};
    len0 = []; len1 = []; len2 = [];
    Tmax = 0;

    for xy = 1:nXY
        if ~xyMask(xy), continue; end

        S = dataloc.d{xy};
        if isempty(S) || ~isstruct(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        if ~isfield(S,'osc') || ~isfield(S.osc, chanName)
            continue;
        end

        X = S.data.(chanName); % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end
        [nCellsXY, Txy] = size(X);

        % ---- class vector ----
        hasClass = isfield(S.osc.(chanName),'oscClass');
        hasBin   = isfield(S.osc.(chanName),'isOsc');

        if hasClass
            cls = uint8(S.osc.(chanName).oscClass(:));
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

        % ---- crop window ----
        iStart = 1; iEnd = Txy;

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

            t0_hr  = (TxFrame - 1) * dt_hours;
            iStart = max(1, round((t0_hr - P.tPre)/dt_hours) + 1);
            iEnd   = min(Txy, round((t0_hr + P.tPost)/dt_hours) + 1);

            if iEnd <= iStart
                continue;
            end
        end

        cropLen = iEnd - iStart + 1;
        Tmax = max(Tmax, cropLen);

        for i = 1:nCellsXY
            tr = X(i, iStart:iEnd);
            if all(~isfinite(tr)), continue; end
            effLen = sum(isfinite(tr));
            if effLen < 5, continue; end

            switch cls(i)
                case 0
                    tr0{end+1} = tr;      %#ok<AGROW>
                    len0(end+1) = effLen; %#ok<AGROW>
                case 1
                    tr1{end+1} = tr;      %#ok<AGROW>
                    len1(end+1) = effLen; %#ok<AGROW>
                case 2
                    tr2{end+1} = tr;      %#ok<AGROW>
                    len2(end+1) = effLen; %#ok<AGROW>
            end
        end
    end

    n0 = numel(tr0); n1 = numel(tr1); n2 = numel(tr2);
    if n0==0 && n1==0 && n2==0
        warning('No traces found for %s.', chanName);
        return;
    end

    % ---- subsample each class by longest ----
    if n0 > P.maxNon
        [~, ord] = sort(len0, 'descend');
        keep = ord(1:P.maxNon);
        tr0 = tr0(keep); len0 = len0(keep); n0 = numel(tr0);
    end
    if n1 > P.maxMed
        [~, ord] = sort(len1, 'descend');
        keep = ord(1:P.maxMed);
        tr1 = tr1(keep); len1 = len1(keep); n1 = numel(tr1);
    end
    if n2 > P.maxHigh
        [~, ord] = sort(len2, 'descend');
        keep = ord(1:P.maxHigh);
        tr2 = tr2(keep); len2 = len2(keep); n2 = numel(tr2);
    end

    % ---- concatenate (class0 -> class1 -> class2) ----
    traces = [tr0(:); tr1(:); tr2(:)];
    clsAll = uint8([zeros(n0,1); ones(n1,1); 2*ones(n2,1)]);   % ✅ FIX: define clsAll here
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
    Mraw  = Mraw(order,:);
    Mdisp = Mdisp(order,:);
    Mpk   = Mpk(order,:);
    clsAll = clsAll(order);

    % boundaries remain counts-based (same class sizes)
    n0 = numel(ord0);
    n1 = numel(ord1);
    n2 = numel(ord2);
    b01 = n0 + 0.5;
    b12 = n0 + n1 + 0.5;

    % ---- time axis ----
    if P.alignToTx
        th = linspace(-P.tPre, P.tPost, Tmax);
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
    set(gca,'YDir','normal');
    colormap(gca, parula);
    colorbar;
    ylabel('cells (0 bottom → 2 top)');
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
