function params = osc_fit_thresholds(dataloc, osc_labels, chanName, dt, Pmin, Pmax, minTrackPts)
% osc_fit_thresholds  Fit thresholds for oscillation detection using manual labels.
%
% Uses a 3D grid search over bandFracMin, ampMin, and cvIPIMax
% optimizing the F1-score for the oscillating class (y==1).

    n = numel(osc_labels.cells);
    bandFrac = nan(n,1);
    cvIPI    = nan(n,1);
    amp      = nan(n,1);
    y        = osc_labels.labels;  % 1=osc, 0=not, -1=unsure

    % ---- compute features for all labeled cells ----
    for k = 1:n
        if ~(y(k)==0 || y(k)==1), continue; end
        c = osc_labels.cells(k);
        S = dataloc.d{c.xy};
        if isempty(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        tr = S.data.(chanName)(c.row,:);
        effLen = sum(isfinite(tr));
        if effLen < minTrackPts
            continue;  % don’t use this example to fit thresholds
        end
        feat = osc_features(tr, dt, Pmin, Pmax, 3);
        bandFrac(k) = feat.bandFrac;
        cvIPI(k)    = feat.cvIPI;
        amp(k)      = feat.amp;
    end

    % ---- keep only valid, labeled examples ----
    mask = (y==0 | y==1) & ~isnan(bandFrac) & ~isnan(amp) & ~isnan(cvIPI);
    bandFrac = bandFrac(mask);
    cvIPI    = cvIPI(mask);
    amp      = amp(mask);
    y        = y(mask);

    if isempty(y)
        error('No valid labeled examples to fit thresholds.');
    end

    % ---- define search grids for thresholds ----
    bfGrid = linspace(min(bandFrac), max(bandFrac), 30);
    aGrid  = linspace(min(amp),      max(amp),      30);
    % cvIPI: lower tends to be more regular; we'll search an *upper* bound
    cvGrid = linspace(min(cvIPI), max(cvIPI), 20);

    bestScore = -inf;
    bestBF    = bfGrid(1);
    bestAmp   = aGrid(1);
    bestCv    = cvGrid(end);
    bestStats = struct('tp',0,'fp',0,'tn',0,'fn',0);

    % ---- grid search to maximize F1 for oscillators ----
    for ib = 1:numel(bfGrid)
        bfThr = bfGrid(ib);
        for ia = 1:numel(aGrid)
            aThr = aGrid(ia);
            for ic = 1:numel(cvGrid)
                cvThr = cvGrid(ic);

                yhat = (bandFrac >= bfThr) & (amp >= aThr) & (cvIPI <= cvThr);

                tp = sum( yhat==1 & y==1 );
                tn = sum( yhat==0 & y==0 );
                fp = sum( yhat==1 & y==0 );
                fn = sum( yhat==0 & y==1 );

                if tp+fn == 0 || tp+fp == 0
                    continue; % no positives or no predicted positives
                end

                sens = tp / (tp+fn);              % recall for osc
                spec = tn / (tn+fp);              % specificity
                prec = tp / (tp+fp);              % precision for osc
                F1   = 2*prec*sens/(prec+sens);   % F1 for osc

                % Optional: enforce some minimum specificity
                if spec < 0.5
                    continue;
                end

                score = F1;

                if score > bestScore
                    bestScore = score;
                    bestBF    = bfThr;
                    bestAmp   = aThr;
                    bestCv    = cvThr;
                    bestStats = struct('tp',tp,'fp',fp,'tn',tn,'fn',fn, ...
                                       'sens',sens,'spec',spec,'prec',prec,'F1',F1);
                end
            end
        end
    end

    params.bandFracMin = bestBF;
    params.ampMin      = bestAmp;
    params.cvIPIMax    = bestCv;

    fprintf('Fitted thresholds (3D grid, F1 for osc):\n');
    fprintf('  bandFracMin = %.4f\n', params.bandFracMin);
    fprintf('  ampMin      = %.4f\n', params.ampMin);
    fprintf('  cvIPIMax    = %.4f\n', params.cvIPIMax);
    fprintf('  tp=%d, fp=%d, tn=%d, fn=%d\n', ...
        bestStats.tp, bestStats.fp, bestStats.tn, bestStats.fn);
    fprintf('  sens=%.2f, spec=%.2f, prec=%.2f, F1=%.2f\n', ...
        bestStats.sens, bestStats.spec, bestStats.prec, bestStats.F1);
end
