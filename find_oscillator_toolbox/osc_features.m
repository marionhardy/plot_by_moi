function feat = osc_features(trace, dt, Pmin, Pmax, smoothWin)
% osc_features  Extract oscillation-related features from one time series.
%
% feat = osc_features(trace, dt, Pmin, Pmax, smoothWin)
%
% Inputs:
%   trace     : 1 x T or T x 1 vector (numeric)
%   dt        : sampling interval (same units as Pmin/Pmax)
%   Pmin      : minimum period of interest
%   Pmax      : maximum period of interest
%   smoothWin : (optional) moving average window (in samples, default 3)
%
% Output:
%   feat struct with fields:
%     .bandFrac : max power in band / total power
%     .cvIPI    : coefficient of variation of inter-peak intervals
%     .amp      : typical oscillation amplitude

    if nargin < 5 || isempty(smoothWin)
        smoothWin = 3;
    end

    x = trace(:)';  % force row
    if any(isnan(x))
        x = fillmissing(x,'linear');
    end
    x = detrend(x);

    if smoothWin > 1 && numel(x) >= smoothWin
        x_smooth = movmean(x, smoothWin);
    else
        x_smooth = x;
    end

    % frequency band of interest
    fs = 1/dt;
    fmin = 1/Pmax;
    fmax = 1/Pmin;

    % power spectrum
    if all(x_smooth == x_smooth(1))
        % flat signal -> avoid numerical weirdness
        Pxx = 0;
        f   = 0;
    else
        [Pxx,f] = periodogram(x_smooth,[],[],fs);
    end

    bandIdx = f >= fmin & f <= fmax;
    if any(bandIdx)
        bandPower = max(Pxx(bandIdx));
    else
        bandPower = 0;
    end
    totalPower = sum(Pxx);
    if totalPower > 0
        bandFrac = bandPower / totalPower;
    else
        bandFrac = 0;
    end

    % peaks & IPIs
    [pks,locs] = findpeaks(x_smooth);
    if numel(locs) >= 3
        IPIs = diff(locs) * dt;
        if mean(IPIs) > 0
            cvIPI = std(IPIs) / mean(IPIs);
        else
            cvIPI = Inf;
        end
        amp   = median(pks) - median(x_smooth);
    else
        cvIPI = Inf;
        amp   = 0;
    end

    feat.bandFrac = bandFrac;
    feat.cvIPI    = cvIPI;
    feat.amp      = amp;
end
