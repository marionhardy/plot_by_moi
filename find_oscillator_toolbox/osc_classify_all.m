function dataloc2 = osc_classify_all(dataloc, chanName, dt, Pmin, Pmax, params, minTrackPts)
% osc_classify_all  Classify oscillations for every cell, every XY, all conditions.
%
% Adds fields:
%   d{xy}.osc.(chanName).isOsc     (binary)
%   d{xy}.osc.(chanName).oscClass  (0=non, 1=medium, 2=high)   % NEW
%   d{xy}.osc.(chanName).oscScore  (0..1)                      % NEW
%   d{xy}.osc.(chanName).bandFrac
%   d{xy}.osc.(chanName).cvIPI
%   d{xy}.osc.(chanName).amp

    dataloc2 = dataloc;

    for xy = 1:numel(dataloc2.d)
        S = dataloc2.d{xy};
        if isempty(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end

        X = S.data.(chanName);   % [nCells x T]
        if isempty(X) || ~isnumeric(X)
            continue;
        end
        nCells = size(X,1);

        isOsc     = false(nCells,1);
        oscClass  = zeros(nCells,1,'uint8');  % NEW
        oscScore  = nan(nCells,1);            % NEW

        bandFrac  = nan(nCells,1);
        cvIPI     = nan(nCells,1);
        amp       = nan(nCells,1);

        for i = 1:nCells
            tr = X(i,:);
            effLen = sum(isfinite(tr));

            if effLen < minTrackPts
                % Too short to judge
                bandFrac(i) = NaN;
                cvIPI(i)    = NaN;
                amp(i)      = NaN;
                isOsc(i)    = false;

                oscClass(i) = uint8(0);  % NEW
                oscScore(i) = NaN;       % NEW
                continue;
            end

            feat = osc_features(tr, dt, Pmin, Pmax, 3);
            bandFrac(i) = feat.bandFrac;
            cvIPI(i)    = feat.cvIPI;
            amp(i)      = feat.amp;

            % Keep binary EXACTLY as before:
            isOsc(i) = osc_decision(feat, params);

            % Add graded outputs (NEW):
            [oscClass(i), oscScore(i)] = osc_decision_3class(feat, params);
        end

        if ~isfield(S,'osc') || ~isstruct(S.osc)
            S.osc = struct();
        end

        S.osc.(chanName).isOsc     = isOsc;
        S.osc.(chanName).oscClass  = oscClass;  % NEW
        S.osc.(chanName).oscScore  = oscScore;  % NEW

        S.osc.(chanName).bandFrac  = bandFrac;
        S.osc.(chanName).cvIPI     = cvIPI;
        S.osc.(chanName).amp       = amp;

        dataloc2.d{xy} = S;
    end
end
