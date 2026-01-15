function [bandFrac, cvIPI, amp, y] = osc_training_features(dataloc, osc_labels, chanName, dt, Pmin, Pmax)
% osc_training_features  Extract features for labeled cells (for QC plots).
%
% [bandFrac, cvIPI, amp, y] = osc_training_features(...)
%
% Inputs:
%   dataloc    : main dataloc struct
%   osc_labels : struct from label_oscillators_all
%   chanName   : channel used for labeling
%   dt, Pmin, Pmax : sampling and period-band parameters
%
% Outputs:
%   bandFrac, cvIPI, amp : feature vectors
%   y                    : labels (0/1)

    n = numel(osc_labels.cells);
    bandFrac = nan(n,1);
    cvIPI    = nan(n,1);
    amp      = nan(n,1);
    y        = osc_labels.labels;

    for k = 1:n
        if ~(y(k)==0 || y(k)==1), continue; end
        c = osc_labels.cells(k);
        S = dataloc.d{c.xy};
        if isempty(S) || ~isfield(S,'data') || ~isfield(S.data, chanName)
            continue;
        end
        tr = S.data.(chanName)(c.row,:);
        feat = osc_features(tr, dt, Pmin, Pmax, 3);
        bandFrac(k) = feat.bandFrac;
        cvIPI(k)    = feat.cvIPI;
        amp(k)      = feat.amp;
    end

    mask = (y==0 | y==1) & ~isnan(bandFrac) & ~isnan(cvIPI) & ~isnan(amp);
    bandFrac = bandFrac(mask);
    cvIPI    = cvIPI(mask);
    amp      = amp(mask);
    y        = y(mask);
end
