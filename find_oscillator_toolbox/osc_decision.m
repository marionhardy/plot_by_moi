function isOsc = osc_decision(feat, params)
% osc_decision  Apply simple threshold logic to osc_features.
%
% isOsc = osc_decision(feat, params)
%
% feat   : struct from osc_features
% params : struct with fields:
%           .bandFracMin
%           .cvIPIMax
%           .ampMin
%
% Output:
%   isOsc : logical true/false

    isOsc = (feat.bandFrac >= params.bandFracMin) && ...
            (feat.cvIPI    <= params.cvIPIMax)    && ...
            (feat.amp      >= params.ampMin);
end
