function [oscClass, oscScore, isOsc] = osc_decision_3class(feat, params, varargin)
% osc_decision_3class
% Uses one-sided thresholds in params:
%   params.bandFracMin (higher better)
%   params.ampMin      (higher better)
%   params.cvIPIMax    (lower better)
%
% Outputs:
%   oscClass: 0=non, 1=medium, 2=high
%   oscScore: 0..1 (confidence-ish)
%   isOsc   : binary (oscClass>=1)

p.ScoreCuts   = [0.45 0.75];   % [mediumCut, highCut]
p.Margins     = struct('bandFrac',0.25,'amp',0.25,'cvIPI',0.25); % fractional margins
p.Weights     = [0.45 0.35 0.20]; % bandFrac, amp, cvIPI

nin = numel(varargin);
if rem(nin,2)~=0, error('Name/value pairs expected.'); end
for k=1:2:nin
    p.(varargin{k}) = varargin{k+1};
end

% --- thresholds ---
bMin = params.bandFracMin;
aMin = params.ampMin;
cMax = params.cvIPIMax;

% Define "medium band" around the thresholds using margins.
% For bandFrac/amp: low->high increases goodness
bLo = bMin*(1 - p.Margins.bandFrac);
bHi = bMin*(1 + p.Margins.bandFrac);

% amp can be near 0 or negative; use additive margin if tiny
if abs(aMin) < 0.05
    aLo = aMin - 0.05;
    aHi = aMin + 0.05;
else
    aLo = aMin*(1 - p.Margins.amp);
    aHi = aMin*(1 + p.Margins.amp);
end

% For cvIPI: lower is better. Define a "good zone" below cMax.
cHi = cMax*(1 + p.Margins.cvIPI); % worse boundary
cLo = cMax*(1 - p.Margins.cvIPI); % better boundary

% --- per-feature scores in [0,1] ---
sBand = clamp01((feat.bandFrac - bLo) / max(eps, (bHi - bLo)));
sAmp  = clamp01((feat.amp      - aLo) / max(eps, (aHi - aLo)));
sCv   = clamp01((cHi - feat.cvIPI)    / max(eps, (cHi - cLo)));  % 1 if low cvIPI

w = p.Weights(:); w = w/sum(w);
oscScore = w(1)*sBand + w(2)*sAmp + w(3)*sCv;

% --- 3-class mapping ---
mCut = p.ScoreCuts(1);
hCut = p.ScoreCuts(2);

if isnan(oscScore)
    oscClass = uint8(0);
elseif oscScore >= hCut
    oscClass = uint8(2);
elseif oscScore >= mCut
    oscClass = uint8(1);
else
    oscClass = uint8(0);
end

isOsc = (oscClass >= 1);

end

function y = clamp01(x)
y = min(1, max(0, x));
end
