function [keep, info] = dl_cycif_regmask(dataloc, xyIdx, varargin)
% DL_CYCIF_REGMASK  Per-round registration QC gate for CycIF endpoint data.
%
% [keep, info] = dl_cycif_regmask(dataloc, xyIdx, ...)
%
% INPUT  (paired values after the two positional args)
%   dataloc       : dataloc struct (v5). Reads dataloc.d{xyIdx}.cycif_reg.
%   xyIdx         : [scalar] xy index to evaluate.
%   'requireReg'  : [logical, default true] if true, rounds with
%                   cycif_reg(r).failed == true are gated OUT.
%   'minAlignFrac': [scalar, default []] if set, additionally gate OUT any
%                   round whose cycif_reg(r).align_frac < this threshold.
%
% OUTPUT
%   keep : [1 x R] logical, true for rounds that PASS QC (usable).
%          R = numel(dataloc.d{xyIdx}.cycif_reg). If cycif_reg is absent,
%          keep = true(1,0) and info.hasReg = false (nothing to gate).
%   info : struct with fields
%          .hasReg     logical, whether cycif_reg existed
%          .nRounds    double,  R
%          .failedIdx  double,  rounds gated for .failed
%          .lowAlign   double,  rounds gated for align_frac threshold
%          .alignFrac  [1 x R], the align_frac per round (NaN if absent)
%
% -------------------------------------------------------------------------
% V5 CHANGE (crit #7): v4's plot_by_MHY had no concept of CycIF
% registration QC -- an endpoint panel would plot ANTIBODY_nuc/_cyt/_cell
% intensities from a round the pipeline had already flagged
% cycif_reg(r).failed == true. The v5 dataloc carries per-round QC
% (failed, align_frac, median_dist_px, ...); this helper turns it into a
% keep-mask so endpoint plots can honor it via a 'requireReg' param.
% Report-only + returns a mask; it does NOT mutate dataloc.
%
% BIOLOGICAL RATIONALE: a failed registration means live-cell tracks and
% fixed-cell nuclei were not reliably matched, so any endpoint intensity
% assigned to a track in that round is effectively a random draw. Plotting
% it as if it were that cell's marker level is worse than dropping it.
% -------------------------------------------------------------------------

ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'requireReg', true, @islogical);
addParameter(ip, 'minAlignFrac', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(ip, varargin{:});
requireReg   = ip.Results.requireReg;
minAlignFrac = ip.Results.minAlignFrac;

info = struct('hasReg', false, 'nRounds', 0, 'failedIdx', [], ...
              'lowAlign', [], 'alignFrac', []);

if xyIdx > numel(dataloc.d) || isempty(dataloc.d{xyIdx}) ...
        || ~isfield(dataloc.d{xyIdx}, 'cycif_reg') || isempty(dataloc.d{xyIdx}.cycif_reg)
    keep = true(1,0);
    return
end

reg = dataloc.d{xyIdx}.cycif_reg;
R = numel(reg);
info.hasReg  = true;
info.nRounds = R;
keep = true(1, R);

% align_frac per round (NaN where field absent)
af = nan(1, R);
for r = 1:R
    if isfield(reg(r), 'align_frac') && ~isempty(reg(r).align_frac)
        af(r) = reg(r).align_frac;
    end
end
info.alignFrac = af;

% Gate on explicit failure flag
if requireReg
    for r = 1:R
        if isfield(reg(r), 'failed') && ~isempty(reg(r).failed) && reg(r).failed
            keep(r) = false;
            info.failedIdx(end+1) = r; %#ok<AGROW>
        end
    end
end

% Gate on align_frac threshold
if ~isempty(minAlignFrac)
    low = af < minAlignFrac;   % NaN < x is false, so absent align_frac is not gated here
    info.lowAlign = find(low);
    keep(low) = false;
end
end
