function fig = cycif_feature_hist(T, varargin)
% cycif_feature_hist  Per-marker distribution diagnostic for setting the
%                     embedding/clustering transform and arcsinh cofactor.
%
%   Plots each marker's intensity distribution on the candidate scales so the
%   transform + cofactor are chosen from the DATA (noise floor, bimodality),
%   not a default. Run ONCE per dataset before committing prep parameters.
%
%   Read it like this:
%     - LINEAR panel: where does the background floor sit? Is there a long
%       right tail (there always is)? Is the marker bimodal (a real ON mode
%       separated from background = the signal you want clustering to see)?
%     - ARCSINH panel: a good cofactor puts the arcsinh KNEE at the noise floor,
%       so background collapses to ~0 and the ON mode spreads out. If the whole
%       distribution looks linear (no compression), cofactor is too high; if the
%       background is stretched into structure, it is too low.
%     - The printed per-marker summary gives robust background/foreground stats
%       to pick cofactor numerically.
%
% Usage examples:
%   cycif_feature_hist(T);                              % arcsinh sweep, all markers
%   cycif_feature_hist(T, 'cofactors',[50 150 500]);    % compare cofactors
%   cycif_feature_hist(T, 'markers',{'iNOS_cyt','PSTAT6_nuc'});
%   fig = cycif_feature_hist(T, 'save_path','.../hist');
%
% INPUT
%   T -- cell table from cycif_build_table
%
% NAME-VALUE OPTIONS
%   'markers'    {1xM}    markers to show; [] = auto-detect       (default [])
%   'cofactors'  vector   arcsinh cofactors to overlay/compare    (default [150])
%   'nbins'      scalar   histogram bins                          (default 60)
%   'max_pts'    scalar   subsample cells for speed (Inf=all)     (default 5e4)
%   'save_path'  string   SVG path prefix; '' = no save           (default '')
%
% OUTPUT
%   fig -- figure handle (grid of markers x scales)
%
% Also prints, per marker, robust background vs foreground estimates and a
% suggested cofactor (= background level, the natural knee).
%
% See also: cycif_prep_features.

p = inputParser();
addParameter(p, 'markers',   []);
addParameter(p, 'cofactors', 150);
addParameter(p, 'nbins',     60);
addParameter(p, 'max_pts',   5e4);
addParameter(p, 'save_path', '');
parse(p, varargin{:});
markers   = p.Results.markers;
cofs      = p.Results.cofactors(:)';
nbins     = p.Results.nbins;
max_pts   = p.Results.max_pts;
save_path = p.Results.save_path;

if isempty(markers)
    cols   = T.Properties.VariableNames;
    is_num = cellfun(@(c) isnumeric(T.(c)) && ~islogical(T.(c)), cols);
    is_id  = ismember(cols, {'xy','cellid','uid','ptx_is_spike'});
    markers = cols(is_num & ~is_id);
end
markers = cellstr(markers);
nMk = numel(markers);

X = table2array(T(:, markers));
if size(X,1) > max_pts
    X = X(randperm(size(X,1), max_pts), :);
end

% --- per-marker numeric summary + suggested cofactor ---------------------
% Background estimated as the mode-side (robust low quantile among positive
% values); foreground as a high quantile. Suggested cofactor ~ background so the
% arcsinh knee sits at the noise floor.
fprintf('\n[cycif_feature_hist] Per-marker summary (linear intensity)\n');
fprintf('%-16s %10s %10s %10s %10s %12s\n', ...
    'marker','bg(p50+)','p90','p99','max','suggest_cof');
sugg = nan(1, nMk);
for m = 1:nMk
    v = X(:, m); v = v(~isnan(v));
    if isempty(v), continue; end
    vpos = v(v > 0);
    if isempty(vpos), vpos = v; end
    bg   = median(vpos);            % background-dominated (most cells are off)
    p90  = quantile(v, 0.90);
    p99  = quantile(v, 0.99);
    sugg(m) = max(bg, eps);
    fprintf('%-16s %10.4g %10.4g %10.4g %10.4g %12.4g\n', ...
        markers{m}, bg, p90, p99, max(v), sugg(m));
end
fprintf(['\n  Suggested cofactor (per marker = background median). A single shared\n' ...
         '  cofactor near the MEDIAN of these (%.4g) is a reasonable global start;\n' ...
         '  per-marker cofactors are possible if backgrounds differ wildly.\n\n'], ...
         median(sugg,'omitnan'));

% --- figure: markers (rows) x [linear | log2 | arcsinh-per-cofactor] -----
nScale = 2 + numel(cofs);   % linear, log2, arcsinh(each cofactor)
scaleLabels = [{'linear','log2'}, arrayfun(@(c) sprintf('asinh/%.4g',c), cofs, 'uni', 0)];

nShow = min(nMk, 12);   % cap rows for legibility
if nShow < nMk
    fprintf('[cycif_feature_hist] Showing first %d of %d markers; pass ''markers'' to select.\n', nShow, nMk);
end

fig = figure('Color','w', 'Position',[80 80 260*nScale 150*nShow]);
tl  = tiledlayout(fig, nShow, nScale, 'TileSpacing','compact', 'Padding','compact');
for r = 1:nShow
    v = X(:, r); v = v(~isnan(v));
    for c = 1:nScale
        nexttile;
        switch c
            case 1, d = v;                       % linear
            case 2, d = log2(max(v,1));          % log2
            otherwise
                cof = cofs(c-2);
                d = asinh(v ./ cof);             % arcsinh
        end
        histogram(d, nbins, 'EdgeColor','none', 'FaceColor',[0.2 0.45 0.75]);
        if r == 1, title(scaleLabels{c}, 'Interpreter','none', 'FontSize',9); end
        if c == 1, ylabel(markers{r}, 'Interpreter','none', 'FontSize',9, 'FontWeight','bold'); end
        set(gca,'YTick',[],'FontSize',7); box off;
    end
end
title(tl, 'Marker distributions by transform  (pick cofactor: knee at background floor)', ...
    'FontSize', 11);

if ~isempty(save_path)
    fpath = [char(save_path) '_feature_hist.svg'];
    fig.Units='inches'; sz=fig.Position(3:4)/96;
    fig.PaperUnits='inches'; fig.PaperSize=sz; fig.PaperPosition=[0 0 sz]; fig.PaperPositionMode='manual';
    print(fig, fpath, '-dsvg', '-painters');
    fprintf('Saved SVG: %s\n', fpath);
end
end
