function Q = cycif_qc_mask_signal(dataloc, varargin)
% CYCIF_QC_MASK_SIGNAL
% Raw-data-led QC for endpoint segmentation. The DISTRIBUTION is the summary
% and every histogram bin is backed by the actual cell crops that populate
% it: click a bin to pop a montage of those cells (raw | raw+outline | label).
% No console fractions are relied on -- the plots and images carry the QC.
%
% Interrogates two things per cell, both from the saved masks + the raw
% segmentation channel (resolved the SAME way ct_save_endpoint_masks did):
%
%   1. SIGNAL PRESENCE  -- snr = (median in-mask) - (median local ring), in
%      raw units. Answers "is there really a nucleus under this mask?".
%      Near-zero/negative snr = spurious or empty detection.
%
%   2. BOUNDARY CORRECTNESS -- bnd = gradient magnitude sampled ON the mask
%      perimeter, normalized by the in-mask signal range. A mask whose edge
%      sits on a real intensity edge scores high; a mask edge drawn across
%      flat interior or empty space scores low (over/under-segmentation, or
%      boundary drawn where there is no object edge).
%
% Both axes get a linked histogram; clicking any bar shows the cells there.
%
% INPUTS
%   dataloc   struct   Populated dataloc.
%
% PAIRED OPTIONS
%   'op'         struct  REQUIRED. celltracer op (channel resolution + donut).
%   'ip'         struct  REQUIRED. celltracer image params (needs .fname).
%   'xy_list'    [1xN]   XYs to audit. Empty -> all with a mask file. ([])
%   'mask_dir'   char    (default dataloc.exp.proc_dir)
%   'mask_name'  char    (default dataloc.exp.name)
%   'max_cells_per_xy' scalar  Subsample cap per XY for speed.   (default inf)
%   'thumbs_per_bin'   scalar  Max crops shown when a bin is clicked. (25)
%   'save_dir'   char    If set, save the two summary figures as PNG.  ('')
%
% OUTPUT
%   Q  struct with:
%     .xy, .lab            [nCell x 1] source XY and label id per cell
%     .snr                 [nCell x 1] signal-presence metric
%     .bnd                 [nCell x 1] boundary-agreement metric
%     .area                [nCell x 1] nucleus area (px^2)
%     .centroid            [nCell x 2] (x,y)
%     .meta                paths + resolved channel, for the click callback
%
% USAGE
%   Q = cycif_qc_mask_signal(dataloc, 'op', op, 'ip', ip);
%   Q = cycif_qc_mask_signal(dataloc, 'op', op, 'ip', ip, 'xy_list', [12 45]);

p = inputParser;
p.addParameter('op', [], @(x) isstruct(x));
p.addParameter('ip', [], @(x) isstruct(x));
p.addParameter('xy_list', [], @(x) isnumeric(x));
p.addParameter('mask_dir', '', @(x) ischar(x) || isstring(x));
p.addParameter('mask_name', '', @(x) ischar(x) || isstring(x));
p.addParameter('max_cells_per_xy', inf, @(x) isnumeric(x) && isscalar(x));
p.addParameter('thumbs_per_bin', 25, @(x) isnumeric(x) && isscalar(x));
p.addParameter('save_dir', '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
op = p.Results.op; ip = p.Results.ip;
assert(~isempty(op) && ~isempty(ip) && isfield(ip,'fname'), ...
    'cycif_qc_mask_signal: op and ip (with .fname) are required.');

mask_dir  = char(p.Results.mask_dir);
mask_name = char(p.Results.mask_name);
if isempty(mask_dir); mask_dir = dataloc.exp.proc_dir; end
if isempty(mask_name); mask_name = dataloc.exp.name; end

nuc_chan = local_resolve_chan(op);
assert(~isnan(nuc_chan), ...
    'cycif_qc_mask_signal: could not resolve segmentation channel from op.');

% XY list
xy_list = p.Results.xy_list;
if isempty(xy_list)
    cand = find(cellfun(@(d) isstruct(d), dataloc.d));
    keep = arrayfun(@(x) exist(local_maskfile(mask_dir,mask_name,x),'file')==2, cand);
    xy_list = cand(keep);
end
xy_list = xy_list(:)';

dao = ImageAccess(ip.fname, ip.indsz);

% ---- Per-cell metric extraction ---------------------------------------
all_xy=[]; all_lab=[]; all_snr=[]; all_bnd=[]; all_area=[]; all_cen=[];
for xy = xy_list
    fmask = local_maskfile(mask_dir, mask_name, xy);
    Mf = load(fmask, 'nuc_labels', 'last_t');
    nuc = Mf.nuc_labels;
    raw = double(dao.GetFrame([Mf.last_t, nuc_chan, xy, 1]));

    labs = unique(nuc(nuc>0));
    if isempty(labs); continue; end
    rp = regionprops(nuc, 'Area','Centroid','PixelIdxList','BoundingBox');

    if isfinite(p.Results.max_cells_per_xy) && numel(labs) > p.Results.max_cells_per_xy
        labs = labs(randperm(numel(labs), p.Results.max_cells_per_xy));
    end

    [snr, bnd, ar, cen] = deal(nan(numel(labs),1), nan(numel(labs),1), ...
                               nan(numel(labs),1), nan(numel(labs),2));
    gmag = local_gradmag(raw);   % precompute once per XY
    for j = 1:numel(labs)
        L = labs(j);
        px = rp(L).PixelIdxList;
        ar(j)  = rp(L).Area;
        cen(j,:) = rp(L).Centroid;

        in_med = median(raw(px));
        ring   = local_ring_pixels(nuc, L, rp(L).BoundingBox, size(raw));
        ring_med = median(raw(ring));
        snr(j) = in_med - ring_med;                       % signal presence

        % boundary agreement: mean gradient on perimeter / in-mask spread
        perim = bwperim(nuc == L);
        gp = gmag(perim);
        spread = max(1, iqr(double(raw(px))));
        bnd(j) = mean(gp) / spread;                       % boundary correctness
    end

    all_xy=[all_xy; repmat(xy,numel(labs),1)]; %#ok<AGROW>
    all_lab=[all_lab; labs(:)]; all_snr=[all_snr; snr]; %#ok<AGROW>
    all_bnd=[all_bnd; bnd]; all_area=[all_area; ar]; all_cen=[all_cen; cen]; %#ok<AGROW>
end

Q = struct('xy',all_xy,'lab',all_lab,'snr',all_snr,'bnd',all_bnd, ...
           'area',all_area,'centroid',all_cen);
Q.meta = struct('mask_dir',mask_dir,'mask_name',mask_name, ...
                'nuc_chan',nuc_chan,'ip_fname',ip.fname,'ip',ip, ...
                'thumbs_per_bin',p.Results.thumbs_per_bin);

if isempty(all_snr)
    warning('cycif_qc_mask_signal: no cells measured.'); return
end

% ---- Linked-histogram figures -----------------------------------------
local_linked_hist(Q, 'snr', 'signal presence  (in-mask - local ring, raw units)', ...
    'left tail = spurious / empty masks', p.Results.save_dir, 'qc_signal_presence');
local_linked_hist(Q, 'bnd', 'boundary agreement  (perimeter gradient / in-mask IQR)', ...
    'left tail = mask edge not on a real intensity edge', p.Results.save_dir, 'qc_boundary');

% ---- Joint view: the two axes together, colored by area ---------------
f = figure('Color','w','Name','QC signal vs boundary','Position',[120 120 620 520]);
scatter(Q.snr, Q.bnd, 14, log10(max(Q.area,1)), 'filled', ...
        'MarkerFaceAlpha',0.5); colormap(parula);
cb=colorbar; cb.Label.String='log_{10} area';
xlabel('signal presence'); ylabel('boundary agreement');
title('bottom-left = suspect (no signal AND no real boundary)');
grid on;
if ~isempty(char(p.Results.save_dir))
    local_savefig(f, char(p.Results.save_dir), 'qc_signal_vs_boundary');
end
end


% =========================================================================
% LINKED HISTOGRAM: click a bar -> montage of the cells in that bin
% =========================================================================
function local_linked_hist(Q, field, xlab, tail_note, save_dir, base)
v = Q.(field);
v = v(~isnan(v) & isfinite(v));
if numel(v) < 2; warning('cycif_qc_mask_signal: too few finite %s values to plot.', field); return; end

f = figure('Color','w','Name',['QC ' field],'Position',[100 100 760 460]);
edges = local_safe_edges(v, 40);
h = histogram(v, edges, 'FaceColor',[0.35 0.5 0.8]);
xlabel(xlab); ylabel('cells');
title({sprintf('%s  (click any bar to see its cells)', field), tail_note}, ...
      'Interpreter','none');

% Attach click behavior: map click x -> bin -> cell indices -> montage.
ax = gca;
S = struct('Q',Q,'field',field,'edges',edges);
ax.ButtonDownFcn = @(src,evt) local_on_click(evt, S);
h.HitTest = 'off';           % let clicks fall through to the axes
h.PickableParts = 'none';

if ~isempty(char(save_dir))
    local_savefig(f, char(save_dir), base);
end
end

function local_on_click(evt, S)
x = evt.IntersectionPoint(1);
b = find(x >= S.edges(1:end-1) & x < S.edges(2:end), 1);
if isempty(b); return; end
v = S.Q.(S.field);
sel = find(v >= S.edges(b) & v < S.edges(b+1));
if isempty(sel); return; end
n = min(numel(sel), S.Q.meta.thumbs_per_bin);
sel = sel(round(linspace(1, numel(sel), n)));
local_show_montage(S.Q, sel, sprintf('%s bin [%.2f, %.2f)  (%d cells)', ...
    S.field, S.edges(b), S.edges(b+1), numel( find(v>=S.edges(b)&v<S.edges(b+1)) )));
end


% =========================================================================
% MONTAGE of selected cells: raw crop with mask outline (one panel per cell)
% Kept single-panel (not the 3-panel triplet) so a full bin of cells stays
% legible; the outline on raw shows both signal presence and boundary fit.
% =========================================================================
function local_show_montage(Q, sel, ttl)
m = Q.meta;
dao = ImageAccess(m.ip_fname, m.ip.indsz);
ncol = 5; nrow = ceil(numel(sel)/ncol);
f = figure('Color','w','Name',ttl,'Position',[160 80 900 min(1000,180*nrow+60)]);
tl = tiledlayout(f, nrow, ncol, 'TileSpacing','tight','Padding','compact');
title(tl, ttl, 'Interpreter','none');

% Group selected cells by XY so each frame loads once.
xys = unique(Q.xy(sel));
pad = 14;
for u = 1:numel(xys)
    xy = xys(u);
    fmask = local_maskfile(m.mask_dir, m.mask_name, xy);
    Mf = load(fmask,'nuc_labels','last_t');
    nuc = Mf.nuc_labels;
    raw = double(dao.GetFrame([Mf.last_t, m.nuc_chan, xy, 1]));
    lo = prctile(raw(:),1); hi = prctile(raw(:),99.5);
    [H,W]=size(nuc);
    idx = sel(Q.xy(sel)==xy);
    for t = idx(:)'
        L = Q.lab(t); ctr = Q.centroid(t,:);
        rr=round(ctr(2)); cc=round(ctr(1));
        r0=max(1,rr-pad); r1=min(H,rr+pad); c0=max(1,cc-pad); c1=min(W,cc+pad);
        rawC=raw(r0:r1,c0:c1); nucC=nuc(r0:r1,c0:c1); thisC=nucC==L;

        ax=nexttile; imagesc(ax,rawC,[lo hi]); colormap(ax,gray); axis(ax,'image','off');
        hold(ax,'on');
        B=bwboundaries(thisC,'noholes');
        for i=1:numel(B); plot(ax,B{i}(:,2),B{i}(:,1),'-','Color',[0.9 0.2 0.2],'LineWidth',1.3); end
        title(ax,sprintf('xy%d s=%.1f b=%.2f',xy,Q.snr(t),Q.bnd(t)),'FontSize',7);
    end
end
end


% =========================================================================
% LOCAL HELPERS
% =========================================================================
function edges = local_safe_edges(v, nbin)
% Return a finite, strictly increasing edge vector even when the data are
% degenerate (all-equal, tied percentiles, or heavy spikes). histogram's
% auto-heuristic errors on non-monotonic/degenerate edges, so guarantee them.
v = v(isfinite(v));
lo = prctile(v, 0.5); hi = prctile(v, 99.5);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = min(v); hi = max(v);                 % fall back to full range
end
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    c = v(1); lo = c - 0.5; hi = c + 0.5;      % all identical -> unit window
end
edges = linspace(lo, hi, max(2, nbin));
end

function ring = local_ring_pixels(nuc, L, bbox, sz)
% Pixels in a dilated collar around nucleus L that belong to NO nucleus
% (true local background), for a per-cell signal-presence baseline.
pad=6;
x0=max(1,floor(bbox(1))-pad); y0=max(1,floor(bbox(2))-pad);
x1=min(sz(2),ceil(bbox(1)+bbox(3))+pad); y1=min(sz(1),ceil(bbox(2)+bbox(4))+pad);
sub = nuc(y0:y1, x0:x1);
this = sub==L;
grown = imdilate(this, strel('disk',pad));
ringmask = grown & (sub==0);          % collar, background only
if ~any(ringmask(:)); ringmask = sub==0; end
[ry,rx]=find(ringmask);
ring = sub2ind(sz, ry+y0-1, rx+x0-1);
end

function g = local_gradmag(raw)
[gx,gy]=gradient(raw); g=hypot(gx,gy);
end

function c = local_resolve_chan(op)
% Same provider-aware resolution as ct_save_endpoint_masks; returns the raw
% acquisition channel to hand to GetFrame (legacy maps plane->op.cind).
c = NaN;
if ~isfield(op,'seg'); return; end
if strcmpi(op.seg.provider,'cellpose')
    if isfield(op.seg,'cp') && isfield(op.seg.cp,'chan_nuc') && ~isempty(op.seg.cp.chan_nuc)
        c = double(op.seg.cp.chan_nuc);
    end
    return
end
if isfield(op.seg,'chan') && ~isempty(op.seg.chan)
    plane=double(op.seg.chan(1));
    if isfield(op,'cind') && ~isempty(op.cind) && plane>=1 && plane<=numel(op.cind)
        c=double(op.cind(plane));
    else
        c=plane;
    end
end
end

function f = local_maskfile(mask_dir, mask_name, xy)
f = fullfile(mask_dir, sprintf('%s_xy%02d_endpoint_masks.mat', mask_name, xy));
end

function local_savefig(f, save_dir, base)
if ~exist(save_dir,'dir'); mkdir(save_dir); end
exportgraphics(f, fullfile(save_dir,[base '.png']), 'Resolution',150);
end