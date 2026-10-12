function figs = dl_state_summary(varargin)
% dl_state_summary  Compact, publication-style summary figures for the state
%   analysis suite: one figure each for stain QC, Q1 (treatment dispersion),
%   Q2 (functional dimensionality) and Q3 (dynamics -> function). Replaces the
%   per-run diagnostic plots (call the analysis functions with 'plots','none').
%
% USAGE
%   figs = dl_state_summary('qc',QC, 'q1',R1, 'q1sum',Q1sum, 'q2',R2, 'q3',R3, 'savepath',saveDir);
%   figs = dl_state_summary('q3',R3, 'q3pc',R3pc, 'q3bg',R3bg, 'save',false);  % Q3 + amplitude + background controls
%   figs = dl_state_summary('q1',dl_treatment_dispersion(Tf,'plots','none'), 'save',false);
%   figs = dl_state_summary('q3focus',R3f, 'save',false);          % focused markers: PLS vs kNN, within-well
%   figs = dl_state_summary('latentmap',LM, 'save',false);         % which sensor each latent dim encodes
%   figs = dl_state_summary('q3',R3, 'q3bg',R3bg, 'poster','partition', 'poster_width_cm',9, 'savepath',saveDir);  % poster only + summary
%   K = table(["TNFa";"IL1b";"IL8";"IL10RA";"HLA_DR"], ["producer";"producer";"producer";"diffuse";"floor"], ...
%             'VariableNames', {'marker','kind'});
%   figs = dl_state_summary('qc',QC, 'q3',R3, 'producers',R3t, 'kinds',K, 'save',false);   % rare-producer view
%
% NAME-VALUE (all optional; a figure is drawn for each input supplied)
%   'qc'        struct   stain QC from section 8 (.summary, .pct_pos_by_stimulus)
%   'q1'        struct   R1 composite (.stimulus .glucose .condition) or one dl_treatment_dispersion result
%   'q1sum'     table    grouping | between_frac | well_frac | disp_ratio | n
%   'q2'        struct   dl_functional_dimensionality result
%   'q3'        struct   dl_dynamics_predict result (marker targets)
%   'q3bg'      struct   dl_dynamics_predict on secondary-only cells (technical-coupling control)
%   'q3pc'      struct   dl_dynamics_predict result with 'targets','pcs','drop_pcs',[] (amplitude check)
%   'q3focus'   struct   dl_dynamics_predict result run with 'nonlinear','knn' (per-marker distribution view)
%   'latentmap' struct   .R2 [L x S] variance of each latent dim explained by each sensor's trace
%                        .cols (S names) .dims (L names) .sets (R.sets of a 'vip' run, optional)
%   'producers' struct   dl_dynamics_predict result run with 'tail' (R.tail): producer frequency by
%                        condition, same-well identity AUC, gate sanity, co-production
%   'kinds'     table    marker | kind ("producer" | "diffuse" | "floor" | "graded"): stain type,
%                        tagged on marker labels in every figure (P / D / F) and colours the QC bars
%   'poster'    'both' | 'partition' | 'comparison' | 'none': poster figures from 'q3' (+ 'q3bg'),
%                        fixed size, vector PDF + SVG                                   ('both')
%                        partition  = where the predictable variance comes from (dynamics only /
%                                     either / treatment only), black tick = secondary-only floor
%                        comparison = treatment vs dynamics (pooled, within well) vs shuffled null
%   'poster_targets' ({} = all, order kept when 'poster_sort' false)  'poster_sort' (true)
%   'poster_labels' ({name, display; ...} TeX)  'poster_alpha' (0.05)  'poster_width_cm' (12)
%   'poster_font_size' (10)  'poster_font' ('Arial')
%   'eta_lo'    scalar   reference line on variance fractions          (0.10)
%   'tag'       char     appended to file names (e.g. winAlign)        ('')
%   'save'      logical  save SVGs                                     (true)
%   'closefigs' logical  close after saving                            (true)
%   'savepath'  char     output folder                                 ('')
%
% OUTPUT
%   figs  figure handles (empty if closed after saving)
%
% NOTES
%   New in DatalocHandler v5. Reads only result structs; no recomputation.

p = inputParser;  p.FunctionName = mfilename;
addParameter(p, 'qc',        []);   addParameter(p, 'q1',    []);
addParameter(p, 'q1sum',     []);   addParameter(p, 'q2',    []);
addParameter(p, 'q3',        []);   addParameter(p, 'q3pc',  []);  addParameter(p, 'q3bg', []);
addParameter(p, 'q3focus',   []);
addParameter(p, 'latentmap', []);
addParameter(p, 'producers', []);
addParameter(p, 'kinds',     [],    @(x) isempty(x) || istable(x));
addParameter(p, 'eta_lo',    0.10);
addParameter(p, 'tag',       '',    @(x) ischar(x) || isstring(x));
addParameter(p, 'save',      true);
addParameter(p, 'closefigs', true);
addParameter(p, 'savepath',  '',    @(x) ischar(x) || isstring(x));
addParameter(p, 'poster',          'both', @(x) any(strcmpi(x, {'both','partition','comparison','none'})));
addParameter(p, 'poster_targets',  {});
addParameter(p, 'poster_sort',     true);
addParameter(p, 'poster_labels',   {},      @iscell);
addParameter(p, 'poster_alpha',    0.05);
addParameter(p, 'poster_width_cm', 12);
addParameter(p, 'poster_font_size', 10);
addParameter(p, 'poster_font',     'Arial', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});  o = p.Results;  o.tag = char(o.tag);  o.savepath = char(o.savepath);
if o.save && isempty(o.savepath), warning('dl_state_summary: no savepath; figures not saved.'); end

figs = gobjects(0);
K = o.kinds;                                                       % marker kind tags (may be empty)
if ~isempty(o.qc), figs(end+1) = i_fig_qc(o.qc, K);              i_save(figs(end), o, 'summary-qc');  end
if ~isempty(o.q1), figs(end+1) = i_fig_q1(o.q1, o.q1sum, o, K);  i_save(figs(end), o, 'summary-q1');  end
if ~isempty(o.q2), figs(end+1) = i_fig_q2(o.q2);                 i_save(figs(end), o, 'summary-q2');  end
if ~isempty(o.q3), figs(end+1) = i_fig_q3(o.q3, o.q3pc, o.q3bg, K); i_save(figs(end), o, 'summary-q3');  end
if ~isempty(o.q3) && ~strcmpi(o.poster, 'none')                   % poster versions of Q3 (PDF + SVG)
    figs = [figs, i_poster(o.q3, o.q3bg, o)];
end
if ~isempty(o.q3focus), figs(end+1) = i_fig_q3focus(o.q3focus); i_save(figs(end), o, 'summary-q3focus');  end
if ~isempty(o.producers)
    assert(isfield(o.producers, 'tail') && ~isempty(o.producers.tail), ...
           'dl_state_summary: ''producers'' needs a dl_dynamics_predict result run with ''tail''.');
    figs(end+1) = i_fig_producers(o.producers.tail, K);            i_save(figs(end), o, 'summary-producers');
end
if ~isempty(o.latentmap), figs(end+1) = i_fig_latentmap(o.latentmap); i_save(figs(end), o, 'summary-latentmap');  end
if o.closefigs && o.save && ~isempty(o.savepath), close(figs);  figs = gobjects(0); end
end


%% =========================================================================
%  Stain QC
%% =========================================================================
function f = i_fig_qc(QC, K)
% The secondary-only control removes autofluorescence + secondary binding, NOT diffuse
% primary binding: SBR and "% above bg p99" measure the BULK shift with antibody, not
% the fraction of positive cells. Producer markers are judged on their bright tail
% (summary-producers); here they are only tagged.
f  = i_newfig('Stain QC', [1150 400]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {'Stain QC against secondary-only background', ...
           'Bulk shift with antibody (not % positive): diffuse stains look bright here, rare producers look dim'}, ...
      'FontWeight', 'bold', 'FontSize', 10.5);

% Bulk shift, coloured by stain kind (grey = <10% of bulk above bg p99 when no kind given)
ax = nexttile(tl);  hold(ax, 'on');
Sm = QC.summary;  [~, ord] = sort(Sm.SBR_median);  Sm = Sm(ord, :);
kd = i_kind(Sm.marker, K);
if all(kd == ""), kd(Sm.pctpos_max_rep < 0.10) = "floor"; end
CD = zeros(height(Sm), 3);
for i = 1:height(Sm)
    switch kd(i)
        case "producer", CD(i, :) = [0.80 0.30 0.30];
        case "diffuse",  CD(i, :) = [0.90 0.60 0.20];
        case "floor",    CD(i, :) = [0.75 0.75 0.75];
        otherwise,       CD(i, :) = [0.20 0.45 0.70];
    end
end
b = barh(ax, max(Sm.SBR_median, 0.05), 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = CD;
set(ax, 'XScale', 'log', 'YTick', 1:height(Sm), 'YTickLabel', i_tagged(Sm.marker, K), 'TickLabelInterpreter', 'none');
xline(ax, 1, ':k');
xlabel(ax, 'bulk median above background (bg SDs)');
title(ax, {'Bulk shift per marker', 'red producer | orange diffuse | grey floor | blue graded'}, 'FontSize', 8.5);
i_style(ax);

% Bulk above background p99, by stimulus
ax = nexttile(tl, [1 2]);
Pt = QC.pct_pos_by_stimulus;
i_heat(ax, Pt{:,:}, i_tagged(Pt.Properties.VariableNames, K), Pt.Properties.RowNames, i_seq(256), [0 100], '%.0f');
title(ax, {'% cells above secondary-only p99 (bulk measure)', ...
           'P = producer: see summary-producers for the bright-tail frequency;  D = diffuse;  F = floor'}, 'FontSize', 8.5);
end


%% =========================================================================
%  Q1: does treatment define state?
%% =========================================================================
function f = i_fig_q1(R1, Q1sum, o, K)
if isfield(R1, 'stimulus'), Rs = R1.stimulus;  Rg = i_getf(R1, 'glucose');  else, Rs = R1;  Rg = []; end
f  = i_newfig('Q1 treatment dispersion', [1250 430]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q1  Does treatment define functional state?   verdict (stimulus): %s', Rs.verdict), ...
      'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');

% (a) variance partition per grouping
ax = nexttile(tl);
if isempty(Q1sum)
    Q1sum = table("grouping", Rs.partition.frac(1), Rs.partition.frac(min(2, end)), Rs.ratio_global, Rs.n, ...
                  'VariableNames', {'grouping','between_frac','well_frac','disp_ratio','n'});
end
P = [Q1sum.between_frac, Q1sum.well_frac, 1 - Q1sum.between_frac - Q1sum.well_frac];
b = barh(ax, flipud(P), 'stacked', 'EdgeColor', 'none', 'BarWidth', 0.7);
cols = [0.20 0.45 0.70; 0.90 0.60 0.20; 0.86 0.86 0.86];
for i = 1:3, b(i).FaceColor = cols(i, :); end
set(ax, 'YTick', 1:height(Q1sum), 'YTickLabel', flipud(Q1sum.grouping), 'TickLabelInterpreter', 'none');
xline(ax, o.eta_lo, ':k');  xlim(ax, [0 1]);
legend(ax, {'between groups', 'wells within group', 'cells'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlabel(ax, 'fraction of functional variance');  title(ax, 'Variance partition');
i_style(ax);

% (b) per-marker eta^2: stimulus, glucose, wells
ax = nexttile(tl);  hold(ax, 'on');
E = Rs.eta2;  [~, ord] = sort(E.eta2);  E = E(ord, :);  y = 1:height(E);
plot(ax, E.eta2_unit, y, 'o', 'Color', [0.90 0.60 0.20], 'MarkerSize', 6, 'LineWidth', 1.2);
plot(ax, E.eta2,      y, 'o', 'MarkerFaceColor', [0.20 0.45 0.70], 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
lg = {'wells within stimulus', 'stimulus'};
if ~isempty(Rg) && all(ismember(E.marker, Rg.eta2.marker))
    [~, ig] = ismember(E.marker, Rg.eta2.marker);
    plot(ax, Rg.eta2.eta2(ig), y, 'd', 'MarkerFaceColor', [0.30 0.65 0.35], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
    lg{end+1} = 'glucose';
end
set(ax, 'YTick', y, 'YTickLabel', i_tagged(E.marker, K), 'TickLabelInterpreter', 'none');  ylim(ax, [0.5 numel(y)+0.5]);
xline(ax, o.eta_lo, ':k');
legend(ax, lg, 'Location', 'southeast', 'Box', 'off');
xlabel(ax, '\eta^2 (fraction of marker variance)');
ttl = 'Per-marker variance explained';
if any(i_kind(E.marker, K) == "producer"), ttl = {ttl, 'P: intensity eta^2 under-reads a change in producer FREQUENCY'}; end
title(ax, ttl, 'FontSize', 8.5);
i_style(ax);

% (c) centroid distances (off-diagonal) vs within-group spread (diagonal)
ax = nexttile(tl);
D = Rs.centroid_dist;  D(1:size(D,1)+1:end) = Rs.group.spread;
i_heat(ax, D, Rs.groups, Rs.groups, i_seq(256), [0 max(D(:))], '%.2f');
title(ax, sprintf('Centroid distance; diagonal = within spread (ratio %.2f)', Rs.ratio_global));
end


%% =========================================================================
%  Q2: structure of functional space
%% =========================================================================
function f = i_fig_q2(R2)
f  = i_newfig('Q2 functional dimensionality', [1150 760]);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Q2  Structure of functional space   verdict: %s (modes: %s)', R2.verdict, R2.modes_source), ...
      'FontWeight', 'bold', 'FontSize', 11, 'Interpreter', 'none');
kept = R2.kept_pcs;  nk = numel(kept);

% (a) scree of kept PCs vs parallel-analysis null
ax = nexttile(tl);  hold(ax, 'on');
lam = R2.eigenvalues(kept);
bar(ax, 1:nk, lam, 'FaceColor', [0.20 0.45 0.70], 'EdgeColor', 'none', 'BarWidth', 0.7);
plot(ax, 1:nk, R2.pa_null95(1:nk), 'k--', 'LineWidth', 1);
set(ax, 'XTick', 1:nk, 'XTickLabel', compose('PC%d', kept));
legend(ax, {'eigenvalue', 'permutation null 95%'}, 'Box', 'off');
ylabel(ax, 'eigenvalue');
title(ax, sprintf('Scree (removed PC%s)  PR = %.1f of %d  (PR/n %.2f)', mat2str(R2.drop_pcs), R2.pr, nk, R2.pr_frac));
i_style(ax);

% (b) loadings + M1-M2 contrast
ax = nexttile(tl);
ns = min(6, nk);  L = R2.loadings(:, kept(1:ns));
xl = compose('PC%d', kept(1:ns));
if any(R2.m1m2_contrast), L = [L, R2.m1m2_contrast / norm(R2.m1m2_contrast)];  xl = [xl, {'M1-M2'}]; end
lbl = cellstr(R2.markers);  ip = ismember(lbl, R2.pulsatile);  lbl(ip) = strcat(lbl(ip), ' *');
cl = max(abs(L(:)));
i_heat(ax, L, xl, lbl, i_div(256), [-cl cl], '%.2f');
title(ax, sprintf('Loadings (* pulsatile); best M1-M2 match PC%d, cos %.2f', R2.m1m2_best_pc, R2.m1m2_best_cos));

% (c) uniqueness
ax = nexttile(tl);
U = R2.uniqueness;  [~, ord] = sort(U.uniqueness);  U = U(ord, :);
b = barh(ax, U.uniqueness, 'FaceColor', 'flat', 'EdgeColor', 'none', 'BarWidth', 0.7);
b.CData = repmat([0.20 0.45 0.70], height(U), 1);  b.CData(U.pulsatile, :) = repmat([0.80 0.30 0.30], nnz(U.pulsatile), 1);
set(ax, 'YTick', 1:height(U), 'YTickLabel', U.marker, 'TickLabelInterpreter', 'none');  xlim(ax, [0 1]);
xlabel(ax, '1 - R^2 from all other markers');  title(ax, 'Uniqueness (red: pulsatile)');
i_style(ax);

% (d) modality: GMM criteria + dip-test tallies
ax = nexttile(tl);  hold(ax, 'on');
G = R2.gmm;
plot(ax, G.K, G.BIC - G.BIC(1), 'o-', 'Color', [0.20 0.45 0.70], 'LineWidth', 1.2, 'MarkerFaceColor', [0.20 0.45 0.70]);
plot(ax, G.K, G.ICL - G.ICL(1), 's-', 'Color', [0.80 0.30 0.30], 'LineWidth', 1.2, 'MarkerFaceColor', [0.80 0.30 0.30]);
yline(ax, 0, ':k');
legend(ax, {sprintf('\\DeltaBIC (min K=%d)', R2.k_bic), sprintf('\\DeltaICL (min K=%d)', R2.k_icl)}, 'Box', 'off', 'Location', 'best');
nsigA = nnz(R2.dip_axes.p_bonf < 0.05);  nsigW = nnz(R2.dip_within.p_bonf < 0.05);
text(ax, 0.03, 0.06, sprintf('dip test: %d/%d pooled axes, %d/%d within-condition tests multimodal', ...
     nsigA, height(R2.dip_axes), nsigW, height(R2.dip_within)), 'Units', 'normalized', 'FontSize', 8);
xlabel(ax, 'GMM components K');  ylabel(ax, 'criterion relative to K = 1');
title(ax, 'Modality (lower = preferred)');
i_style(ax);
end


%% =========================================================================
%  Q3: dynamics -> function
%% =========================================================================
function f = i_fig_q3(R3, R3pc, R3bg, K)
hasPC = ~isempty(R3pc);
f  = i_newfig('Q3 dynamics predict function', [1250 420 + 330*hasPC]);
tl = tiledlayout(f, 1 + hasPC, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = 'Q3  Do live dynamics predict function beyond treatment?';
if isfield(R3, 'resolution') && ~(R3.resolution.q2_resolved && R3.resolution.cmi_resolved)
    ttl = sprintf('%s   (UNRESOLVED: nperm < %d)', ttl, R3.resolution.nperm_needed);
end
title(tl, ttl, 'FontWeight', 'bold', 'FontSize', 11);
Q = R3.q2;  [~, ord] = sort(Q.q2_dyn, 'descend');  Q = Q(ord, :);  C = R3.cmi(ord, :);  x = 1:height(Q);

% (a) Q2 horse race
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, x, Q.q2_tx,          'o', 'MarkerFaceColor', [0.55 0.55 0.55], 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
plot(ax, x, Q.q2_dyn,         'o', 'MarkerFaceColor', [0.80 0.30 0.30], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
plot(ax, x, Q.q2_txdyn,       'o', 'MarkerFaceColor', [0.50 0.30 0.70], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
plot(ax, x, Q.q2_best_marker, 'o', 'Color', [0.30 0.65 0.35], 'MarkerSize', 7, 'LineWidth', 1.2);
plot(ax, x, Q.q2_null95,      '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);
lgq = {'treatment', 'dynamics', 'treatment + dynamics', 'best single marker', 'dynamics null 95%'};
if ~isempty(R3bg)                         % same targets measured in secondary-only cells
    [hit, ib] = ismember(Q.target, R3bg.q2.target);
    xb = x(hit);  plot(ax, xb(:), R3bg.q2.q2_dyn(ib(hit)), 'x', 'Color', [0.1 0.1 0.1], 'MarkerSize', 8, 'LineWidth', 1.5);
    lgq{end+1} = 'dynamics, secondary-only cells';
end
yline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
set(ax, 'XTick', x, 'XTickLabel', i_tagged(Q.target, K), 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
legend(ax, lgq, ...
       'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
ylabel(ax, 'Q^2Y (outer CV, new wells)');
ttl = 'Predictability by view';
if any(i_kind(Q.target, K) == "producer"), ttl = {ttl, 'P: intensity Q^2 mostly reflects the bulk; producer identity in summary-producers'}; end
title(ax, ttl, 'FontSize', 8.5);
i_style(ax);

% (b) conditional MI vs within-treatment null
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, x, C.null95, '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);
plot(ax, x, C.cmi, 'o', 'MarkerFaceColor', [0.80 0.30 0.30], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
unr = Q.verdict == "unresolved (raise nperm)";
if any(unr), xu = x(unr);  text(ax, xu(:), C.cmi(unr), '  +', 'FontSize', 9, 'FontWeight', 'bold'); end
set(ax, 'XTick', x, 'XTickLabel', i_tagged(Q.target, K), 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
legend(ax, {'null 95%', 'I(dyn; y | condition)'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 7);
ylabel(ax, 'nats');  title(ax, 'Conditional MI (+ = unresolved)', 'Interpreter', 'none');
i_style(ax);

% (c) single-marker cross-prediction
ax = nexttile(tl);
if ~isempty(R3.single_q2) && any(~isnan(R3.single_q2(:)))
    i_heat(ax, R3.single_q2(:, ord), cellstr(Q.target), R3.single_names, i_seq(256), [0 max(0.01, max(R3.single_q2(:)))], '%.2f');
    xlabel(ax, 'target');  ylabel(ax, 'predictor');
    title(ax, 'Single-marker Q^2Y (shared amplitude inflates these)');
else
    axis(ax, 'off');
end

% (d) amplitude check: dynamics -> functional PCs incl. PC1
if ~isempty(R3pc)
    ax = nexttile(tl, [1 3]);  hold(ax, 'on');
    P = R3pc.q2;  xp = 1:height(P);
    bar(ax, xp, P.q2_dyn, 0.5, 'FaceColor', [0.80 0.30 0.30], 'EdgeColor', 'none');
    plot(ax, xp, P.q2_null95, '_', 'Color', 'k', 'MarkerSize', 18, 'LineWidth', 1.5);
    plot(ax, xp, P.q2_tx, 'o', 'MarkerFaceColor', [0.55 0.55 0.55], 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
    set(ax, 'XTick', xp, 'XTickLabel', P.target, 'TickLabelInterpreter', 'none');
    legend(ax, {'dynamics', 'null 95%', 'treatment'}, 'Box', 'off', 'Location', 'northeast');
    ylabel(ax, 'Q^2Y');
    title(ax, 'Amplitude check: what do dynamics predict? (fPC1 = uniform-sign amplitude axis)');
    i_style(ax);
end
end


%% =========================================================================
%  Q3 focus: does the latent order cells within a well? (no thresholding)
%% =========================================================================
function f = i_fig_q3focus(R)
% Rows = targets. Cols: observed distribution by within-well predicted quintile
% (PLS, kNN), then Q2 pooled vs within-well against nulls.
nT  = numel(R.pred);  hasNL = isfield(R, 'nl') && any(~isnan(R.nl.q2_knn));
nc  = 2 + hasNL;
f   = i_newfig('Q3 focus', [380*nc 230*nT + 60]);
tl  = tiledlayout(f, nT, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {'Do sensor dynamics order cells within a well?', ...
           'Observed intensity, split by quintile of outer-CV prediction ranked within each well'}, ...
      'FontWeight', 'bold', 'FontSize', 11);
cq  = [linspace(0.80, 0.10, 5)' linspace(0.86, 0.30, 5)' linspace(0.93, 0.60, 5)'];   % light -> dark
for t = 1:nT
    P = R.pred(t);  y = P.y;
    preds = {P.yhat};  pl = {'PLS'};
    if hasNL, preds{end+1} = P.yhat_nl;  pl{end+1} = 'kNN'; end
    for m = 1:numel(preds)
        ax = nexttile(tl);  hold(ax, 'on');
        qb = i_wquint(preds{m}, P.W);
        g  = linspace(prctile(y, 0.5), prctile(y, 99.5), 200);
        for q = 1:5
            yq = y(qb == q);  if numel(yq) < 20, continue; end
            d  = ksdensity(yq, g);  d = 0.9 * d / max(d);
            fill(ax, [g fliplr(g)], [d zeros(size(d))] + (q - 1), cq(q,:), 'EdgeColor', cq(q,:)*0.8, ...
                 'FaceAlpha', 0.75, 'LineWidth', 0.75);
            plot(ax, median(yq), q - 1, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
        end
        set(ax, 'YTick', (0:4) + 0.35, 'YTickLabel', {'Q1 low','Q2','Q3','Q4','Q5 high'});
        ylim(ax, [-0.1 5]);  xlim(ax, g([1 end]));
        if t == nT, xlabel(ax, 'observed intensity (per-rep z)'); end
        title(ax, sprintf('%s  |  %s', R.targets{t}, pl{m}), 'Interpreter', 'none');
        i_style(ax);  ax.YGrid = 'off';
    end
    % Q2 summary: pooled vs within-well, PLS vs kNN, against null 95th
    ax = nexttile(tl);  hold(ax, 'on');
    N  = R.nl(t, :);
    v  = [N.q2_lin N.q2w_lin; N.q2_knn N.q2w_knn];
    n9 = [N.q2_lin_null95 N.q2w_lin_null95; N.q2_knn_null95 N.q2w_knn_null95];
    if ~hasNL, v = v(1,:);  n9 = n9(1,:); end
    b = bar(ax, v', 'grouped', 'EdgeColor', 'none', 'BarWidth', 0.8);
    cols = [0.80 0.30 0.30; 0.20 0.45 0.70];
    for m = 1:numel(b)
        b(m).FaceColor = cols(m,:);
        xe = b(m).XEndPoints;
        hn = plot(ax, xe, n9(m,:), '_', 'Color', 'k', 'MarkerSize', 12, 'LineWidth', 1.5);
    end
    yline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
    set(ax, 'XTick', 1:2, 'XTickLabel', {'pooled', 'within-well'});
    ylabel(ax, 'Q^2 (outer CV)');
    if t == 1, legend(ax, [b(:)', hn], [pl, {'null 95%'}], 'Box', 'off', 'Location', 'northeast', 'FontSize', 7); end
    title(ax, sprintf('%s  (kNN k = %g)', R.targets{t}, N.knn_k), 'Interpreter', 'none');
    i_style(ax);
end
end

function qb = i_wquint(yh, W)
% Quintile of the prediction ranked within each well (NaN for wells < 5 cells)
qb = nan(size(yh));
if isempty(W), W = ones(size(yh)); end
for w = unique(W(:))'
    i = find(W == w & ~isnan(yh));
    if numel(i) < 5, continue; end
    qb(i) = ceil(5 * tiedrank(yh(i)) / numel(i));
end
end

%% =========================================================================
%  Rare producers: how many (treatment) and which (dynamics)
%% =========================================================================
function f = i_fig_producers(TR, K)
tg = string(TR.targets(:))';  nt = numel(tg);
cR = [0.80 0.30 0.30];  cB = [0.20 0.45 0.70];  cG = [0.55 0.55 0.55];
f  = i_newfig('Rare producers', [1350 760]);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, {sprintf('Rare producers  (gate: replicate median + %g robust SD of stained cells; wells are the unit for frequency)', TR.kref), ...
           'Does treatment set HOW MANY cells produce?  Do live dynamics pick WHICH cells produce, within a well?'}, ...
      'FontWeight', 'bold', 'FontSize', 10.5);

% (a, b) producer frequency per level of the first two breakdown variables
vars = strings(0, 1);  if ~isempty(TR.freq), vars = unique(TR.freq.var, 'stable'); end
for v = 1:2
    ax = nexttile(tl);
    if v > numel(vars) || isempty(TR.freq), axis(ax, 'off'); continue; end
    Fv = TR.freq(TR.freq.var == vars(v), :);  lv = unique(Fv.level, 'stable');
    M  = nan(numel(lv), nt);
    for t = 1:nt
        Ft = Fv(Fv.target == tg(t), :);  [~, ix] = ismember(Ft.level, lv);  M(ix, t) = Ft.pct_mean;
    end
    E  = TR.freq_effect(TR.freq_effect.var == vars(v), :);  [~, ie] = ismember(tg, E.target);
    xl = tg;  ok = ie > 0;
    if any(ok)
        tgo = tg(ok);  e2 = E.eta2(ie(ok));  st = i_star(E.p(ie(ok)), TR.alpha);
        xl(ok) = compose('%s  eta2 %.2f%s', tgo(:), e2(:), st(:))';
    end
    i_heat(ax, M, cellstr(xl), cellstr(lv), i_seq(256), [0 max(1, max(M(:)))], '%.1f');
    title(ax, {sprintf('%% producers per well, by %s', vars(v)), 'eta2 = share of well-to-well frequency variance (* p < alpha)'}, ...
          'FontSize', 8.5, 'Interpreter', 'none');
end

% (c) identity: same-well AUC
ax = nexttile(tl);  hold(ax, 'on');
Tb = TR.table(TR.table.k == TR.kref, :);  [~, ix] = ismember(tg, Tb.target);  Tb = Tb(ix, :);  x = 1:nt;
h = gobjects(0);  lg = {};
h(end+1) = plot(ax, x, Tb.auc_tx, 'o', 'MarkerFaceColor', cG, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);   lg{end+1} = 'treatment frequency (pooled)';
h(end+1) = plot(ax, x, Tb.auc_lin, 'o', 'MarkerFaceColor', cR, 'MarkerEdgeColor', 'none', 'MarkerSize', 7);  lg{end+1} = 'dynamics PLS (pooled)';
h(end+1) = plot(ax, x, Tb.auc_lin_w, 'd', 'Color', cR, 'MarkerSize', 9, 'LineWidth', 1.4);                    lg{end+1} = 'dynamics PLS, same well';
if any(~isnan(Tb.auc_knn_w))
    h(end+1) = plot(ax, x, Tb.auc_knn_w, 'd', 'Color', cB, 'MarkerSize', 9, 'LineWidth', 1.4);                lg{end+1} = 'dynamics kNN, same well';
end
h(end+1) = plot(ax, x, Tb.auc_lin_w_null95, '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5);          lg{end+1} = 'same-well shuffle 95%';
yline(ax, 0.5, '-', 'Color', [0.7 0.7 0.7]);
set(ax, 'XTick', x, 'XTickLabel', i_tagged(tg, K), 'TickLabelInterpreter', 'none');  xlim(ax, [0.5 nt + 0.5]);
ylabel(ax, 'AUC, producer vs non-producer');  legend(ax, h, lg, 'Box', 'off', 'FontSize', 7, 'Location', 'best');
title(ax, {'Which cells produce?', 'Same-well AUC: treatment and well removed; 0.5 = no information'}, 'FontSize', 8.5);
i_style(ax);

% (d) gate sanity across k
ax = nexttile(tl);  hold(ax, 'on');
ks = TR.k;  nK = numel(ks);  w = 0.8 / nK;  cK = gray(nK + 2);  cK = cK(2:end-1, :);
for kk = 1:nK
    Tk = TR.table(TR.table.k == ks(kk), :);  [~, ix] = ismember(tg, Tk.target);  xx = (1:nt) - 0.4 + w * (kk - 0.5);
    bar(ax, xx, Tk.pct_tail(ix), w, 'FaceColor', cK(kk, :), 'EdgeColor', 'none');
    if TR.has_bg, plot(ax, xx, max(Tk.bg_fpr(ix), 1e-3), '_', 'Color', cR, 'MarkerSize', 8, 'LineWidth', 1.6); end
end
set(ax, 'XTick', 1:nt, 'XTickLabel', tg, 'YScale', 'log', 'TickLabelInterpreter', 'none');  xlim(ax, [0.5 nt + 0.5]);
ylabel(ax, '% cells above gate');
title(ax, {sprintf('Gate sanity: k = %s (light to dark)', mat2str(ks)), 'red = secondary-only cells passing the same gate'}, 'FontSize', 8.5);
i_style(ax);

% (e) co-production
ax = nexttile(tl);
C  = TR.copos;  Lg = log2(C.enrich);  Lg(1:size(Lg, 1) + 1:end) = NaN;
cl = max(1, max(abs(Lg(:)), [], 'omitnan'));
i_heat(ax, Lg, i_tagged(C.names, K), i_tagged(C.names, K), i_div(256), [-cl cl], '%.1f');
title(ax, {'Same cells? log_2 co-call enrichment (top tails)', 'producer x producer high = co-producers; with D/F/graded = artefact suspect'}, 'FontSize', 8.5);

% (f) verdicts
ax = nexttile(tl);  axis(ax, 'off');
V = TR.verdict;  txt = {sprintf('VERDICT (k = %g)', TR.kref), ''};
for i = 1:height(V)
    txt{end+1} = char(V.target(i)); %#ok<AGROW>
    txt = [txt, strcat({'   '}, i_wrap(char(V.frequency(i)), 48)), strcat({'   '}, i_wrap(char(V.identity(i)), 48))]; %#ok<AGROW>
    if V.gate(i) ~= "clean", txt{end+1} = ['   gate: ' char(V.gate(i))]; end %#ok<AGROW>
end
text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8.5, 'Interpreter', 'none');
end

function s = i_star(p, a)
s = repmat("", size(p));  s(p < a) = " *";
end

function c = i_wrap(s, w)
% Word-wrap a char vector into a row cellstr of lines <= w characters
words = strsplit(s, ' ');  c = {};  ln = '';
for i = 1:numel(words)
    if ~isempty(ln) && numel(ln) + numel(words{i}) + 1 > w, c{end+1} = ln; ln = words{i}; %#ok<AGROW>
    elseif isempty(ln), ln = words{i};
    else, ln = [ln ' ' words{i}]; %#ok<AGROW>
    end
end
if ~isempty(ln), c{end+1} = ln; end
end

function kd = i_kind(names, K)
% Stain kind per marker name ("" when unknown)
names = string(names(:));  kd = repmat("", numel(names), 1);
if isempty(K), return; end
[hit, loc] = ismember(lower(regexprep(names, '_(cyt|nuc)$', '')), lower(string(K.marker)));
kd(hit) = lower(string(K.kind(loc(hit))));
end

function lab = i_tagged(names, K)
% Marker labels with kind tags: (P) producer, (D) diffuse, (F) floor
kd  = i_kind(names, K);  lab = string(names(:));
tagm = containers.Map({'producer','diffuse','floor'}, {' (P)', ' (D)', ' (F)'});
for i = 1:numel(lab)
    if strlength(kd(i)) > 0 && isKey(tagm, char(kd(i))), lab(i) = lab(i) + tagm(char(kd(i))); end
end
lab = cellstr(lab);
end


%% =========================================================================
%  Latent map: which sensor does each latent dimension encode?
%% =========================================================================
function f = i_fig_latentmap(LM)
L  = numel(LM.dims);  S = numel(LM.cols);
ord = 1:L;  edges = [];  gl = {};
if isfield(LM, 'sets') && ~isempty(LM.sets) && isfield(LM.sets, 'members')
    ord = [];
    for g = 1:numel(LM.sets.members)
        ord = [ord, find(ismember(LM.dims, LM.sets.members{g}))]; %#ok<AGROW>
        edges(end+1) = numel(ord); %#ok<AGROW>
    end
    ord = [ord, setdiff(1:L, ord)];
    gl  = LM.sets.names;
end
f  = i_newfig('Latent map', [max(700, 34 * L + 260) 360]);
ax = axes(f);
M  = LM.R2(ord, :)';
imagesc(ax, M, [0 1]);  colormap(ax, i_seq(256));  cb = colorbar(ax);  cb.Box = 'off';  cb.Label.String = 'R^2';
for i = 1:S
    for j = 1:L
        text(ax, j, i, sprintf('%.2f', M(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 6.5, ...
             'Color', [0 0 0] + (M(i, j) > 0.6));
    end
end
for e = edges(1:end-1), xline(ax, e + 0.5, 'k-', 'LineWidth', 1.5); end
for g = 1:numel(edges)
    x0 = 1;  if g > 1, x0 = edges(g-1) + 1; end
    text(ax, (x0 + edges(g)) / 2, 0.35, strtok(gl{g}, ':'), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
end
set(ax, 'XTick', 1:L, 'XTickLabel', regexprep(LM.dims(ord), '^latent_', 'L'), 'XTickLabelRotation', 90, ...
        'YTick', 1:S, 'YTickLabel', LM.cols);
ylim(ax, [0 S + 0.5]);
title(ax, {'Which sensor does each latent dimension encode?', ...
           'R^2 of each latent dim from each sensor''s own trace PCs (linear); columns ordered by VIP group'}, ...
      'FontSize', 9.5);
i_style(ax);
end

%% =========================================================================
%  Helpers
%% =========================================================================
function f = i_newfig(name, sz)
f = figure('Color', 'w', 'Name', name, 'Units', 'pixels');
f.Position(3:4) = sz;
end

function i_style(ax)
set(ax, 'FontName', 'Helvetica', 'FontSize', 9, 'Box', 'off', 'TickDir', 'out', ...
        'LineWidth', 0.75, 'Color', 'none', 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
ax.YGrid = 'on';  ax.GridAlpha = 0.08;
end

function i_heat(ax, M, xl, yl, cm, cl, fmt)
% Annotated heatmap; text colour flips on dark cells
imagesc(ax, M, cl);  colormap(ax, cm);  cb = colorbar(ax);  cb.Box = 'off';  cb.TickDirection = 'out';
[r, c] = size(M);  mid = cl(1) + 0.6 * diff(cl);
for i = 1:r
    for j = 1:c
        if isnan(M(i,j)), continue; end
        tc = [0.15 0.15 0.15];  if M(i,j) > mid || M(i,j) < cl(1) + 0.15*diff(cl) && cl(1) < 0, tc = [1 1 1]; end
        text(ax, j, i, sprintf(fmt, M(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', tc);
    end
end
set(ax, 'XTick', 1:c, 'XTickLabel', xl, 'XTickLabelRotation', 45, 'YTick', 1:r, 'YTickLabel', yl, ...
        'TickLabelInterpreter', 'none', 'FontSize', 8, 'TickLength', [0 0], 'Box', 'off');
end

function cm = i_seq(n)
% White -> blue
cm = [linspace(1, 0.12, n)' linspace(1, 0.33, n)' linspace(1, 0.60, n)'];
end

function cm = i_div(n)
% Blue - white - red
h  = floor(n/2);
cm = [linspace(0.13,1,h)' linspace(0.40,1,h)' linspace(0.67,1,h)'; ...
      linspace(1,0.70,n-h)' linspace(1,0.09,n-h)' linspace(1,0.17,n-h)'];
end

function v = i_getf(s, f)
if isfield(s, f), v = s.(f); else, v = []; end
end

function i_save(f, o, ptype)
% <yyyy-mm-dd>__<plottype>_<tag>.svg ; overwrite on re-run
if ~o.save || isempty(o.savepath), return; end
if ~exist(o.savepath, 'dir'), mkdir(o.savepath); end
tg = '';  if ~isempty(o.tag), tg = ['_' matlab.lang.makeValidName(o.tag)]; end
fn = fullfile(o.savepath, sprintf('%s__%s%s.svg', char(datetime('today', 'Format', 'yyyy-MM-dd')), ptype, tg));
print(f, fn, '-dsvg', '-painters');
fprintf('  saved %s\n', fn);
end

%% =========================================================================
%  Poster figures (Q3): sources of predictable variance + prediction vs chance
%% =========================================================================
function [figs, P] = i_poster(R3, R3bg, o)
% (A) Partition of the Q2 of treatment + dynamics (commonality on held-out Q2):
%     dynamics only = Q2(tx+dyn) - Q2(tx); treatment only = Q2(tx+dyn) - Q2(dyn);
%     overlap = Q2(tx) + Q2(dyn) - Q2(tx+dyn). Negative components (CV noise /
%     suppression) are drawn as 0 and listed in the footnote.
% (B) Treatment vs dynamics (pooled, within well) against the shuffled-null band.
% Fixed physical size (cm), vector PDF + SVG. Colours: validated categorical slots
% (blue = dynamics, orange = treatment); the overlap is neutral grey.
% ---- palette (validated categorical slots 1-2, light surface) -------------------
C.dyn  = [ 42 120 214] / 255;            % blue   #2a78d6  dynamics
C.tx   = [235 104  52] / 255;            % orange #eb6834  treatment
C.ovl  = [185 183 174] / 255;            % neutral grey #b9b7ae  overlap
C.null = [228 227 222] / 255;            % recessive band #e4e3de  shuffled null
C.ink  = [ 11  11  11] / 255;            % primary text #0b0b0b
C.ink2 = [ 82  81  78] / 255;            % secondary text #52514e
C.grid = [232 231 226] / 255;            % grid #e8e7e2
C.surf = [1 1 1];

%% ---- Data ------------------------------------------------------------------
Q = R3.q2;
if ~isempty(o.poster_targets)
    [tf, loc] = ismember(string(o.poster_targets), Q.target);
    if any(~tf), warning('%s poster_targets not in q3: %s', '[dl_state_summary]', strjoin(string(o.poster_targets(~tf)), ', ')); end
    Q = Q(loc(tf), :);
end
nT = height(Q);
assert(nT > 0, '%s no poster targets to plot.', '[dl_state_summary]');
qw = nan(nT, 1);                                                   % within-well dynamics Q2
if isfield(R3, 'nl') && ismember('q2w_lin', R3.nl.Properties.VariableNames)
    [~, k] = ismember(Q.target, R3.nl.target);  qw(k > 0) = R3.nl.q2w_lin(k(k > 0));
end
qb = nan(nT, 1);                                                   % secondary-only dynamics Q2
if ~isempty(R3bg) && isfield(R3bg, 'q2')
    [~, k] = ismember(Q.target, R3bg.q2.target);  qb(k > 0) = R3bg.q2.q2_dyn(k(k > 0));
end
lab = i_pst_labels(Q.target, o.poster_labels);

% Partition (commonality on held-out Q2)
dOnly = Q.q2_txdyn - Q.q2_tx;
tOnly = Q.q2_txdyn - Q.q2_dyn;
ovl   = Q.q2_tx + Q.q2_dyn - Q.q2_txdyn;
P = table(Q.target, dOnly, ovl, tOnly, Q.q2_txdyn, any([dOnly ovl tOnly] < 0, 2), ...
          'VariableNames', {'target','dyn_only','overlap','tx_only','total','clipped'});
figs = gobjects(0);
sig = Q.p_delta_holm < o.poster_alpha & Q.dq2 > 0;                        % dynamics adds beyond treatment

%% ---- (A) Partition: stacked horizontal bars ----------------------------------------
if any(strcmpi(o.poster, {'both','partition'}))
    ord = (1:nT)';  if o.poster_sort, [~, ord] = sort(P.total, 'descend'); end
    Sg  = max([P.dyn_only P.overlap P.tx_only], 0);  Sg = Sg(ord, :);
    hgt = 2.2 + 0.62 * nT + 2.2;                                    % title/legend + rows + axis/footnote (cm)
    f = i_pst_fig(o, hgt);  ax = i_pst_axes(f, o, C, hgt);
    y = 1:nT;
    hb = barh(ax, y, Sg, 0.62, 'stacked', 'EdgeColor', C.surf, 'LineWidth', 1.5);   % surface gaps between segments
    hb(1).FaceColor = C.dyn;  hb(2).FaceColor = C.ovl;  hb(3).FaceColor = C.tx;
    xmax = max([sum(Sg, 2); P.total; qb(ord)], [], 'omitnan');  xpad = 0.22 * xmax;
    for i = 1:nT                                                     % total, direct label
        t = P.total(ord(i));
        s = sprintf('%.3f', t);  if sig(ord(i)), s = [s '*']; end %#ok<AGROW>
        text(ax, max(sum(Sg(i, :)), t) + 0.012 * xmax, y(i), s, 'Color', C.ink2, 'FontSize', o.poster_font_size - 1, ...
             'FontName', o.poster_font, 'VerticalAlignment', 'middle');
    end
    hl = gobjects(1, 3);  hl(1:3) = hb;  lg = {'dynamics only', 'either (overlap)', 'treatment only'};
    if any(~isnan(qb))
        qq = qb(ord);  ok = ~isnan(qq);
        hl(4) = plot(ax, qq(ok), y(ok), '|', 'Color', C.ink, 'MarkerSize', 13, 'LineWidth', 1.6);
        lg{4} = 'technical floor (secondary-only)';
    end
    set(ax, 'YTick', y, 'YTickLabel', lab(ord), 'YDir', 'reverse', 'YLim', [0.4 nT + 0.6], ...
            'XLim', [0, xmax + xpad]);
    xlabel(ax, 'Q^2  (fraction of marker variance predicted, held-out wells)', 'Color', C.ink2);
    title(ax, 'Sources of predictable endpoint variance', 'FontWeight', 'bold', 'Color', C.ink, ...
          'FontSize', o.poster_font_size + 1);
    legend(ax, hl, lg, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
           'FontSize', o.poster_font_size - 1, 'TextColor', C.ink2, 'NumColumns', min(numel(lg), 2));
    fn = {'* dynamics add beyond treatment (Holm p < 0.05). Number = Q^2 of treatment + dynamics.'};
    if any(P.clipped)
        fn{end+1} = sprintf('Negative components drawn as 0 (max %.3f): %s.', ...
                            -min(min([P.dyn_only P.overlap P.tx_only], [], 2)), strjoin(lab(P.clipped), ', '));
    end
    i_pst_foot(f, fn, o, C);
    figs(end+1) = f;
    i_pst_save(f, o, 'q3-poster-partition');
end

%% ---- (B) Comparison: dot plot against the null band ----------------------------------
if any(strcmpi(o.poster, {'both','comparison'}))
    ord = (1:nT)';  if o.poster_sort, [~, ord] = sort(Q.q2_dyn, 'descend'); end
    hgt = 2.2 + 0.62 * nT + 2.2;
    f = i_pst_fig(o, hgt);  ax = i_pst_axes(f, o, C, hgt);
    y = (1:nT)';
    vals = [Q.q2_tx Q.q2_dyn qw Q.q2_null95];  vals = vals(ord, :);
    xl = [min(0, min(vals(:), [], 'omitnan')) - 0.005, max(vals(:), [], 'omitnan') * 1.08 + 0.005];
    for i = 1:nT                                                     % null band, per row
        nb = max(vals(i, 4), 0);
        patch(ax, [0 nb nb 0], y(i) + [-0.3 -0.3 0.3 0.3], C.null, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hN = patch(ax, nan(1, 4), nan(1, 4), C.null, 'EdgeColor', 'none');
    xline(ax, 0, '-', 'Color', C.ink2, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    hT = plot(ax, vals(:, 1), y, 'o', 'MarkerSize', 8, 'MarkerFaceColor', C.tx,  'MarkerEdgeColor', C.surf, 'LineWidth', 1);
    hD = plot(ax, vals(:, 2), y, 'o', 'MarkerSize', 9, 'MarkerFaceColor', C.dyn, 'MarkerEdgeColor', C.surf, 'LineWidth', 1);
    hl = gobjects(1, 0);  hl(1) = hT;  hl(2) = hD;  lg = {'treatment', 'dynamics'};
    if any(~isnan(qw))
        hl(end+1) = plot(ax, vals(:, 3), y, 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', C.dyn, 'LineWidth', 1.4);
        lg{end+1} = 'dynamics, within well';
    end
    hl(end+1) = hN;  lg{end+1} = 'shuffled null (95%)';
    set(ax, 'YTick', y, 'YTickLabel', lab(ord), 'YDir', 'reverse', 'YLim', [0.4 nT + 0.6], 'XLim', xl);
    xlabel(ax, 'Q^2  (fraction of marker variance predicted, held-out wells)', 'Color', C.ink2);
    title(ax, 'Endpoint prediction vs chance', 'FontWeight', 'bold', 'Color', C.ink, ...
          'FontSize', o.poster_font_size + 1);
    legend(ax, hl, lg, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', ...
           'FontSize', o.poster_font_size - 1, 'TextColor', C.ink2, 'NumColumns', 2);
    i_pst_foot(f, {'Within well: observed and predicted centred on their well mean (removes treatment and plating).'}, o, C);
    figs(end+1) = f;
    i_pst_save(f, o, 'q3-poster-comparison');
end
end

function f = i_pst_fig(o, hgt)
% Figure at its printed size (cm), white surface
f = figure('Color', 'w', 'Units', 'centimeters', 'Name', 'Q3 poster');
f.Position(3:4) = [o.poster_width_cm, hgt];
f.PaperUnits = 'centimeters';  f.PaperSize = [o.poster_width_cm, hgt];  f.PaperPosition = [0 0 o.poster_width_cm hgt];
end

function ax = i_pst_axes(f, o, C, hgt)
% Recessive axes: x grid only, no box, outward ticks, one font. Placed in cm:
% 2.8 cm for row labels, 2.2 cm below (axis label + footnote), 2.2 cm above (title + legend)
ax = axes(f, 'Units', 'centimeters', 'Position', [2.8, 2.2, o.poster_width_cm - 2.8 - 0.8, hgt - 4.4]);
hold(ax, 'on');
set(ax, 'FontName', o.poster_font, 'FontSize', o.poster_font_size, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.008 0.008], ...
        'XGrid', 'on', 'YGrid', 'off', 'GridColor', C.grid, 'GridAlpha', 1, 'XColor', C.ink2, 'YColor', C.ink2, ...
        'TickLabelInterpreter', 'tex', 'Layer', 'bottom', 'LineWidth', 0.8);
end

function i_pst_foot(f, lines, o, C)
% Footnote under the axes, secondary ink
annotation(f, 'textbox', [0.02 0.0 0.96 0.9 / f.Position(4)], 'String', lines, 'EdgeColor', 'none', ...
           'FontName', o.poster_font, 'FontSize', o.poster_font_size - 2, 'Color', C.ink2, 'Interpreter', 'tex', ...
           'VerticalAlignment', 'bottom', 'FitBoxToText', 'off');
end

function lab = i_pst_labels(names, over)
% Display names for markers (TeX); user overrides win
names = cellstr(names);
def = {'HLA_DR','HLA-DR'; 'PDL1','PD-L1'; 'IL1b','IL-1\beta'; 'TNFa','TNF\alpha'; 'IL8','IL-8'; ...
       'IL10RA','IL-10RA'; 'PSTAT3','pSTAT3'; 'PSTAT6','pSTAT6'; 'CD206','CD206'; 'CD68','CD68'; ...
       'CD86','CD86'; 'iNOS','iNOS'};
if ~isempty(over), def = [over; def]; end
lab = strrep(names, '_', '\_');
for i = 1:numel(names)
    k = find(strcmp(def(:, 1), names{i}), 1);
    if ~isempty(k), lab{i} = def{k, 2}; end
end
end

function i_pst_save(f, o, ptype)
% <yyyy-mm-dd>__<ptype>_<tag>.pdf + .svg (vector) ; overwrite on re-run
if ~o.save || isempty(o.savepath), return; end
if ~exist(o.savepath, 'dir'), mkdir(o.savepath); end
tg = '';  if ~isempty(o.tag), tg = ['_' matlab.lang.makeValidName(o.tag)]; end
base = fullfile(o.savepath, sprintf('%s__%s%s', char(datetime('today', 'Format', 'yyyy-MM-dd')), ptype, tg));
exportgraphics(f, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
print(f, [base '.svg'], '-dsvg', '-painters');
fprintf('  saved %s.pdf / .svg\n', base);
end