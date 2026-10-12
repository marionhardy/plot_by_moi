function ms_pc_histograms(S, clustCol, varargin)
%MS_PC_HISTOGRAMS  Per-axis marginal histograms of an embedding, one coloured
%   histogram per cluster. Shows how clusters separate/overlap ALONG each axis
%   (the 1-D marginals of the 2-D scatter). Version-safe: uses histcounts +
%   manual line drawing (no histogram() — avoids MATLAB string-arg quirks).
%
% USAGE:
%   ms_pc_histograms(S, 'clust_latent');
%   ms_pc_histograms(S, 'clust_latent', 'Space','pca', 'Layout','joint');
%   ms_pc_histograms(S, 'clust_latent', 'Space','umap', 'Layout','joint');   % titled UMAP
%   ms_pc_histograms(S, clustID_vector, 'Bins',60, 'Style','bar');
%
% INPUT:
%   S         object from dl_build_live_cycif_tsc (needs the embedding cols)
%   clustCol  cellT column name (or [nCells x 1] label vector)
%
% NAME-VALUE:
%   'Space'   'pc'|'pca'|'umap'|'cyumap'|... which embedding    ['pc']
%   'Bins'    number of bins                                    [50]
%   'Norm'    'count' | 'pdf' | 'probability'                   ['probability']
%   'DropNoise' drop label 0/""/"noise"                         [true]
%   'Style'   'bar' (translucent) | 'step' (outline)            ['step']
%   'Layout'  'tiled' (two panels) | 'sep' (two figures) |
%             'joint' (scatter + PC1 top marginal + PC2 right marginal) ['tiled']
%
% CHANGELOG
%   2026-10-12 (v5): the 'joint' figure name and title said "PCA" for every Space, so the
%   UMAP pass was labelled as PCA (the data drawn were correct). They now name the space
%   actually drawn. Label-only change; interface and outputs unchanged.

q.Space='pc'; q.Bins=50; q.Norm='probability'; q.DropNoise=true;
q.Style='fill'; q.Layout='tiled';
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'space',     q.Space=lower(varargin{k+1}); k=k+2;
        case 'bins',      q.Bins=varargin{k+1}; k=k+2;
        case 'norm',      q.Norm=lower(varargin{k+1}); k=k+2;
        case 'dropnoise', q.DropNoise=logical(varargin{k+1}); k=k+2;
        case 'style',     q.Style=lower(varargin{k+1}); k=k+2;
        case 'layout',    q.Layout=lower(varargin{k+1}); k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

% resolve the two embedding axis columns (same map as ms_color_embedding)
switch q.Space
    case {'pc','pca'},     c1='PC1';     c2='PC2';     spName='PCA';
    case 'umap',           c1='umap1';   c2='umap2';   spName='UMAP';
    case 'cyumap',         c1='cyumap1'; c2='cyumap2'; spName='CyCIF UMAP';
    case {'cypc','cypca'}, c1='cyPC1';   c2='cyPC2';   spName='CyCIF PCA';
    case 'tsne',           c1='tsne1';   c2='tsne2';   spName='t-SNE';
    otherwise,             c1=[q.Space '1']; c2=[q.Space '2']; spName=upper(q.Space);
end
assert(all(ismember({c1,c2},S.cellT.Properties.VariableNames)), ...
    'embedding columns %s/%s not in cellT — run ms_embedding first.', c1, c2);
AX = {S.cellT.(c1), S.cellT.(c2)};  axname = {c1,c2};

% labels
if isnumeric(clustCol)||iscategorical(clustCol), lab=clustCol; else, lab=S.cellT.(char(clustCol)); end
labs = string(lab);
cl = unique(labs);
if q.DropNoise, cl(cl=="0"|cl==""|lower(cl)=="noise")=[]; end
nClust = numel(cl);
cols = lines(nClust);

% shared bin edges per axis
edges = cell(1,2);
for a=1:2
    x = AX{a}(isfinite(AX{a}));
    edges{a} = linspace(min(x), max(x), q.Bins+1);
end

titleStr = char(string(clustCol));

if strcmp(q.Layout,'joint')
    % scatter centre + axis-2 marginal (LEFT) + axis-1 marginal (BOTTOM), filled.
    figure('Color','w','Name',[spName ' joint by ' titleStr]);
    tl = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');

    % LEFT marginal (axis 2): rows 1-3, col 1 — value on Y, density on X (flipped axis)
    axLt = nexttile(tl,1,[3 1]);
    local_draw_marg(axLt, AX{2}, labs, cl, nClust, cols, edges{2}, q.Norm, true, 'fill');
    set(axLt,'XDir','reverse');                    % density grows toward the scatter
    ylabel(axLt,axname{2}); xlabel(axLt,q.Norm);
    ylim(axLt,[edges{2}(1) edges{2}(end)]);

    % scatter centre: rows 1-3, cols 2-4
    axSc = nexttile(tl,2,[3 3]); hold(axSc,'on');
    for i=1:nClust
        m = labs==cl(i);
        scatter(axSc, AX{1}(m), AX{2}(m), 6, cols(i,:), 'filled', ...
            'MarkerFaceAlpha',0.5, 'DisplayName',sprintf('c%s',cl(i)));
    end
    set(axSc,'YTickLabel',[]);                     % axis 2 labelled on the left marginal
    xlabel(axSc,axname{1});
    xlim(axSc,[edges{1}(1) edges{1}(end)]); ylim(axSc,[edges{2}(1) edges{2}(end)]);
    legend(axSc,'Location','bestoutside','Interpreter','none');

    % BOTTOM marginal (axis 1): row 4, cols 2-4
    axBot = nexttile(tl,14,[1 3]);
    local_draw_marg(axBot, AX{1}, labs, cl, nClust, cols, edges{1}, q.Norm, false, 'fill');
    set(axBot,'YDir','reverse');                   % density grows toward the scatter
    xlabel(axBot,axname{1}); ylabel(axBot,q.Norm);
    xlim(axBot,[edges{1}(1) edges{1}(end)]);

    title(tl, sprintf('%s + marginals (%s)', spName, titleStr), 'Interpreter','none');
    return
end

% tiled / sep: two marginal panels
if strcmp(q.Layout,'tiled')
    figure('Color','w','Name',[spName ' marginals by ' titleStr]);
    tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
end
for a=1:2
    if strcmp(q.Layout,'tiled'), ax=nexttile(tl);
    else, figure('Color','w','Name',[axname{a} ' by ' titleStr]); ax=gca; end
    local_draw_marg(ax, AX{a}, labs, cl, nClust, cols, edges{a}, q.Norm, false, q.Style);
    xlabel(ax,axname{a}); ylabel(ax,q.Norm);
    title(ax,[axname{a} ' distribution per cluster']);
    if a==1, legend(ax,'Location','best','Interpreter','none'); end
end
if strcmp(q.Layout,'tiled')
    title(tl, sprintf('Per-axis cluster marginals, %s (%s)', spName, titleStr), 'Interpreter','none');
end
end

% ---- local drawing helper (histcounts + manual line; no histogram()) -----
function local_draw_marg(ax, xall, labs, cl, nClust, cols, ed, norm, horiz, style)
    if nargin<10, style='step'; end
    hold(ax,'on');
    ctrs = ed(1:end-1) + diff(ed)/2;
    for i=1:nClust
        x = xall(labs==cl(i)); x = x(isfinite(x));
        switch norm
            case 'count',       cnts = histcounts(x, ed);
            case 'probability', cnts = histcounts(x, ed, 'Normalization','probability');
            case 'pdf',         cnts = histcounts(x, ed, 'Normalization','pdf');
            otherwise,          cnts = histcounts(x, ed);
        end
        dn = sprintf('c%s (n=%d)', cl(i), numel(x));
        % staircase outline coords (value axis vs density axis)
        vx = [ed(1), repelem(ed(2:end-1),2), ed(end)];   % along the value axis
        vy = repelem(cnts, 2);                            % density
        switch style
            case 'fill'   % translucent filled staircase, closed to baseline
                if horiz   % value on Y, density on X (left/right marginal)
                    px = [0, vy, 0];  py = [ed(1), vx, ed(end)];
                else       % value on X, density on Y (top/bottom marginal)
                    px = [ed(1), vx, ed(end)];  py = [0, vy, 0];
                end
                patch(ax, px, py, cols(i,:), 'FaceAlpha',0.3, ...
                    'EdgeColor',cols(i,:), 'LineWidth',1, 'DisplayName',dn);
            case 'bar'
                bar(ax, ctrs, cnts, 1, 'FaceColor',cols(i,:), 'FaceAlpha',0.35, ...
                    'EdgeColor','none', 'DisplayName',dn);
            otherwise      % 'step' outline
                if horiz
                    plot(ax, vy, vx, '-', 'Color',cols(i,:), 'LineWidth',1.3, 'DisplayName',dn);
                else
                    plot(ax, vx, vy, '-', 'Color',cols(i,:), 'LineWidth',1.3, 'DisplayName',dn);
                end
        end
    end
end