function ms_color_embedding(S, colorBy, varargin)
%CT_COLOR_EMBEDDING  Colour the latent embedding (PC or UMAP) by any platemap
%                    field, CyCIF marker, or cluster label; optionally FACET
%                    into panels by a grouping column.
%
%   Uses the row-aligned object from dl_build_live_cycif_tsc: S.cellT (spine +
%   PC1/PC2 [+umap1/2]), S.cycifT (markers). Any per-cell vector of length
%   nCells can colour the scatter; any categorical cellT column can facet it.
%
% USAGE:
%   ms_color_embedding(S, 'condition')
%   ms_color_embedding(S, 'CD206','Log2',true)
%   ms_color_embedding(S, clustID)                              % numeric labels
%   ms_color_embedding(S, 'CD206','Log2',true,'FacetBy','ptx_pool_label')
%   ms_color_embedding(S, 'tx1_label','FacetBy','ptx_pool_label','Space','umap')
%
% INPUT:
%   S        struct from dl_build_live_cycif_tsc
%   colorBy  char/string cellT column (categorical) | cycifT marker
%            (continuous) | numeric [nCells x 1] vector (labels)
%
% NAME-VALUE:
%   'Space'   'pc' (default) | 'umap'
%   'Log2'    log2(max(x,1)) a continuous marker before colouring [false]
%   'FacetBy' cellT column to split into panels ['' = single plot]
%   'PtSize'  marker size [6]
%   'Title'   custom super-title ['']

q.Space='pc'; q.Log2=false; q.FacetBy=''; q.PtSize=6; q.Title='';
k=1; while k<=numel(varargin)
    switch lower(varargin{k})
        case 'space',   q.Space=lower(varargin{k+1}); k=k+2;
        case 'log2',    q.Log2=logical(varargin{k+1}); k=k+2;
        case 'facetby', q.FacetBy=char(varargin{k+1}); k=k+2;
        case 'ptsize',  q.PtSize=varargin{k+1}; k=k+2;
        case 'title',   q.Title=varargin{k+1}; k=k+2;
        otherwise, error('unknown param "%s"',varargin{k});
    end
end

% --- embedding coordinates (Space names the column pair, e.g. pc, umap, cyumap, cypc) ---
sp = lower(q.Space);
% map friendly names to column stems
switch sp
    case {'pc','pca'},     c1='PC1';     c2='PC2';      % latent PCA -> PC1/PC2
    case 'umap',           c1='umap1';   c2='umap2';
    case 'cyumap',         c1='cyumap1'; c2='cyumap2';   % CyCIF-marker UMAP
    case {'cypc','cypca'}, c1='cyPC1';   c2='cyPC2';     % CyCIF-marker PCA
    case 'tsne',           c1='tsne1';   c2='tsne2';
    otherwise      % treat Space as a literal prefix: <sp>1 / <sp>2
        c1=[q.Space '1']; c2=[q.Space '2'];
end
assert(all(ismember({c1,c2},S.cellT.Properties.VariableNames)), ...
    'embedding columns %s/%s not in cellT — run ms_embedding first.', c1, c2);
XY=[S.cellT.(c1) S.cellT.(c2)]; axlbl={c1,c2};
n=size(XY,1);

% --- colour vector + kind ---
if isnumeric(colorBy)
    assert(numel(colorBy)==n,'numeric colorBy must be [nCells x 1]');
    cvec=colorBy(:); kind='label'; cname='label';
elseif ismember(char(colorBy),S.cellT.Properties.VariableNames)
    cvec=S.cellT.(char(colorBy)); cname=char(colorBy); kind='cat';
elseif ismember(char(colorBy),S.cycifT.Properties.VariableNames)
    cvec=S.cycifT.(char(colorBy)); cname=char(colorBy);
    if q.Log2, cvec=log2(max(cvec,1)); cname=[cname ' (log2)']; end
    kind='cont';
else
    error('colorBy "%s" not found in cellT or cycifT, and is not numeric.',char(colorBy));
end

% --- facet levels ---
use_facet = ~isempty(q.FacetBy);
if use_facet
    assert(ismember(q.FacetBy,S.cellT.Properties.VariableNames), ...
        'FacetBy "%s" not a cellT column.', q.FacetBy);
    fvar = S.cellT.(q.FacetBy);
    flev = unique(fvar,'stable');
    nF = numel(flev);
else
    nF = 1;
end

% --- shared square limits across all panels (so facets are comparable) ---
pad = 0.03*range(XY,1);
lim = [min(XY)-pad; max(XY)+pad];   % 2x2: rows=[min;max], cols=[x y]

% --- colour scale shared across panels for continuous ---
if strcmp(kind,'cont')
    clim_shared = prctile(cvec(isfinite(cvec)),[1 99]);   % robust, shared
end

figure('Color','w','Name',sprintf('%s by %s%s', upper(q.Space), cname, ...
    ternary(use_facet,[' | facet ' q.FacetBy],'')));
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');

for fi = 1:nF
    if use_facet
        sel = fvar==flev(fi); ptitle = string(flev(fi));
    else
        sel = true(n,1); ptitle = cname;
    end
    ax = nexttile; hold(ax,'on');
    switch kind
        case 'cont'
            scatter(ax, XY(sel,1),XY(sel,2), q.PtSize, cvec(sel), 'filled');
            colormap(ax,parula); caxis(ax, clim_shared);
        case {'cat','label'}
            g = categorical(cvec(sel));
            gscatter(ax, XY(sel,1),XY(sel,2), g, [], '.', q.PtSize);
            if fi==1
                lg=legend(ax,'Location','bestoutside','Interpreter','none');
            elseif ~isempty(findobj(ax,'Type','Legend'))
                legend(ax,'off');   % one legend suffices for shared categories
            end
    end
    xlim(ax,lim(:,1)'); ylim(ax,lim(:,2)'); axis(ax,'square');
    xlabel(ax,axlbl{1}); ylabel(ax,axlbl{2});
    title(ax, ptitle, 'Interpreter','none');
end

if strcmp(kind,'cont')
    cb=colorbar(ax); cb.Layout.Tile='east'; cb.Label.String=cname;   % one shared colorbar
end
if isempty(q.Title)
    q.Title=sprintf('Latent %s coloured by %s%s', upper(q.Space), cname, ...
        ternary(use_facet,[' (faceted by ' q.FacetBy ')'],''));
end
title(tl, q.Title, 'Interpreter','none');
end

% inline ternary helper
function out = ternary(cond,a,b); if cond, out=a; else, out=b; end; end
