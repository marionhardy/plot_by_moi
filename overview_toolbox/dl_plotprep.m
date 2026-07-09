function prep = dl_plotprep(plotby, dataloc, p)
% DL_PLOTPREP  Shared data-prep core for the v5 plot suite.
%
% prep = dl_plotprep(plotby, dataloc, p)
%
% Consolidates ALL platemap-parsing, condition-indexing, replicate-
% combining, and time-window resolution that v4's plot_by_MHY did inline
% (and that the first v5 port duplicated across three files). The three
% renderers (plot_mean, plot_heatmap, plot_sorted_stacks) call this once
% and consume the returned struct -- they contain NO platemap logic.
%
% INPUT
%   plotby  : 'treatment' | 'celltype' | 'custom'
%   dataloc : v5 dataloc struct (also accepts v4/legacy; see dl_check_*).
%   p       : parsed parameter struct from the calling renderer. Fields
%             used here: channel, combinexys, subset, exclude, specsubset,
%             nogene, tx_order, ifchan, aftertreatment, tmaxaftertx,
%             tbeforetx, plotfromzero, ncells, tktm.
%
% OUTPUT  prep struct with fields:
%   .dataloc    : dataloc, possibly with .d/.z/.IFd pooled if combinexys.
%   .facetby    : cellstr of facet-level fieldnames (one FIGURE each).
%   .groupby    : cellstr of group fieldnames (tiles within a figure).
%   .legname    : cellstr, display label per group (spaces, not underscores).
%   .idx        : struct, condition-name -> xy-index vector (post-combine).
%   .goodxy     : logical vector, xys carrying plottable data.
%   .linetp     : cell, per-xy treatment-timepoint vector (for xlines).
%   .xymat      : cell, per-xy index vector (post-combine numbering).
%   .channel    : validated channel spec (cellstr, 1+ for overlay).
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): this is the single biggest
% structural change of the refactor. v4 interleaved this logic with
% rendering across a 560-line nested switch; the first v5 port copied it
% verbatim into each plot_ file, which is where the stridx-3D and
% figure-per-facet regressions crept in. Centralizing here means the
% parsing is written once, tested once (see dl_plotprep diagnostics in the
% companion Obsidian note), and the renderers stay thin. Named pmd layer
% access via dl_pmd_layers (crit #4); replicate combine via dl_combine_xys
% (crit #6). NOT backwards-mutating: operates on a local copy of dataloc.
% -------------------------------------------------------------------------

if isempty(dataloc.platemapd.pmd)
    error('dl_plotprep:NoPlatemap', ...
        ['dataloc has no platemap data. Run dl_ensure_platemap first, ' ...
         'and dl_check_platemap_compat to verify NA sanitization.']);
end

%% ---- 1. Cell-line names per well (Cell x Gene concatenation) ----
pmd = dataloc.platemapd.pmd;
nct = ceil(size(pmd.Cell,3)/3);
ci  = 3*((1:nct)-1) + 1;                 % Cell-name layers (every 3rd)

if ~p.nogene
    genecat = cell(size(pmd.Gene(:,:,1)));
    for s = 1:size(pmd.Gene,3)
        genecat = cellfun(@(x,y)[x,y], genecat, pmd.Gene(:,:,s), 'Un', 0);
    end
    cnames = cellfun(@(x,y)[x,'_',y], pmd.Cell(:,:,ci), ...
        repmat(genecat,[1,1,numel(ci)]), 'Un', 0);
else
    cnames = cell(size(pmd.Cell(:,:,1)));
    for s = 1:numel(ci)
        cnames = cellfun(@(x,y)[x,y], cnames, pmd.Cell(:,:,s), 'Un', 0);
    end
    nancell = find(~cell2mat(cellfun(@ischar,cnames,'Un',0)));
    for s = 1:numel(nancell); cnames{nancell(s)} = ''; end
end
cnames = regexprep(cnames, {'\W','^[\d_]*(\w)'}, {'','$1'});

celltypes = unique(cnames(cellfun(@ischar,cnames)));
gn = cellfun(@isvarname, celltypes);
celltypes = celltypes(gn);
celltypes = celltypes(~cellfun(@isempty, celltypes));

cellfn = cellfun(@(x)x{1}, regexp(celltypes,'@','split'), 'un', 0);
cellfn = arrayfun(@(x)regexprep(cellfn{x},'[^a-zA-Z0-9]','_'), 1:numel(cellfn), 'un', 0);
cellfn = arrayfun(@(x)regexprep(cellfn{x},'(^[\d_]+\w)','x$1'), 1:numel(cellfn), 'un', 0)';

idx = struct();
for s = 1:numel(celltypes)
    idx.(celltypes{s}) = [pmd.xy{strcmp(celltypes{s},cnames)}];
end

%% ---- 2. Treatment names per well (named layer access, crit #4) ----
L = dl_pmd_layers();
catTxnames = intersect(L.txFields, fieldnames(pmd), 'stable');

Txcat  = cell(size(pmd.Cell,1), size(pmd.Cell,2));   % 8x12 seed
linetp = cell(size(pmd.Cell,1), size(pmd.Cell,2));

for s = 1:numel(catTxnames)
    nsub = L.nsub(size(pmd.(catTxnames{s}),3));
    subTxCat = [];
    for iSub = 1:nsub
        tid = false(1, nsub*L.block);
        tid((iSub-1)*L.block + L.nameParts) = true;      % Name+Conc+Units
        subSlice = cellfun(@(x) local_concat(x), ...
            cat(3, num2cell(pmd.(catTxnames{s})(:,:,tid), 3)), 'Un', 0);
        if isempty(subTxCat); subTxCat = subSlice;
        else; subTxCat = cellfun(@(a,b) local_join(a,b), subTxCat, subSlice, 'Un', 0);
        end
    end
    Txcat  = cat(3, Txcat, subTxCat);
    linetp = cat(3, linetp, cellfun(@(x)cat(2,x{:}), ...
        cat(3, num2cell(pmd.(catTxnames{s})(:,:,L.idx(1,'Time')), 3)), 'Un', 0));
end
linetp = num2cell(linetp(:,:,2:end), 3);

for s = 1:numel(Txcat); if strcmp(Txcat{s},''); Txcat{s} = NaN; end; end
stridx  = cellfun(@ischar, Txcat);                       % 3-D here
spacer  = repmat({'+'}, size(Txcat)); spacer(~stridx) = {[]};

% Flatten sub-treatment layers into one label per well
Txcat2 = Txcat(:,:,2);
for s = 2:size(Txcat,3)-1
    Txcat2 = cat(3, Txcat2, spacer(:,:,s+1), Txcat(:,:,s+1));
end
Txcat2 = cellfun(@(x)cat(2,x{:}), cat(3,num2cell(Txcat2,3)), 'Un', 0);
stridx = cellfun(@ischar, Txcat2);                       % 2-D: matches pmd.xy

txs        = unique(Txcat2(stridx));
treatments = unique(Txcat2(stridx));
for s = 1:numel(txs)
    txs{s} = regexprep(txs{s}, char(0), '');
    txs{s} = regexprep(txs{s}, '(^|\s|+)', '');
    txs{s} = regexprep(txs{s}, '[^a-zA-Z0-9]', '_');
end
for s = 1:numel(Txcat2)
    if ischar(Txcat2{s})
        Txcat2{s} = regexprep(Txcat2{s}, char(0), '');
        Txcat2{s} = regexprep(Txcat2{s}, '(^|\s|+)', '');
        Txcat2{s} = regexprep(Txcat2{s}, '[^a-zA-Z0-9]', '_');
    end
end

%% ---- 3. Subset / exclude filtering ----
[txs, treatments, cellfn] = local_filter(txs, treatments, cellfn, p.subset,  p.specsubset, false);
[txs, treatments, cellfn] = local_filter(txs, treatments, cellfn, p.exclude, false,        true);

%% ---- 4. tx_order sort (order legend by a chosen treatment column) ----
splitTxs = cellfun(@(x)strsplit(x,'+'), treatments, 'Un', 0);
if ~isempty(splitTxs)
    sz = max(cellfun(@numel, splitTxs));
    ss = strings(numel(treatments), sz);
    for kk = 1:numel(treatments)
        for jj = 1:numel(splitTxs{kk}); ss(kk,jj) = splitTxs{kk}{jj}; end
    end
    if numel(splitTxs) > 2
        [~,I] = sort(ss);
        txs = txs(I(:,p.tx_order)); treatments = treatments(I(:,p.tx_order));
    end
end
treatments = arrayfun(@(x)regexprep(treatments{x},char(0),''), 1:numel(treatments), 'un', 0)';
treatments = arrayfun(@(x)regexprep(treatments{x},'(^|\s)',''), 1:numel(treatments), 'un', 0)';
treatments = arrayfun(@(x)regexprep(treatments{x},'++','+'),    1:numel(treatments), 'un', 0)';
treatments = arrayfun(@(x)regexprep(treatments{x},'_',' '),     1:numel(treatments), 'un', 0)';

xymat = pmd.xy(stridx);
for s = 1:numel(txs)
    idx.(txs{s}) = cat(2, xymat{strcmp(txs{s}, Txcat2(stridx))});
end
linetp = linetp(stridx);

%% ---- 5. goodxy: which xys carry data (fabricate cellindex + warn) ----
goodxy = false(1, max([pmd.xy{:}]));
for ii = 1:numel(dataloc.d)
    d = dataloc.d{ii};
    if isempty(d) || ~isfield(d,'data') || isempty(fieldnames(d.data)); continue; end
    if ~isfield(d,'cellindex') || isempty(d.cellindex)
        fn = fieldnames(d.data);
        nRows = size(d.data.(fn{1}),1);
        dataloc.d{ii}.cellindex = (1:nRows)';
        warning('dl_plotprep:FabricatedCellIndex', ...
            'XY %d has no cellindex (legacy dataloc?) -- fabricated 1:%d.', ii, nRows);
    end
    goodxy(ii) = true;
end

%% ---- 6. Combine replicate xys if requested (crit #6) ----
if p.combinexys
    [dataloc, idx, linetp, xymat, goodxy] = ...
        dl_combine_xys(dataloc, txs, cellfn, idx, xymat, linetp, goodxy, p);
end

%% ---- 7. Facet / group assignment ----
switch plotby
    case 'treatment'; facetby = cellfn; groupby = txs;    legname = treatments;
    case 'celltype';  facetby = txs;    groupby = cellfn; legname = cellfn;
    case 'custom';    facetby = p.facetby; groupby = p.groupby; legname = cellfn;
end
legname = arrayfun(@(x)regexprep(legname{x},'_',' '), 1:numel(legname), 'un', 0)';

%% ---- 8. Pack ----
prep = struct();
prep.dataloc = dataloc;
prep.facetby = facetby;
prep.groupby = groupby;
prep.legname = legname;
prep.idx     = idx;
prep.goodxy  = goodxy;
prep.linetp  = linetp;
prep.xymat   = xymat;
prep.channel = p.channel;
end % dl_plotprep


% ===================== local helpers =====================

function out = local_concat(x)
% Concatenate the char parts (Name+Conc+Units) of one well's layer cell.
chars = x(cellfun(@ischar, x));
if isempty(chars); out = ''; else; out = strtrim(cat(2, chars{:})); end
end

function out = local_join(a, b)
% Join two sub-treatment labels with '+' (skipping empties).
a = strtrim(a); b = strtrim(b);
if isempty(a) && isempty(b); out = '';
elseif isempty(a); out = b;
elseif isempty(b); out = a;
else; out = [a '+' b];
end
end

function [txs, treatments, cellfn] = local_filter(txs, treatments, cellfn, terms, specAll, isExclude)
% Apply subset (keep) or exclude (drop) term matching to tx & cell lists.
if isempty(terms); return; end
if ~iscell(terms); terms = {terms}; end
txHit = false(numel(txs),   numel(terms));
ctHit = false(numel(cellfn),numel(terms));
for k = 1:numel(terms)
    txHit(:,k) = contains(txs,    terms{k}, 'IgnoreCase', true);
    ctHit(:,k) = contains(cellfn, terms{k}, 'IgnoreCase', true);
end
if specAll; txHit = all(txHit,2); ctHit = all(ctHit,2);
else;       txHit = any(txHit,2); ctHit = any(ctHit,2);
end
if isExclude
    if any(txHit); txs = txs(~txHit); treatments = treatments(~txHit); end
    if any(ctHit); cellfn = cellfn(~ctHit); end
else
    if any(txHit); txs = txs(txHit);  treatments = treatments(txHit);  end
    if any(ctHit); cellfn = cellfn(ctHit); end
end
end
