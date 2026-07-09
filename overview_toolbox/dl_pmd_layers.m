function L = dl_pmd_layers()
% DL_PMD_LAYERS  Named layer-index map for platemapd.pmd treatment fields.
%
% OUTPUT
%   L : struct of named indices into the 3rd dim of a pmd treatment field.
%       Each treatment field (pTx, Tx1..Tx4) stores sub-treatments as
%       consecutive 5-layer blocks: [Name, Conc, Units, Time, TimeUnit].
%       A field with K sub-treatments has size(...,3) == 5*K.
%
%       L.block      = 5                      % layers per sub-treatment
%       L.Name       = 1                      % offset within a block
%       L.Conc       = 2
%       L.Units      = 3
%       L.Time       = 4
%       L.TimeUnit   = 5
%       L.nameParts  = [1 2 3]                % Name+Conc+Units (label text)
%       L.txFields   = {'pTx','Tx1','Tx2','Tx3','Tx4'}
%
%   Helper closures (avoid magic arithmetic at call sites):
%       L.idx(k, which)  -> absolute layer index for sub-treatment k,
%                           where `which` is one of the field names above
%                           (e.g. L.idx(2,'Time') -> 9 for the 2nd block).
%       L.nsub(field3dim)-> number of sub-treatments given size(...,3).
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): v4's plot_by_MHY reached into
% pmd.(txField)(:,:,4) and (:,:,(iSub-1)*5+(1:3)) with bare integers. Per
% lab convention (named indexing over integer indexing), those slices are
% now resolved through this lookup so the layer semantics documented in the
% dataloc schema (pTx layer 4 = Time, etc.) live in exactly one place.
% Adding a 6th per-sub-treatment layer upstream now means editing only this
% file, not every (:,:,N) in the plotters.
% -------------------------------------------------------------------------

L.block    = 5;
L.Name     = 1;
L.Conc     = 2;
L.Units    = 3;
L.Time     = 4;
L.TimeUnit = 5;
L.nameParts = [L.Name, L.Conc, L.Units];
L.txFields  = {'pTx','Tx1','Tx2','Tx3','Tx4'};

L.idx  = @(k, which) (k-1)*L.block + L.(which);
L.nsub = @(dim3) dim3 / L.block;
end
