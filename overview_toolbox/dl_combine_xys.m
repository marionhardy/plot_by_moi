function [dataloc, idx2, linetp, xymat, goodxy, p] = dl_combine_xys(dataloc, txs, cellfn, idx, xymat, linetp, goodxy, p)
% DL_COMBINE_XYS  Pool replicate XYs sharing celltype x treatment into
%                 combined pseudo-XYs, for plotting.
%
% [dataloc, idx2, linetp, xymat, goodxy, p] = ...
%     dl_combine_xys(dataloc, txs, cellfn, idx, xymat, linetp, goodxy, p)
%
% INPUT  (positional, matches the v4 Combine_XYs interface exactly)
%   dataloc : dataloc struct; .d/.z/.IFd are pooled in place.
%   txs     : cellstr of treatment fieldnames (idx keys).
%   cellfn  : cellstr of celltype fieldnames (idx keys).
%   idx     : struct mapping condition name -> original xy list.
%   xymat   : cell of per-xy index vectors.
%   linetp  : cell of per-xy treatment-timepoint vectors (for xlines).
%   goodxy  : logical vector, which original xys hold data.
%   p       : param struct (uses p.ifchan).
%
% OUTPUT
%   dataloc : .d/.z/.IFd replaced with combined pseudo-XYs. Each combined
%             d{k} gains a .source = [origXY, origCellID] map so downstream
%             tools (e.g. ct_cell_video) can recover provenance.
%   idx2    : idx rebuilt against the NEW combined xy numbering.
%   linetp, xymat, goodxy : rebuilt to match combined numbering.
%   p       : passed through unchanged (returned for interface symmetry).
%
% -------------------------------------------------------------------------
% V5 CHANGE (DatalocHandler_v4 -> v5): this was a private subfunction
% (Combine_XYs) buried inside plot_by_MHY, silently duplicated into each v5
% plot_ function during the first refactor pass. Per the crit #6 decision
% it is now a SINGLE shared, documented helper in the plot layer (not three
% private copies, and not moved fully upstream). Behaviour is byte-for-byte
% identical to the v4 subfunction -- only the location, name, and this
% header changed.
%
% ARCHITECTURAL NOTE (flagged, not changed): combining replicates is an
% analysis-layer transform living in the plot layer. It rewrites the xy
% coordinate system so every index after this call refers to combined
% pseudo-XYs, not platemap xys -- hence the .source patch. A cleaner future
% design would be an upstream dl_combine_replicates producing a real
% dataloc; deferred by decision to keep it in the plot layer for now.
% -------------------------------------------------------------------------

fns = fieldnames(idx);
for ifn = 1:numel(fns); idx2.(fns{ifn}) = []; end
linetp2 = {};

totalXYs = 0;
for iCellline = 1:numel(cellfn)
    for iTxs = 1:numel(txs)
        XYs2Comb = intersect(idx.(txs{iTxs}),idx.(cellfn{iCellline}));
        if XYs2Comb > 0; totalXYs = totalXYs + 1; end
    end
end

PlotD = struct('d',[],'z',[],'IFd',[],'txx',[]);
PlotD.d = cell([1,totalXYs]); PlotD.z = cell([1,totalXYs]); PlotD.IFd = cell([1,totalXYs]);
source_per_xy = cell(1, totalXYs);
newXYNums = 1:totalXYs;
xyCounter = 0;
sizeD = size(dataloc.d,2);

for iCellline = 1:numel(cellfn)
    for iTxs = 1:size(txs,1)
        XYs2Comb = intersect(idx.(txs{iTxs}),idx.(cellfn{iCellline}));
        XYs2Comb = XYs2Comb((XYs2Comb <= sizeD));
        goodxys = arrayfun(@(x)~isempty(dataloc.d{x}),XYs2Comb);
        XYs2Comb = XYs2Comb(goodxys);
        numXYs = length(XYs2Comb);
        if ~isempty(XYs2Comb)
            fieldNames = fieldnames(dataloc.d{XYs2Comb(1)}.data);
            if any(contains(fieldNames, 'Corr')); fieldNames = fieldNames(~contains(fieldNames, 'Corr')); end
            NumChans = length(fieldNames);
            firstxy = true;
            for iXY = 1:numXYs
                ThisXY = XYs2Comb(iXY);
                if goodxy(ThisXY) && ThisXY <= length(dataloc.d)
                    if firstxy
                        xyCounter = xyCounter + 1; newXYNum = newXYNums(xyCounter);
                        idx2.(txs{iTxs}) = [idx2.(txs{iTxs}), newXYNum];
                    end
                    if ~isempty(dataloc.d{ThisXY}) && isfield(dataloc.d{ThisXY}, 'data')
                        if ~isfield(PlotD.d{newXYNum},'data'); PlotD.d{newXYNum}.data = struct(); end
                        if ~isfield(PlotD.z{newXYNum},'data'); PlotD.z{newXYNum}.data = struct(); end
                        if ~isfield(PlotD.IFd{newXYNum},'data'); PlotD.IFd{newXYNum}.data = struct(); end
                        for iChan = 1:NumChans
                            ThisChan = fieldNames{iChan};
                            if isfield(dataloc.d{ThisXY}.data, ThisChan)
                                HasZ = false;
                                if isfield(dataloc, 'z') && ~isempty(dataloc.z) && ThisXY <= length(dataloc.d)
                                    if isfield(dataloc.z{ThisXY}, 'data') && ~isempty(dataloc.z{ThisXY}.data)
                                        HasZ = isfield(dataloc.z{ThisXY}.data, ThisChan);
                                    end
                                end
                                HasIFd = false;
                                if isfield(dataloc, 'IFd') && ~isempty(dataloc.IFd)
                                    HasIFd = isfield(dataloc.IFd{ThisXY}.data, p.ifchan);
                                end
                                if ~isfield(PlotD.d{newXYNum}.data, ThisChan)
                                    PlotD.d{newXYNum}.data.(ThisChan) = dataloc.d{ThisXY}.data.(ThisChan);
                                    if iChan == 1
                                        nRows = size(dataloc.d{ThisXY}.data.(ThisChan), 1);
                                        source_per_xy{newXYNum} = [source_per_xy{newXYNum}; repmat(ThisXY, nRows, 1), (1:nRows)'];
                                    end
                                    if firstxy
                                        linetp2{newXYNum} = linetp{cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0))};
                                        firstxy = false;
                                    end
                                    if HasZ; PlotD.z{newXYNum}.data.(ThisChan) = dataloc.z{ThisXY}.data.(ThisChan); end
                                    if HasIFd && iChan < 2; PlotD.IFd{newXYNum}.data.(p.ifchan) = dataloc.IFd{ThisXY}.data.(p.ifchan); end
                                else
                                    PlotD.d{newXYNum}.data.(ThisChan) = [PlotD.d{newXYNum}.data.(ThisChan); dataloc.d{ThisXY}.data.(ThisChan)];
                                    if HasZ; PlotD.z{newXYNum}.data.(ThisChan) = [PlotD.z{newXYNum}.data.(ThisChan); dataloc.z{ThisXY}.data.(ThisChan)]; end
                                    if HasIFd && iChan < 2; PlotD.IFd{newXYNum}.data.(p.ifchan) = [PlotD.IFd{newXYNum}.data.(p.ifchan); dataloc.IFd{ThisXY}.data.(p.ifchan)]; end
                                end
                            end
                        end
                    end
                end
            end
        end
        if numel(XYs2Comb) > 0 && exist('newXYNum', 'var'); idx2.(cellfn{iCellline}) = [idx2.(cellfn{iCellline}), newXYNum]; end
    end
end

dataloc.d = PlotD.d;
for iSrc = 1:numel(source_per_xy)
    if ~isempty(dataloc.d{iSrc}) && isfield(dataloc.d{iSrc},'data')
        dataloc.d{iSrc}.source = source_per_xy{iSrc};
    end
end
dataloc.IFd = PlotD.IFd;
dataloc.z = PlotD.z;
goodxy = ~cellfun(@isempty, dataloc.d)';
linetp = linetp2';
xymat = cell(numel(newXYNums),1);
for ixymat = 1:numel(newXYNums); xymat{ixymat} = newXYNums(ixymat); end
end
