function R = osc_audit(dataloc, GroupStruct, varargin)
% OSC_AUDIT — single-file sanity/quality check for oscillator labels.
% Usage:
%   R = osc_audit(f2Data{1}, CellGroups.<exp>, ...
%       'channel','HYLIGHT','pulsethresh',0.5,'hours',[0 12], ...
%       'aftertreatment',1,'xy',[],'showplots',false);

%% Args
ip = inputParser; ip.CaseSensitive=false;
addParameter(ip,'channel','HYLIGHT');
addParameter(ip,'pulsethresh',0.5,@isscalar);           % pulses/hour threshold
addParameter(ip,'hours',[0 12],@(v) isnumeric(v)&&numel(v)==2); % window rel. to anchor tx
addParameter(ip,'aftertreatment',1,@isscalar);          % which tx to anchor to
addParameter(ip,'xy',[],@(v) isempty(v)||isvector(v));  % restrict XYs
addParameter(ip,'showplots',false,@islogical);
parse(ip,varargin{:}); p=ip.Results; chan=char(p.channel);

%% Sets & timing
osc = unique(GroupStruct.Oscillatory(:));
non = unique(GroupStruct.NonOscillatory(:));
fph = local_frames_per_hour(dataloc);   % frames/hour
allIDs = local_all_ids(dataloc);

%% A) Basic checks
overlap = intersect(osc,non);
missing_osc = setdiff(osc,allIDs);
missing_non = setdiff(non,allIDs);

fprintf('--- Basic checks ---\n');
fprintf('Osc IDs: %d | Non IDs: %d | Overlap: %d\n', numel(osc), numel(non), numel(overlap));
if ~isempty(overlap)
    fprintf('  OVERLAP EXAMPLE (first 10): %s\n', mat2str(overlap(1:min(10,end))));
end
fprintf('Missing in data — Osc: %d | Non: %d\n', numel(missing_osc), numel(missing_non));

%% B) Per-XY label coverage (union-aware unlabeled)
goodXY = find(cellfun(@(D) ~isempty(D) && isfield(D,'cellindex') && ~isempty(D.cellindex), dataloc.d));
if ~isempty(p.xy), goodXY = intersect(goodXY(:), p.xy(:)); end
mapTbl = [];
U = union(osc,non);
for xy = goodXY(:)'
    ids = dataloc.d{xy}.cellindex(:);
    nOsc = nnz(ismember(ids, osc));
    nNon = nnz(ismember(ids, non));
    nLab = nnz(ismember(ids, U));
    nUnl = numel(ids) - nLab;
    if (nOsc+nNon)>0
        mapTbl = [mapTbl; xy, nOsc, nNon, nUnl]; %#ok<AGROW>
    end
end
fprintf('XY coverage rows: %d (has ≥1 labeled cell)\n', size(mapTbl,1));
if ~isempty(mapTbl)
    headN = min(10,size(mapTbl,1));
    fprintf('  XY  Osc  Non  Unlabeled (first %d)\n', headN);
    fprintf('  %3d %4d %4d %9d\n', mapTbl(1:headN,:).');
end

%% C) Pulse-rate audit (recompute from z)
tot_mism_osc = 0; tot_mism_non = 0;
freq_osc = []; freq_non = [];
for xy = goodXY(:)'
    D = dataloc.d{xy};
    if ~isfield(D,'data') || ~isfield(D.data,chan), continue; end
    T = size(D.data.(chan),2);
    ids = D.cellindex(:);

    % anchor tx frame
    txFrames = local_tx_frames_for_xy(dataloc, xy); % vector (may be empty)
    if isempty(txFrames), anchor = 1;
    else, anchor = max(1, round(txFrames(min(p.aftertreatment, numel(txFrames))))); end

    % window in frames
    t0 = max(1, anchor + floor(p.hours(1)*fph));
    t1 = min(T, anchor + floor(p.hours(2)*fph));
    if t1 <= t0, t0=1; t1=T; end
    hours = max(eps,(t1-t0+1)/fph);

    % z present?
    hasZ = isfield(dataloc,'z') && numel(dataloc.z)>=xy && ~isempty(dataloc.z{xy}) ...
        && isfield(dataloc.z{xy},'data') && isfield(dataloc.z{xy}.data,chan);
    if ~hasZ, continue; end
    Z = dataloc.z{xy}.data.(chan);

    % pulses/hour per cell
    n = numel(ids); freq = nan(n,1);
    for i = 1:n
        pk = [];
        if isfield(Z(i),'pkpos') && ~isempty(Z(i).pkpos), pk = Z(i).pkpos(:); end
        inwin = (pk>=t0) & (pk<=t1);
        freq(i) = numel(pk(inwin))/hours;
    end

    isOsc = ismember(ids,osc);
    isNon = ismember(ids,non);

    freq_osc = [freq_osc; freq(isOsc)]; %#ok<AGROW>
    freq_non = [freq_non; freq(isNon)]; %#ok<AGROW>

    tot_mism_osc = tot_mism_osc + nnz(isOsc & (freq <  p.pulsethresh));
    tot_mism_non = tot_mism_non + nnz(isNon & (freq >= p.pulsethresh));
end

fprintf('Pulse audit (%.2f pulses/hour, window [%g %g] h, tx#%d anchor):\n', ...
    p.pulsethresh, p.hours(1), p.hours(2), p.aftertreatment);
fprintf('  Labeled Osc but below thresh: %d\n', tot_mism_osc);
fprintf('  Labeled Non but above/equal:  %d\n', tot_mism_non);

mOsc = mean(freq_osc,'omitnan'); sOsc = std(freq_osc,[],'omitnan');
mNon = mean(freq_non,'omitnan'); sNon = std(freq_non,[],'omitnan');
fprintf('  Mean freq — Osc: %.3f ± %.3f (n=%d) | Non: %.3f ± %.3f (n=%d)\n', ...
    mOsc,sOsc,nnz(~isnan(freq_osc)), mNon,sNon,nnz(~isnan(freq_non)));

%% D) Optional plots
if p.showplots && (~isempty(freq_osc) || ~isempty(freq_non))
    figure('Color','w'); hold on;
    if ~isempty(freq_osc), h1 = histogram(freq_osc,'Normalization','pdf','FaceAlpha',0.5); h1.DisplayName='Osc'; end
    if ~isempty(freq_non), h2 = histogram(freq_non,'Normalization','pdf','FaceAlpha',0.5); h2.DisplayName='Non'; end
    xline(p.pulsethresh,'k--','DisplayName','Thresh');
    xlabel('pulses/hour'); ylabel('pdf'); legend('Location','best');
    title(sprintf('%s pulse freq', chan));
end

%% Return
R = struct();
R.overlap_ids = overlap;
R.missing_osc = missing_osc;
R.missing_non = missing_non;
R.xy_counts = mapTbl;              % [xy, nOsc, nNon, nUnlabeled]
R.freq_osc = freq_osc;
R.freq_non = freq_non;
R.mismatch_osc_below = tot_mism_osc;
R.mismatch_non_above = tot_mism_non;
R.params = p;
end

% ============================ Helpers ============================

function fph = local_frames_per_hour(dataloc)
f = 6;
if isfield(dataloc,'movieinfo') && isfield(dataloc.movieinfo,'tsamp') && ~isempty(dataloc.movieinfo.tsamp)
    f = dataloc.movieinfo.tsamp;
end
fph = 60/max(eps,f);
end

function ids = local_all_ids(dataloc)
ids = [];
for i=1:numel(dataloc.d)
    D = dataloc.d{i};
    if ~isempty(D) && isfield(D,'cellindex') && ~isempty(D.cellindex)
        ids = [ids; D.cellindex(:)]; %#ok<AGROW>
    end
end
ids = unique(ids);
end

function txFrames = local_tx_frames_for_xy(dataloc, xy)
% Collect treatment-frame times for one XY from platemap (plane 4).
txFrames = [];
if ~isfield(dataloc,'platemapd') || ~isfield(dataloc.platemapd,'pmd'), return; end
pmd = dataloc.platemapd.pmd;
if ~isfield(pmd,'xy') || isempty(pmd.xy) || ~isfield(pmd,'Cell') || isempty(pmd.Cell), return; end

[H,W,~] = size(pmd.Cell);
r=NaN; c=NaN;
for rr=1:H
    for cc=1:W
        xys = local_safe_xylist(pmd, rr, cc);
        if any(isfinite(xys) & (round(double(xys))==xy)), r=rr; c=cc; break; end
    end
    if ~isnan(r), break; end
end
if isnan(r), return; end

txFields = {'pTx','Tx1','Tx2','Tx3','Tx4'};
for s = 1:numel(txFields)
    f = txFields{s};
    if ~isfield(pmd,f) || isempty(pmd.(f)), continue; end
    A = pmd.(f);
    if size(A,3) >= 4
        val = A(r,c,4);
        if iscell(val),   val = val{1}; end
        if isstring(val), val = char(val); end
        if ischar(val)
            num = str2double(val);
            if ~isnan(num), txFrames(end+1) = double(num); end %#ok<AGROW>
        elseif isnumeric(val)
            v = val(:); v = v(isfinite(v));
            if ~isempty(v), txFrames(end+1) = double(v(1)); end %#ok<AGROW>
        end
    end
end
txFrames = unique(txFrames(:)','stable');
end

function xys = local_safe_xylist(pmd, r, c)
xys = [];
try
    xys = pmd.xy{r,c};
catch
    k = sub2ind(size(pmd.Cell), r, c);
    if k<=numel(pmd.xy), xys = pmd.xy{k}; end
end
end
