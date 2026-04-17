function T = cycif_build_table(dataloc, comp)
% cycif_build_table  Flatten dataloc into a per-cell table with markers and
%                    parsed condition metadata. Raw intensities stored;
%                    any transformation is applied at plot time.
%
% INPUT
%   dataloc   — standard dataloc struct (post ct_cycif_link)
%   comp      — compartment string: 'nuc' | 'cyt' | 'cell'  (default 'nuc')
%
% OUTPUT
%   T  — table, one row per cell, columns:
%          xy, cellid, [markers...],
%          ptx_name, ptx_conc, ptx_units, ptx_time, ptx_timeunit, ptx_label,
%          tx1_name, tx1_conc, tx1_units, tx1_time, tx1_timeunit,
%          tx2_name, tx2_conc, tx2_units, tx2_time, tx2_timeunit,
%          tx1_label  (tx1_name + tx2_name combined, e.g. 'IFNg+LPS')
%
% IDENTITY
%   Cells are uniquely identified by the (xy, cellid) pair, where cellid
%   matches dataloc.d{xy}.cellindex — enabling join-back to live traces.

if nargin < 2 || isempty(comp), comp = 'nuc'; end

pmd = dataloc.platemapd.pmd;

% --- Build flat XY -> grid-position lookup --------------------------------
xy_grid = nan(8,12);
for r = 1:8
    for c = 1:12
        v = pmd.xy{r,c};
        if ~isempty(v) && ~isnan(v), xy_grid(r,c) = v; end
    end
end

% --- Discover markers: union across ALL valid XYs -------------------------
% Required because different wells receive different antibody subsets.
% Wells missing a given marker will get NaN (handled in intensity loop).
suffix  = ['_' comp];
sfx_len = numel(suffix);

assert(any(cellfun(@(d) ~isempty(d) && isfield(d,'data'), dataloc.d)), ...
       'No valid dataloc.d entries found.');

mk_fields = {};
for xy_s = 1:numel(dataloc.d)
    d_s = dataloc.d{xy_s};
    if isempty(d_s) || ~isfield(d_s,'data'), continue; end
    fns     = fieldnames(d_s.data);
    new_fns = fns(cellfun(@(f) numel(f) > sfx_len && ...
                           strcmp(f(end-sfx_len+1:end), suffix), fns));
    mk_fields = union(mk_fields, new_fns, 'stable');
end
markers = cellfun(@(f) f(1:end-sfx_len), mk_fields, 'UniformOutput', false);
nMk     = numel(markers);
fprintf('[cycif_build_table] %d markers (%s): %s\n', nMk, comp, strjoin(markers,', '));

% --- Accumulate rows ------------------------------------------------------
nXY       = numel(dataloc.d);
row_data  = {};
row_cids  = {};
row_meta  = {};

for xy = 1:nXY
    d = dataloc.d{xy};
    if isempty(d) || ~isfield(d,'data'), continue; end
    % No early exit on a specific marker — wells have heterogeneous panels.
    % Missing markers get NaN in the intensity loop below.

    % Grid position — needed for both staininfo masking and condition metadata
    [gr, gc] = find(xy_grid == xy, 1);
    if isempty(gr)
        warning('XY %d not found in pmd.xy — skipping.', xy); continue;
    end

    % Determine nCells from the first marker present in this XY
    present = mk_fields(cellfun(@(f) isfield(d.data,f), mk_fields));
    if isempty(present), continue; end   % no cycIF data at all for this XY
    nCells = numel(d.data.(present{1}));
    if nCells < 1, continue; end

    % Intensity matrix [nCells x nMk]
    % Only populate markers that staininfo explicitly assigns to this well.
    % Markers measured in the same stain round but assigned to other wells
    % (e.g. VEGF in iNOS rows) reflect channel background, not protein
    % expression — they are set to NaN here to encode experimental design.
    assigned = i_assigned_markers(dataloc.platemapd.staininfo, gr, gc);

    mat = NaN(nCells, nMk);
    for m = 1:nMk
        if isfield(d.data, mk_fields{m}) && ismember(markers{m}, assigned)
            mat(:,m) = d.data.(mk_fields{m});
        end
    end

    % Cell identity — use cellindex if present and conformant, else 1:N
    if isfield(d, 'cellindex') && numel(d.cellindex) == nCells
        cid = d.cellindex(:);
    else
        cid = (1:nCells)';
    end

    cond     = i_parse_cond(pmd, gr, gc);
    meta_row = repmat(cond, nCells, 1);
    meta_row.xy = repmat(xy, nCells, 1);

    row_data{end+1} = mat;       %#ok<AGROW>
    row_cids{end+1} = cid;       %#ok<AGROW>
    row_meta{end+1} = meta_row;  %#ok<AGROW>
end

assert(~isempty(row_data), 'No data collected — check comp/dataloc.');

% --- Assemble final table -------------------------------------------------
Tmarkers = array2table(vertcat(row_data{:}), 'VariableNames', markers);
Tmeta    = vertcat(row_meta{:});
cellid   = vertcat(row_cids{:});

% Column order: xy | cellid | markers | condition metadata
T = [table(Tmeta.xy, cellid, 'VariableNames', {'xy','cellid'}), ...
     Tmarkers, ...
     Tmeta(:, setdiff(Tmeta.Properties.VariableNames, {'xy'}, 'stable'))];

fprintf('[cycif_build_table] %d cells | %d unique tx1_labels\n', ...
        height(T), numel(unique(T.tx1_label)));
end


% =========================================================================
%  PRIVATE HELPERS
% =========================================================================

function cond = i_parse_cond(pmd, r, c)
% Extract and format condition metadata at plate position (r,c).

    function s = safe(arr, r, c, k)
        if size(arr,3) >= k
            v = arr{r,c,k};
            if isempty(v) || (isnumeric(v) && isnan(v)), s = ''; else, s = char(v); end
        else
            s = '';
        end
    end

% pTx
ptx_name  = safe(pmd.pTx, r,c,1);
ptx_conc  = safe(pmd.pTx, r,c,2);
ptx_units = safe(pmd.pTx, r,c,3);
ptx_time  = safe(pmd.pTx, r,c,4);
ptx_tunit = safe(pmd.pTx, r,c,5);

% ptx_pool_label: time-stripped for pooling e.g. "Glucose 17mM"
% ptx_label:      full label with timing e.g. "Glucose 17mM (-2h)"
% Wells where timeunit == 'tp' are spike treatments — flagged via ptx_is_spike.
if ~isempty(ptx_name)
    ptx_pool_label = ptx_name;
    if ~isempty(ptx_conc), ptx_pool_label = [ptx_pool_label ' ' ptx_conc ptx_units]; end

    ptx_label = ptx_pool_label;
    if ~isempty(ptx_time), ptx_label = [ptx_label ' (' ptx_time ptx_tunit ')']; end
else
    ptx_label      = 'none';
    ptx_pool_label = 'none';
end

% Flag spike treatments (timeunit 'tp' = added at a timepoint mid-experiment,
% not a standard pre-treatment). These pool under the same ptx_pool_label
% but are worth excluding or examining separately.
ptx_is_spike = strcmp(strtrim(ptx_tunit), 'tp');

% Tx1 (treatment 1, slots 1-5)
tx1_name  = safe(pmd.Tx1, r,c,1);
tx1_conc  = safe(pmd.Tx1, r,c,2);
tx1_units = safe(pmd.Tx1, r,c,3);
tx1_time  = safe(pmd.Tx1, r,c,4);
tx1_tunit = safe(pmd.Tx1, r,c,5);

% Tx2 (treatment 2, slots 6-10)
tx2_name  = safe(pmd.Tx1, r,c,6);
tx2_conc  = safe(pmd.Tx1, r,c,7);
tx2_units = safe(pmd.Tx1, r,c,8);
tx2_time  = safe(pmd.Tx1, r,c,9);
tx2_tunit = safe(pmd.Tx1, r,c,10);

% Combined Tx label e.g. 'IFNg+LPS' | 'IL-4' | 'none'
parts = {};
if ~isempty(tx1_name), parts{end+1} = tx1_name; end
if ~isempty(tx2_name), parts{end+1} = tx2_name; end
tx1_label = strjoin(parts, '+');
if isempty(tx1_label), tx1_label = 'none'; end

cond = table( ...
    string(ptx_name),       string(ptx_conc),  string(ptx_units), ...
    string(ptx_time),       string(ptx_tunit), string(ptx_label), ...
    string(ptx_pool_label), ptx_is_spike, ...
    string(tx1_name),  string(tx1_conc),  string(tx1_units), ...
    string(tx1_time),  string(tx1_tunit), ...
    string(tx2_name),  string(tx2_conc),  string(tx2_units), ...
    string(tx2_time),  string(tx2_tunit), ...
    string(tx1_label), ...
    'VariableNames', { ...
        'ptx_name','ptx_conc','ptx_units','ptx_time','ptx_timeunit','ptx_label', ...
        'ptx_pool_label','ptx_is_spike', ...
        'tx1_name','tx1_conc','tx1_units','tx1_time','tx1_timeunit', ...
        'tx2_name','tx2_conc','tx2_units','tx2_time','tx2_timeunit', ...
        'tx1_label'});
end

function assigned = i_assigned_markers(staininfo, plate_row, plate_col)
% Return the set of marker names legitimately assigned to this plate position
% across all stain rounds, parsed from staininfo.markers{plate_row, plate_col}.
% Format: 'ANTIBODY-CHANNEL' semicolon-delimited, e.g. 'iNOS-CHERRY; CD206-ORANGE'
assigned = {};
for r = 1:numel(staininfo)
    if plate_row > size(staininfo(r).markers, 1) || ...
       plate_col > size(staininfo(r).markers, 2)
        continue;
    end
    mstr = strtrim(staininfo(r).markers{plate_row, plate_col});
    if isempty(mstr), continue; end
    entries = strtrim(strsplit(mstr, ';'));
    for ie = 1:numel(entries)
        entry = strtrim(entries{ie});
        if isempty(entry), continue; end
        % Parse ANTIBODY-CHANNEL — last hyphen-delimited token is channel
        parts = strsplit(entry, '-');
        if numel(parts) < 2, continue; end
        ab = strtrim(strjoin(parts(1:end-1), '-'));
        if ~isempty(ab)
            assigned{end+1} = ab; %#ok<AGROW>
        end
    end
end
end