function T = cycif_build_table(dataloc, comp)
% cycif_build_table  Flatten dataloc into a per-cell table with markers and
%                    parsed condition metadata. Raw intensities stored;
%                    any transformation is applied at plot time.
%
% INPUT
%   dataloc   — standard dataloc struct (post ct_cycif_link)
%   comp      — compartment specification, one of:
%                 'nuc' | 'cyt' | 'cell'        (single compartment, char)
%                 {'nuc','cyt'}                  (multiple, cell array)
%                 {'nuc','cyt','ratio_cn'}       (adds cyt/nuc ratio)
%                 {'nuc','cyt','ratio_cn','ratio_nc'}  (both directions)
%                 default: 'nuc'
%
%              Ratios require BOTH _nuc and _cyt fields on the marker;
%              markers with only one compartment measured get NaN in the
%              ratio column. A small epsilon (1e-9) is added to the
%              denominator to avoid Inf/NaN from divide-by-zero.
%
% OUTPUT
%   T  — table, one row per cell, columns:
%          xy, cellid, [markers, one per compartment...],
%          ptx_name, ptx_conc, ptx_units, ptx_time, ptx_timeunit, ptx_label,
%          ptx_pool_label, ptx_is_spike,
%          tx1_name ... txN_name (dynamic, one set per treatment slot),
%          tx1_label  (all treatment names combined, e.g. 'IFNg+LPS+HC')
%          condition  (ptx_pool_label + ' | ' + tx1_label, for cellxgene color-by)
%
%   Marker column names: '<marker>_<compartment>' (e.g. 'GFP_nuc', 'GFP_cyt',
%   'GFP_ratio_cn'). Compartment order in the table follows the order in
%   comp.
%
% IDENTITY
%   Cells are uniquely identified by the (xy, cellid) pair, where cellid
%   matches dataloc.d{xy}.cellindex — enabling join-back to live traces.

if nargin < 2 || isempty(comp), comp = 'nuc'; end

% Normalize comp to a cell array of compartment strings
if ischar(comp) || (isstring(comp) && isscalar(comp))
    comps = {char(comp)};
elseif iscell(comp)
    comps = cellfun(@char, comp, 'UniformOutput', false);
else
    error('cycif_build_table:badComp', ...
        'comp must be a char, string, or cell array of strings.');
end

% Validate each compartment token
valid_comps = {'nuc', 'cyt', 'cell', 'ratio_cn', 'ratio_nc'};
for k = 1:numel(comps)
    if ~any(strcmp(comps{k}, valid_comps))
        error('cycif_build_table:badComp', ...
            'Unknown compartment "%s". Valid: %s.', comps{k}, strjoin(valid_comps, ', '));
    end
end
n_comps = numel(comps);

% Ratios (if any) need eps to avoid divide-by-zero
RATIO_EPS = 1e-9;

pmd = dataloc.platemapd.pmd;

% --- Build flat XY -> grid-position lookup --------------------------------
xy_grid = nan(8,12);
for r = 1:8
    for c = 1:12
        v = pmd.xy{r,c};
        if ~isempty(v) && ~isnan(v), xy_grid(r,c) = v; end
    end
end

% --- Discover markers: union across ALL valid XYs and ALL compartments ---
% Required because different wells receive different antibody subsets and
% we may need multiple compartments per marker.
%
% Strip _nuc/_cyt/_cell suffix from every field to get bare marker names.
% Ratios don't have their own storage -- they're derived from _nuc and _cyt.

assert(any(cellfun(@(d) ~isempty(d) && isfield(d,'data'), dataloc.d)), ...
       'No valid dataloc.d entries found.');

storage_suffixes = {'_nuc', '_cyt', '_cell'};
markers = {};
for xy_s = 1:numel(dataloc.d)
    d_s = dataloc.d{xy_s};
    if isempty(d_s) || ~isfield(d_s,'data'), continue; end
    fns = fieldnames(d_s.data);
    for kk = 1:numel(fns)
        f = fns{kk};
        for ss = 1:numel(storage_suffixes)
            sfx = storage_suffixes{ss};
            if numel(f) > numel(sfx) && strcmp(f(end-numel(sfx)+1:end), sfx)
                markers = union(markers, {f(1:end-numel(sfx))}, 'stable');
                break;
            end
        end
    end
end
nMk = numel(markers);

% Build the ordered list of output columns: for each marker, one column per
% requested compartment. Column names are '<marker>_<compartment>'.
n_cols_markers = nMk * n_comps;
col_names = cell(1, n_cols_markers);
for m = 1:nMk
    for cc = 1:n_comps
        col_names{(m-1)*n_comps + cc} = [markers{m} '_' comps{cc}];
    end
end

fprintf('[cycif_build_table] %d markers (%s): %s\n', ...
        nMk, strjoin(comps, ','), strjoin(markers, ', '));

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

    % Determine nCells from the first storage field present in this XY.
    % Any '<marker>_<nuc|cyt|cell>' field can serve as the size probe.
    nCells = 0;
    for m = 1:nMk
        for ss = 1:numel(storage_suffixes)
            fld = [markers{m} storage_suffixes{ss}];
            if isfield(d.data, fld)
                nCells = numel(d.data.(fld));
                break;
            end
        end
        if nCells > 0; break; end
    end
    if nCells < 1, continue; end   % no cycIF data at all for this XY

    % Intensity matrix [nCells x (nMk * n_comps)]
    % Only populate markers that staininfo explicitly assigns to this well.
    % Markers measured in the same stain round but assigned to other wells
    % (e.g. VEGF in iNOS rows) reflect channel background, not protein
    % expression — they are set to NaN here to encode experimental design.
    assigned = i_assigned_markers(dataloc.platemapd.staininfo, gr, gc);

    mat = NaN(nCells, nMk * n_comps);
    for m = 1:nMk
        if ~ismember(markers{m}, assigned); continue; end
        % Pull nuc and cyt once each if either is needed (for ratios or direct)
        nuc_field = [markers{m} '_nuc'];
        cyt_field = [markers{m} '_cyt'];
        has_nuc = isfield(d.data, nuc_field);
        has_cyt = isfield(d.data, cyt_field);
        for cc = 1:n_comps
            col_idx = (m-1)*n_comps + cc;
            switch comps{cc}
                case {'nuc', 'cyt', 'cell'}
                    fld = [markers{m} '_' comps{cc}];
                    if isfield(d.data, fld)
                        mat(:, col_idx) = d.data.(fld);
                    end
                case 'ratio_cn'
                    if has_nuc && has_cyt
                        mat(:, col_idx) = d.data.(cyt_field) ./ ...
                                          (d.data.(nuc_field) + RATIO_EPS);
                    end
                case 'ratio_nc'
                    if has_nuc && has_cyt
                        mat(:, col_idx) = d.data.(nuc_field) ./ ...
                                          (d.data.(cyt_field) + RATIO_EPS);
                    end
            end
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
Tmarkers = array2table(vertcat(row_data{:}), 'VariableNames', col_names);
Tmeta    = vertcat(row_meta{:});
cellid   = vertcat(row_cids{:});

% Column order: xy | cellid | markers | condition metadata
T = [table(Tmeta.xy, cellid, 'VariableNames', {'xy','cellid'}), ...
     Tmarkers, ...
     Tmeta(:, setdiff(Tmeta.Properties.VariableNames, {'xy'}, 'stable'))];

% Add combined condition label for cellxgene color-by.
% e.g. 'Glucose 17mM | IFNg+LPS+HC'
T.condition = T.ptx_pool_label + ' | ' + T.tx1_label;

fprintf('[cycif_build_table] %d cells | %d unique tx1_labels | %d unique conditions\n', ...
        height(T), numel(unique(T.tx1_label)), numel(unique(T.condition)));
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

% --- Dynamic Tx slots -- read every pmd.Tx<N> section ---
% Historically only pmd.Tx1 was read. But iman_readdatasheet splits multiple
% "Tx" blocks in the xlsx into pmd.Tx1, pmd.Tx2, ..., pmd.TxN. Iterate over
% all of them and concatenate their per-well slots. Column names follow a
% global counter (tx1_*, tx2_*, ...); when each Tx<N> section holds one
% slot, tx<N>_ maps to Tx<N> naturally.
tx_field_names = fieldnames(pmd);
tx_field_names = tx_field_names( ...
    ~cellfun('isempty', regexp(tx_field_names, '^Tx\d+$', 'once')));
% Sort by numeric suffix so Tx2 comes after Tx1 (not lexicographic)
[~, tx_order] = sort(cellfun(@(s) sscanf(s, 'Tx%d'), tx_field_names));
tx_field_names = tx_field_names(tx_order);

tx_names  = {};
tx_concs  = {};
tx_units  = {};
tx_times  = {};
tx_tunits = {};
for it = 1:numel(tx_field_names)
    Tx = pmd.(tx_field_names{it});
    n_slots_here = floor(size(Tx, 3) / 5);
    for k = 1:n_slots_here
        base = (k-1)*5;
        tx_names{end+1}  = strtrim(safe(Tx, r,c, base+1)); %#ok<AGROW>
        tx_concs{end+1}  = strtrim(safe(Tx, r,c, base+2)); %#ok<AGROW>
        tx_units{end+1}  = strtrim(safe(Tx, r,c, base+3)); %#ok<AGROW>
        tx_times{end+1}  = strtrim(safe(Tx, r,c, base+4)); %#ok<AGROW>
        tx_tunits{end+1} = strtrim(safe(Tx, r,c, base+5)); %#ok<AGROW>
    end
end
nTx1_slots = numel(tx_names);   % keep the same name for downstream loop

% Combined label: all non-empty treatment names joined with '+'
non_empty = tx_names(~cellfun(@isempty, tx_names));
non_empty = cellfun(@strtrim, non_empty, 'UniformOutput', false);
tx1_label = strjoin(non_empty, '+');
if isempty(tx1_label), tx1_label = 'none'; end

% Build output table: fixed columns + dynamic tx slot columns
var_names = {'ptx_name','ptx_conc','ptx_units','ptx_time','ptx_timeunit', ...
             'ptx_label','ptx_pool_label','ptx_is_spike'};
var_vals  = {string(ptx_name), string(ptx_conc), string(ptx_units), ...
             string(ptx_time), string(ptx_tunit), ...
             string(ptx_label), string(ptx_pool_label), ptx_is_spike};

for k = 1:nTx1_slots
    prefix = sprintf('tx%d', k);
    var_names = [var_names, {[prefix '_name'],[prefix '_conc'], ...
                              [prefix '_units'],[prefix '_time'],[prefix '_timeunit']}]; %#ok<AGROW>
    var_vals  = [var_vals,  {string(tx_names{k}), string(tx_concs{k}), ...
                              string(tx_units{k}), string(tx_times{k}), string(tx_tunits{k})}]; %#ok<AGROW>
end
var_names{end+1} = 'tx1_label';
var_vals{end+1}  = string(tx1_label);

cond = table(var_vals{:}, 'VariableNames', var_names);
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