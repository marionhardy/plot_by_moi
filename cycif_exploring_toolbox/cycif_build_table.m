function T = cycif_build_table(dataloc, comp, varargin)
% cycif_build_table  Flatten dataloc into a per-cell table with markers and
%                    parsed condition metadata. Raw intensities stored;
%                    any transformation is applied at plot time.
%
% INPUT
%   dataloc   — standard dataloc struct (post ct_cycif_link)
%   comp      — compartment string: 'nuc' | 'cyt' | 'cell'  (default 'nuc')
%
% NAME-VALUE PARAMETERS
%   'strict_platemap' — logical (default true). When true, only markers
%          explicitly assigned via the platemap staininfo grid (in
%          'ANTIBODY-CHANNEL' format) are populated for each well;
%          unassigned markers are NaN even if measured. This enforces
%          the platemap as design document.
%          When false, all measured markers in dataloc.d{xy}.data are
%          populated. Use this for fixed-mode workflows built with
%          ct_cycif_link_fixed's chan_map fallback, where the platemap
%          grid may be empty by design.
%
% OUTPUT
%   T  — table, one row per cell, columns:
%          xy, cellid, [markers...],
%          ptx_name, ptx_conc, ptx_units, ptx_time, ptx_timeunit, ptx_label,
%          ptx_pool_label, ptx_is_spike,
%          tx1_name ... txN_name (dynamic, one set per treatment slot),
%          tx1_label  (all treatment names combined, e.g. 'IFNg+LPS+HC')
%          condition  (ptx_pool_label + ' | ' + tx1_label, for cellxgene color-by)
%
% IDENTITY
%   Cells are uniquely identified by the (xy, cellid) pair, where cellid
%   matches dataloc.d{xy}.cellindex — enabling join-back to live traces.

if nargin < 2 || isempty(comp), comp = 'nuc'; end

% --- Parse name-value pairs ---
ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'strict_platemap', true, @islogical);
parse(ip, varargin{:});
strict_platemap = ip.Results.strict_platemap;

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
    % strict_platemap=true : only populate markers explicitly assigned via
    %                        platemap staininfo (design-document mode).
    %                        Markers measured in other rounds/wells but not
    %                        assigned here get NaN.
    % strict_platemap=false: populate any measured marker present in
    %                        d.data. Used for fixed-mode workflows where
    %                        chan_map (not the platemap grid) drives
    %                        marker assignment.
    if strict_platemap
        assigned = i_assigned_markers(dataloc.platemapd.staininfo, gr, gc);
    end

    mat = NaN(nCells, nMk);
    for m = 1:nMk
        if ~isfield(d.data, mk_fields{m}); continue; end
        if strict_platemap && ~ismember(markers{m}, assigned); continue; end
        mat(:,m) = d.data.(mk_fields{m});
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

% Dynamic Tx slots — reads all co-treatment blocks from pmd.Tx1
nTx1_slots = floor(size(pmd.Tx1, 3) / 5);
tx_names  = cell(1, nTx1_slots);
tx_concs  = cell(1, nTx1_slots);
tx_units  = cell(1, nTx1_slots);
tx_times  = cell(1, nTx1_slots);
tx_tunits = cell(1, nTx1_slots);
for k = 1:nTx1_slots
    base = (k-1)*5;
    tx_names{k}  = strtrim(safe(pmd.Tx1, r,c, base+1));
    tx_concs{k}  = strtrim(safe(pmd.Tx1, r,c, base+2));
    tx_units{k}  = strtrim(safe(pmd.Tx1, r,c, base+3));
    tx_times{k}  = strtrim(safe(pmd.Tx1, r,c, base+4));
    tx_tunits{k} = strtrim(safe(pmd.Tx1, r,c, base+5));
end

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