function report = dl_check_platemap_compat(dataloc)
% DL_CHECK_PLATEMAP_COMPAT  Report-only legacy/v5 compatibility check for dataloc.
%
% report = dl_check_platemap_compat(dataloc)
%
% Run this yourself in the console before plotting, e.g.:
%   report = dl_check_platemap_compat(dataloc);
%
% Deliberately NOT called automatically from plot_mean / plot_heatmap /
% plot_sorted_stacks: a plotting function silently mutating your dataloc,
% or silently re-reading files, is a side effect you didn't ask for.
% This is a pre-flight step you run and read, on your own terms.
%
% WHAT IT CHECKS (see "2026-07-06 dataloc legacy vs v5 compatibility.md"
% in the project notebook for full background):
%
%   1. dataloc.version -- schema completeness axis. Legacy = missing or
%      '0'-'4' (dl_config_paths defaults missing version to '4';
%      dl_load_dataloc upgrades '0'-'4' -> '5' via dl_upgrade_dataloc).
%      Mostly harmless to plotting (adds .exp/.seg/.filterp scaffolding
%      the plot_ functions never read), reported for completeness only.
%
%   2. Literal 'NA' contamination in dataloc.platemapd.pmd -- the axis
%      that actually matters. dl_ensure_platemap.m's local_sanitize_na
%      scrubs the v5 unimaged-well 'NA' convention out of pmd, but ONLY
%      when it re-reads the platemap xlsx (empty pmd, stale mtime, or
%      opts.updatemap=true) -- this is NOT part of dl_upgrade_dataloc.
%      A dataloc.version=='5' can still carry unsanitized 'NA' if it was
%      last loaded before the fix landed, or the cached .pmd was never
%      invalidated. 'NA' is ischar, so it survives the Cell/Gene
%      concatenation filter in plot_by_MHY/plot_mean/plot_sorted_stacks
%      exactly like a real annotation -- this is the confirmed root cause
%      of the repeating-token "wall of text" titles.
%      Detection predicate mirrors local_sanitize_na EXACTLY:
%      ischar(v) && strcmp(strtrim(v),'NA'), applied per pmd field.
%
%   3. Missing/mismatched dataloc.d{xy}.cellindex -- "(optional)" per
%      ct_datafilt_mhy.m's own docstring. plot_ functions now fabricate +
%      warn per-xy at plot time (matching dl_window_stats.m's
%      cws_get_cellindex convention); this just gives you the count
%      up front so warnings at plot time aren't a surprise.
%
% INPUT
%   dataloc : the dataloc struct, any version.
%
% OUTPUT
%   report  : struct with fields
%     .version         (char)   dataloc.version, or 'missing (legacy)'
%     .isLegacySchema  (logical) true if version missing or < '5'
%     .naCounts        (struct)  per-pmd-field count of contaminated cells
%     .naTotal         (double)  total contaminated cells across all fields
%     .missingCellIndexXYs (double vector) xy indices lacking cellindex
%   Printed to console regardless; returned for programmatic use (e.g.
%   gating a script on report.naTotal == 0).
%
% NOTHING IN THIS FUNCTION MODIFIES dataloc. Report-only, by design.

report = struct('version','', 'isLegacySchema', false, ...
    'naCounts', struct(), 'naTotal', 0, 'missingCellIndexXYs', []);

fprintf('--- dl_check_platemap_compat: %s ---\n', ...
    conditional_str(isfield(dataloc,'file') && isfield(dataloc.file,'base'), ...
                     dataloc, 'base', '(no file.base)'));

%% 1. Version / schema axis
if ~isfield(dataloc,'version') || isempty(dataloc.version)
    report.version = 'missing (legacy)';
    report.isLegacySchema = true;
else
    report.version = char(dataloc.version);
    report.isLegacySchema = ~strcmp(report.version,'5');
end
fprintf('  dataloc.version : %s%s\n', report.version, ...
    conditionalSuffix(report.isLegacySchema, '  <-- pre-v5 schema'));

%% 2. NA contamination axis (mirrors local_sanitize_na predicate exactly)
if ~isfield(dataloc,'platemapd') || isempty(dataloc.platemapd) || ~isfield(dataloc.platemapd,'pmd') || isempty(dataloc.platemapd.pmd)
    fprintf('  platemapd.pmd   : ABSENT/EMPTY -- cannot check NA contamination.\n');
else
    pmd = dataloc.platemapd.pmd;
    fn = fieldnames(pmd);
    for k = 1:numel(fn)
        val = pmd.(fn{k});
        if ~iscell(val); continue; end
        isNA = cellfun(@(v) ischar(v) && strcmp(strtrim(v),'NA'), val);
        n = sum(isNA(:));
        if n > 0
            report.naCounts.(fn{k}) = n;
            report.naTotal = report.naTotal + n;
        end
    end
    if report.naTotal == 0
        fprintf('  NA contamination: none found -- pmd looks sanitized.\n');
    else
        fprintf('  NA contamination: %d cell(s) across %d field(s):\n', ...
            report.naTotal, numel(fieldnames(report.naCounts)));
        naFields = fieldnames(report.naCounts);
        for k = 1:numel(naFields)
            fprintf('    %-12s : %d\n', naFields{k}, report.naCounts.(naFields{k}));
        end
        fprintf(['  ACTION: force a platemap re-read to sanitize, e.g.\n' ...
                 '    dataloc = dl_ensure_platemap(dataloc, struct(''updatemap'',true));\n' ...
                 '  This function does NOT do that for you.\n']);
    end
end

%% 3. cellindex axis
if isfield(dataloc,'d') && ~isempty(dataloc.d)
    missing = [];
    for ii = 1:numel(dataloc.d)
        if isempty(dataloc.d{ii}) || ~isfield(dataloc.d{ii},'data') || isempty(fieldnames(dataloc.d{ii}.data))
            continue
        end
        if ~isfield(dataloc.d{ii},'cellindex') || isempty(dataloc.d{ii}.cellindex)
            missing(end+1) = ii; %#ok<AGROW>
        end
    end
    report.missingCellIndexXYs = missing;
    if isempty(missing)
        fprintf('  cellindex       : present for all non-empty xys.\n');
    else
        fprintf('  cellindex       : MISSING for %d xy(s): %s\n', ...
            numel(missing), mat2str(missing));
        fprintf('    (plot_ functions will fabricate 1:nRows + warn per-xy at plot time)\n');
    end
else
    fprintf('  dataloc.d       : ABSENT/EMPTY -- cannot check cellindex.\n');
end

fprintf('-----------------------------------------\n');
end

function s = conditional_str(cond, dataloc, fieldName, fallback)
if cond; s = dataloc.file.(fieldName); else; s = fallback; end
end

function s = conditionalSuffix(cond, txt)
if cond; s = txt; else; s = ''; end
end