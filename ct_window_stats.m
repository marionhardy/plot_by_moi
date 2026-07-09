function [cellStats, xyStats, dataTable] = ct_window_stats(dataloc, varargin)
% ct_window_stats  Extract mean sensor activity in a time window and test
%                  for differences between conditions (one factor at a time).
%
% Accepts a single dataloc OR a cell array of datalocs from separate
% experiments. When multiple datalocs are provided, extraction runs
% independently per experiment (each uses its own platemap), then data
% are pooled before stats. XY indices are offset to remain globally unique.
% An ExpID column tracks experiment identity in all output tables.
%
% Design mirrors plot_by_moi:
%   facetby  = the factor held fixed   (e.g. 'Oligo' → one test per oligo dose)
%   groupby  = the factor being compared within each facet (e.g. 'Gluc')
%
% One Kruskal-Wallis + Dunn post-hoc (BH-corrected) is run per facet per channel.
% Two parallel analyses are produced:
%   Cell-level : all cells pooled per condition  (high n, note pseudoreplication)
%   XY-level   : one mean per well               (conservative, proper experimental n)
%
% USAGE
%   % Single experiment
%   [cs, xs, dt] = ct_window_stats(dataloc, 'channel', {'RAMPKAR2'}, ...
%       'facetby', {'Oligo','Vehicle'}, 'groupby', {'Gluc'}, ...
%       'window_hr', [1 12], 'baseline', 'pre_tx_mean');
%
%   % Multiple experiments pooled
%   [cs, xs, dt] = ct_window_stats({f1Data{2}, f2Data{1}, f3Data{1}}, ...
%       'expnames', {'Exp1','Exp2','Exp3'}, ...
%       'channel',  {'RAMPKAR2','PercevalHR'}, ...
%       'facetby',  {'Oligo','Vehicle'}, 'groupby', {'Gluc'}, ...
%       'window_hr', [1 12], 'baseline', 'pre_tx_mean', ...
%       'exclude',  {'MK8722','NoIns','bad'});
%
% INPUTS
%   dataloc              - v5 dataloc struct OR cell array of v5 dataloc structs
%
% REQUIRED PARAMETER
%   window_hr            - [t_start t_end] hours relative to Tx (e.g. [1 12])
%
% OPTIONAL PARAMETERS
%   expnames             - cell of strings, one per dataloc    (default: {'Exp1',...})
%   channel              - cell of channel names               (default: {'RAMPKAR2'})
%   facetby              - cell of tokens for fixed factor     (default: {})
%   groupby              - cell of tokens for compared factor  (default: {})
%   baseline_window_hr   - [t_start t_end] pre-Tx window       (default: [-5 0])
%   baseline             - 'none' | 'first_non_nan' |
%                          'pre_tx_mean' | 'dff_first' |
%                          'dff_pre_tx'                         (default: 'pre_tx_mean')
%   exclude              - tokens to exclude from levels        (default: {})
%   alpha                - significance threshold               (default: 0.05)
%   savedir              - output folder                        (default: pwd/stats)
%   savename             - base filename without extension      (default: auto)
%   verbose              - print summary to console             (default: true)
%
% OUTPUTS
%   cellStats  - table of KW + pairwise Dunn results, cell-level
%   xyStats    - table of KW + pairwise Dunn results, XY-mean level
%   dataTable  - long table of raw per-cell extracted values (includes ExpID)

%% ---- Input parsing -------------------------------------------------------
ip = inputParser; ip.CaseSensitive = false;
addParameter(ip, 'expnames',          {},            @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'channel',           {'RAMPKAR2'},  @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'facetby',           {},            @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'groupby',           {},            @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'window_hr',         [],            @(x) isnumeric(x)&&numel(x)==2);
addParameter(ip, 'baseline_window_hr',[-5 0],        @(x) isnumeric(x)&&numel(x)==2);
addParameter(ip, 'baseline',          'pre_tx_mean', @(x) ischar(x)||isstring(x));
addParameter(ip, 'exclude',           {},            @(x) iscell(x)||ischar(x)||isstring(x));
addParameter(ip, 'alpha',             0.05,          @(x) isnumeric(x)&&isscalar(x));
addParameter(ip, 'savedir',           '',            @(x) ischar(x)||isstring(x));
addParameter(ip, 'savename',          '',            @(x) ischar(x)||isstring(x));
addParameter(ip, 'verbose',           true,          @islogical);
parse(ip, varargin{:}); p = ip.Results;

if isempty(p.window_hr)
    error('ct_window_stats: ''window_hr'' is required ([t_start t_end] hours rel. to Tx).');
end

% Wrap single dataloc in cell for uniform handling
if isstruct(dataloc), dataloc = {dataloc}; end
nExp = numel(dataloc);

% Build experiment names
expnames = cellstr(p.expnames);
if numel(expnames) < nExp
    for e = numel(expnames)+1:nExp
        expnames{e} = sprintf('Exp%d', e);
    end
end

% Normalise list inputs to row cell arrays
chans     = cellstr(p.channel);   if iscolumn(chans),     chans     = chans';     end
p.facetby = cellstr(p.facetby);   if iscolumn(p.facetby), p.facetby = p.facetby'; end
p.groupby = cellstr(p.groupby);   if iscolumn(p.groupby), p.groupby = p.groupby'; end
p.exclude = cellstr(p.exclude);   if iscolumn(p.exclude), p.exclude = p.exclude'; end
baseMode  = lower(strtrim(char(p.baseline)));

%% ---- Pre-allocate output storage ----------------------------------------
rowAccum = {};   % accumulates one row per cell; converted to table after loop

%% ---- Per-experiment extraction loop -------------------------------------
xyOffset = 0;   % cumulative XY offset to keep indices globally unique

for iExp = 1:nExp
    dl       = dataloc{iExp};
    expID    = expnames{iExp};

    % Timing (per experiment)
    if isfield(dl,'movieinfo') && ~isempty(dl.movieinfo.tsamp)
        looptime = dl.movieinfo.tsamp;
    else
        looptime = 6;
        warning('ct_window_stats [%s]: movieinfo.tsamp missing, assuming 6 min/frame.', expID);
    end
    dt_hr = looptime / 60;

    % Resolve platemap levels (per experiment — handles different layouts)
    pmd = dl.platemapd.pmd;

    facetFields = cws_resolve_fields(pmd, p.facetby, {'pTx','Tx1','Tx2'});
    groupFields = cws_resolve_fields(pmd, p.groupby, {'Tx1','Tx2','pTx'});

    facetLevels = cws_find_levels(pmd, p.facetby, facetFields, dt_hr);
    if isempty(facetLevels)
        facetLevels = cws_find_levels(pmd, {''}, facetFields, dt_hr);
    end
    groupLevels = cws_find_levels(pmd, p.groupby, groupFields, dt_hr);
    if isempty(groupLevels)
        groupLevels = cws_find_levels(pmd, {''}, groupFields, dt_hr);
    end

    % Apply exclusions
    if ~isempty(p.exclude)
        exL = lower(p.exclude);
        facetLevels = facetLevels(~arrayfun(@(x) any(contains(lower(x.title),exL)), facetLevels));
        groupLevels = groupLevels(~arrayfun(@(x) any(contains(lower(x.title),exL)), groupLevels));
    end
    if isempty(facetLevels)
        warning('ct_window_stats [%s]: no facet levels found — skipping.', expID);
        xyOffset = xyOffset + numel(dl.d);
        continue;
    end

    % Global Tx time for this experiment
    allTimes = [[facetLevels.time_hr], [groupLevels.time_hr]];
    allTimes = allTimes(isfinite(allTimes) & allTimes > 0);
    tx_hr    = min([allTimes, NaN]);
    if ~isfinite(tx_hr)
        warning('ct_window_stats [%s]: no treatment time found — anchoring to frame 1.', expID);
        tx_hr = 0;
    end

    % Extraction
    for ch = 1:numel(chans)
        chan = chans{ch};

        for f = 1:numel(facetLevels)
            facetXY    = facetLevels(f).xys(:)';
            facetTitle = facetLevels(f).title;

            % Build active groups for this facet
            activeGroups = struct('title',{},'xys',{});
            for g = 1:numel(groupLevels)
                xyg = intersect(facetXY, groupLevels(g).xys(:)');
                if ~isempty(xyg)
                    activeGroups(end+1).title = groupLevels(g).title; %#ok<AGROW>
                    activeGroups(end).xys     = xyg;
                end
            end
            if isempty(activeGroups)
                activeGroups(1).title = facetTitle;
                activeGroups(1).xys   = facetXY;
            end

            for g = 1:numel(activeGroups)
                xys = activeGroups(g).xys(:)';

                for xy = xys
                    if xy < 1 || xy > numel(dl.d), continue; end
                    S = dl.d{xy};
                    if isempty(S)||~isstruct(S)||~isfield(S,'data'), continue; end
                    if ~isfield(S.data, chan), continue; end

                    M  = S.data.(chan);
                    ci = cws_get_cellindex(S, size(M,1));
                    if isempty(M) || size(M,1) ~= numel(ci), continue; end

                    T   = size(M,2);
                    txF = max(1, min(T, round(tx_hr / dt_hr) + 1));

                    M = cws_apply_baseline(M, baseMode, txF, ...
                        p.baseline_window_hr, dt_hr, T);

                    wF1 = max(1,   min(T, txF + round(p.window_hr(1) / dt_hr)));
                    wF2 = max(wF1, min(T, txF + round(p.window_hr(2) / dt_hr)));

                    cellMu = mean(M(:, wF1:wF2), 2, 'omitnan');
                    valid  = isfinite(cellMu);
                    cellMu = cellMu(valid);
                    ciVal  = ci(valid);
                    if isempty(cellMu), continue; end

                    % XY index offset makes it globally unique across experiments
                    xyGlobal = xy + xyOffset;

                    for r = 1:numel(cellMu)
                        rowAccum{end+1} = {expID, chan, facetTitle, ...
                            activeGroups(g).title, xyGlobal, ciVal(r), cellMu(r)}; %#ok<AGROW>
                    end
                end % xy
            end % group
        end % facet
    end % channel

    % Advance XY offset by the number of XY slots in this experiment
    xyOffset = xyOffset + numel(dl.d);

end % experiment loop

if isempty(rowAccum)
    error('ct_window_stats: no data extracted. Check channel names, window, and XY filters.');
end

%% ---- Build long data table ----------------------------------------------
dataTable = cell2table(vertcat(rowAccum{:}), ...
    'VariableNames', {'ExpID','Channel','Facet','Group','XY','CellID','WindowMean'});
dataTable.ExpID      = string(dataTable.ExpID);
dataTable.Channel    = string(dataTable.Channel);
dataTable.Facet      = string(dataTable.Facet);
dataTable.Group      = string(dataTable.Group);
if iscell(dataTable.XY),         dataTable.XY         = cell2mat(dataTable.XY);         end
if iscell(dataTable.CellID),     dataTable.CellID     = cell2mat(dataTable.CellID);     end
if iscell(dataTable.WindowMean), dataTable.WindowMean = cell2mat(dataTable.WindowMean); end

%% ---- Statistical tests --------------------------------------------------
% One KW + Dunn per channel × facet combination.
% Cell-level  : every cell is an observation (high n — note pseudoreplication).
% Replicate-level: one mean per experiment per condition (n = number of experiments).

cellStatRows = {};
repStatRows  = {};

chanList  = unique(dataTable.Channel, 'stable');
facetList = unique(dataTable.Facet,   'stable');

for ch = 1:numel(chanList)
    for f = 1:numel(facetList)
        mask = dataTable.Channel == chanList(ch) & ...
               dataTable.Facet   == facetList(f);
        sub  = dataTable(mask, :);
        if isempty(sub), continue; end

        groups = unique(sub.Group, 'stable');
        nG     = numel(groups);
        if nG < 2
            if p.verbose
                fprintf('[%s | %s] Only 1 group — skipping stats.\n', ...
                    chanList(ch), facetList(f));
            end
            continue;
        end

        % ---- Cell-level KW -----------------------------------------------
        cellStatRows = cws_run_kw_dunn(sub, groups, ...
            chanList(ch), facetList(f), 'cell', p.alpha, cellStatRows, p.verbose);

        % ---- Replicate-level: one-sample t-test vs 0 per group ----------
        % Collapse to one mean per ExpID per condition, then test vs 0.
        repMu = varfun(@(x) mean(x,'omitnan'), sub, ...
            'GroupingVariables', {'ExpID','Group'}, ...
            'InputVariables',    'WindowMean');
        repMu.Properties.VariableNames{'Fun_WindowMean'} = 'WindowMean';

        repStatRows = cws_run_ttest_vs_baseline(repMu, groups, ...
            chanList(ch), facetList(f), p.alpha, repStatRows, p.verbose);

    end
end

%% ---- Assemble stats tables ----------------------------------------------
statVars = {'Channel','Facet','Level','Group1','Group2', ...
            'n1','n2','KW_p','Dunn_z','Dunn_p_raw','Dunn_p_BH','Significant'};

if ~isempty(cellStatRows)
    cellStats = cell2table(vertcat(cellStatRows{:}), 'VariableNames', statVars);
else
    cellStats = cell2table({}, 'VariableNames', statVars);
    warning('ct_window_stats: no cell-level stats produced.');
end

if ~isempty(repStatRows)
    xyStats = cell2table(vertcat(repStatRows{:}), ...
        'VariableNames', {'Channel','Facet','Group','n','Mean','SEM','t','df','p','Significant'});
else
    xyStats = cell2table({}, ...
        'VariableNames', {'Channel','Facet','Group','n','Mean','SEM','t','df','p','Significant'});
    warning('ct_window_stats: no replicate-level stats produced.');
end

%% ---- Save to Excel ------------------------------------------------------
saveDir = char(p.savedir);
if isempty(saveDir), saveDir = fullfile(pwd,'stats'); end
if ~exist(saveDir,'dir'), mkdir(saveDir); end

saveName = char(p.savename);
if isempty(saveName)
    chanStr  = regexprep(strjoin(chans,'_'),    '[^a-zA-Z0-9_]','_');
    winStr   = sprintf('win%.1f_%.1fhr', p.window_hr(1), p.window_hr(2));
    saveName = sprintf('ct_window_stats_%s_%s', chanStr, winStr);
end
saveName = regexprep(saveName, '[^a-zA-Z0-9_\-]','_');
xlsPath  = fullfile(saveDir, [saveName '.xlsx']);

writetable(dataTable,  xlsPath, 'Sheet','RawData');
writetable(cellStats,  xlsPath, 'Sheet','Stats_CellLevel');
writetable(xyStats,    xlsPath, 'Sheet','Stats_ReplicateLevel');
fprintf('ct_window_stats: saved → %s\n', xlsPath);

end % main function

%% =========================================================================
%% Local helpers
%% =========================================================================

function [statRows] = cws_run_kw_dunn(tbl, groups, chan, facet, level, alpha, statRows, verbose)
% Run KW on all groups, then pairwise Dunn with BH correction within this facet.
    nG = numel(groups);

    % Collect group vectors
    gVecs = cell(nG,1);
    for g = 1:nG
        gVecs{g} = tbl.WindowMean(tbl.Group == groups(g));
    end

    % KW on all groups together
    allVals  = cellfun(@(x) x(:), gVecs, 'UniformOutput', false);
    allVals  = vertcat(allVals{:});
    allLabels = [];
    for g = 1:nG
        allLabels = [allLabels; repmat(g, numel(gVecs{g}), 1)]; %#ok<AGROW>
    end
    [kw_p, ~, kw_stats] = kruskalwallis(allVals, allLabels, 'off');

    % Pairwise Dunn test (manual: z-statistic from KW rank sums)
    nPairs   = nG*(nG-1)/2;
    p_raw    = nan(nPairs,1);
    z_stat   = nan(nPairs,1);
    pairIdx  = nan(nPairs,2);
    idx = 0;
    N = numel(allVals);

    % Rank all values together (for Dunn)
    [~, rankAll] = sort(allVals);
    ranks = zeros(N,1); ranks(rankAll) = 1:N;

    % Handle ties: assign average ranks
    [uVals, ~, ic] = unique(allVals);
    for u = 1:numel(uVals)
        tied = find(ic == u);
        if numel(tied) > 1, ranks(tied) = mean(ranks(tied)); end
    end

    % Group rank sums and sizes
    Ri = zeros(nG,1); ni = zeros(nG,1);
    offset = 0;
    for g = 1:nG
        ni(g) = numel(gVecs{g});
        Ri(g) = sum(ranks(offset+1 : offset+ni(g)));
        offset = offset + ni(g);
    end

    % Tie correction factor
    [~, ~, ic2] = unique(allVals);
    tCounts = histcounts(ic2, 1:max(ic2)+1);
    tieFactor = sum(tCounts.^3 - tCounts) / (12*(N-1));

    for i = 1:nG-1
        for j = i+1:nG
            idx = idx + 1;
            pairIdx(idx,:) = [i j];
            SE = sqrt( (N*(N+1)/12 - tieFactor) * (1/ni(i) + 1/ni(j)) );
            z_stat(idx) = (Ri(i)/ni(i) - Ri(j)/ni(j)) / SE;
            p_raw(idx)  = 2 * (1 - normcdf(abs(z_stat(idx))));
        end
    end

    % BH correction within this facet×channel block
    p_bh = cws_bh_correct(p_raw);

    % Accumulate rows
    for idx = 1:nPairs
        i = pairIdx(idx,1); j = pairIdx(idx,2);
        sig = p_bh(idx) < alpha;
        statRows{end+1} = {char(chan), char(facet), level, ...
            char(groups(i)), char(groups(j)), ...
            ni(i), ni(j), ...
            kw_p, z_stat(idx), p_raw(idx), p_bh(idx), sig}; %#ok<AGROW>
    end

    if verbose
        fprintf('\n=== KW: %s | %s | %s-level  (p=%.4f) ===\n', ...
            char(chan), char(facet), level, kw_p);
        for idx = 1:nPairs
            i = pairIdx(idx,1); j = pairIdx(idx,2);
            sig = p_bh(idx) < alpha;
            fprintf('  %-20s vs %-20s  z=%6.3f  p_raw=%.4f  p_BH=%.4f  %s\n', ...
                char(groups(i)), char(groups(j)), ...
                z_stat(idx), p_raw(idx), p_bh(idx), ...
                repmat('*', 1, sig));
        end
    end
end

function p_bh = cws_bh_correct(p_raw)
% Benjamini-Hochberg FDR correction.
    m = numel(p_raw);
    [ps, I] = sort(p_raw(:));
    p_bh_sorted = min(1, ps .* m ./ (1:m)');
    % Enforce monotonicity (from largest to smallest)
    for k = m-1:-1:1
        p_bh_sorted(k) = min(p_bh_sorted(k), p_bh_sorted(k+1));
    end
    p_bh = nan(size(p_raw));
    p_bh(I) = p_bh_sorted;
end

function M = cws_apply_baseline(M, mode, txF, bwin_hr, dt_hr, T)
% Apply baseline correction identical to plot_by_moi modes.
    switch mode
        case 'none'
            return;

        case 'first_non_nan'
            % First non-NaN value per row
            b = nan(size(M,1),1);
            for r = 1:size(M,1)
                j = find(~isnan(M(r,:)),1,'first');
                if ~isempty(j), b(r) = M(r,j); end
            end
            M = M - b;

        case 'pre_tx_mean'
            % Mean over [txF + bwin_hr(1), txF + bwin_hr(2)]
            b1 = max(1,   min(T, txF + round(bwin_hr(1)/dt_hr)));
            b2 = max(b1,  min(T, txF + round(bwin_hr(2)/dt_hr)));
            b  = mean(M(:, b1:b2), 2, 'omitnan');
            M  = M - b;

        case 'first_tp'
            M = M - M(:,1);

        case 'dff_first'
            b = M(:,1);
            M = (M - b) ./ b;

        case 'dff_pre_tx'
            b1 = max(1,   min(T, txF + round(bwin_hr(1)/dt_hr)));
            b2 = max(b1,  min(T, txF + round(bwin_hr(2)/dt_hr)));
            b  = mean(M(:, b1:b2), 2, 'omitnan');
            M  = (M - b) ./ b;

        otherwise
            warning('ct_window_stats: unknown baseline mode ''%s'' — skipping.', mode);
    end
end

function ci = cws_get_cellindex(S, nRows)
% Return cellindex as column vector; fabricate 1:nRows if absent.
    if isfield(S,'cellindex') && numel(S.cellindex) == nRows
        ci = S.cellindex(:);
    elseif isfield(S,'cellindex') && ~isempty(S.cellindex)
        warning('ct_window_stats: cellindex length mismatch — using row indices.');
        ci = (1:nRows)';
    else
        ci = (1:nRows)';
    end
end

function fields = cws_resolve_fields(pmd, requested, defaults)
    fields = {};
    for i = 1:numel(requested)
        if ~isempty(requested{i}) && isfield(pmd, requested{i})
            fields{end+1} = requested{i}; %#ok<AGROW>
        end
    end
    if isempty(fields)
        for i = 1:numel(defaults)
            if isfield(pmd, defaults{i})
                fields{end+1} = defaults{i}; %#ok<AGROW>
            end
        end
    end
end

function levels = cws_find_levels(pmd, tokens, fields, dt_hr)
% Same logic as pss_find_levels in plot_sorted_stacks.
    tokens = cellstr(tokens);
    [R,C]  = size(pmd.xy);
    gMap   = containers.Map();

    for t = 1:numel(tokens)
        tok = tokens{t}; if isempty(tok), continue; end
        for f = 1:numel(fields)
            fname = fields{f}; if ~isfield(pmd,fname), continue; end
            A = pmd.(fname);
            for r = 1:R; for c = 1:C
                nm = cws_safe_str(A,r,c,1);
                if isempty(nm)||~contains(lower(nm),lower(tok)), continue; end
                dose = cws_safe_num(A,r,c,2);
                unit = cws_safe_str(A,r,c,3);
                time_hr = NaN;
                if any(strcmp(fname,{'Tx1','Tx2','Tx3','Tx4'}))
                    t_raw  = cws_safe_num(A,r,c,4);
                    t_unit = lower(strtrim(cws_safe_str(A,r,c,5)));
                    if isfinite(t_raw)
                        switch t_unit
                            case {'h','hr','hrs','hour','hours'}, time_hr = t_raw;
                            otherwise,                            time_hr = t_raw * dt_hr;
                        end
                    end
                end
                if isnan(dose), key=sprintf('%s_NaN_%s',nm,unit); dispNm=nm;
                else, key=sprintf('%s_%.6f_%s',nm,dose,unit);
                      dispNm=sprintf('%s %.4g %s',nm,dose,unit);
                end
                xy_here = pmd.xy{r,c}; if isempty(xy_here), continue; end
                if isKey(gMap,key)
                    e=gMap(key); e.xys=[e.xys,xy_here(:)'];
                    if isfinite(time_hr), e.times(end+1)=time_hr; end
                    gMap(key)=e;
                else
                    e.title=dispNm; e.xys=xy_here(:)'; e.dose=dose; e.times=[];
                    if isfinite(time_hr), e.times=time_hr; end
                    gMap(key)=e;
                end
            end; end
        end
    end

    keys = gMap.keys;
    if isempty(keys)
        levels = struct('title',{},'xys',{},'dose',{},'time_hr',{}); return;
    end
    levels = repmat(struct('title','','xys',[],'dose',NaN,'time_hr',NaN),1,numel(keys));
    for k = 1:numel(keys)
        e = gMap(keys{k});
        levels(k).title   = e.title;
        levels(k).xys     = unique(e.xys);
        levels(k).dose    = e.dose;
        levels(k).time_hr = min([e.times, NaN]);
    end
    doses = [levels.dose]; isNaN = isnan(doses);
    [~,I] = sort(doses(~isNaN));
    order = [find(~isNaN(I)), find(isNaN)]; %#ok<FNDSB>
    if numel(order)==numel(levels), levels=levels(order); end
end

function statRows = cws_run_ttest_vs_baseline(repTbl, groups, chan, facet, alpha, statRows, verbose)
% One-sample t-test vs 0 (baseline) for each group independently.
% repTbl has columns: ExpID, Group, WindowMean — one row per replicate.
% Reports: mean, SEM, t-statistic, df, p-value, significance flag.

    if verbose
        fprintf('\n=== t-test vs baseline: %s | %s ===\n', char(chan), char(facet));
    end

    for g = 1:numel(groups)
        vals = repTbl.WindowMean(repTbl.Group == groups(g));
        vals = vals(isfinite(vals));
        n    = numel(vals);

        if n < 2
            warning('ct_window_stats: %s | %s | %s — fewer than 2 replicates, skipping t-test.', ...
                char(chan), char(facet), char(groups(g)));
            continue;
        end

        mu  = mean(vals);
        sem = std(vals) / sqrt(n);
        df  = n - 1;

        % One-sample t-test vs 0
        [~, pval, ~, tstat] = ttest(vals, 0, 'Alpha', alpha);
        t   = tstat.tstat;
        sig = pval < alpha;

        statRows{end+1} = {char(chan), char(facet), char(groups(g)), ...
            n, mu, sem, t, df, pval, sig}; %#ok<AGROW>

        if verbose
            fprintf('  %-25s  n=%d  mean=%7.4f  SEM=%6.4f  t(%d)=%6.3f  p=%.4f  %s\n', ...
                char(groups(g)), n, mu, sem, df, t, pval, repmat('*',1,sig));
        end
    end
end

function s = cws_safe_str(A,r,c,k)
    s=''; if k>size(A,3), return; end
    try; v=A{r,c,k};
        if ischar(v),   s=strtrim(v);       return; end
        if isstring(v), s=strtrim(char(v));  return; end
    catch; end
end

function x = cws_safe_num(A,r,c,k)
    x=NaN; if k>size(A,3), return; end
    try; v=A{r,c,k};
        if isnumeric(v)&&isscalar(v), x=double(v); return; end
        if ischar(v)||isstring(v)
            vv=str2double(strtrim(char(v))); if isfinite(vv), x=vv; end
        end
    catch; end
end