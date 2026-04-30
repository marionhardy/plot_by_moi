function S = cycif_flow_stats(T, gates, varargin)
%CYCIF_FLOW_STATS  Population statistics per gate × group (FlowJo-style).
%
% S = cycif_flow_stats(T, gates, Name, Value, ...)
%
% INPUT
%   T     : flow table from cycif_flow_build
%   gates : struct array of gates (from cycif_flow_gate)
%
% NAME-VALUE
%   'GroupBy' : column name or cellstr of column names   (default 'xy')
%   'Chans'   : cellstr of columns to report per gate × group
%               (median + IQR on raw, un-transformed values)  (default {})
%
% OUTPUT
%   S : table with columns
%         <GroupBy...>, Gate, N, N_parent, Pct_parent, Pct_total,
%         <chan>_median, <chan>_iqr  (one pair per channel in Chans)

p = inputParser;
p.addParameter('GroupBy', 'xy');
p.addParameter('Chans', {});
p.parse(varargin{:}); o = p.Results;
if ischar(o.GroupBy) || isstring(o.GroupBy); o.GroupBy = cellstr(o.GroupBy); end

% All gate masks (parents composed in)
allMasks = cycif_flow_apply(T, gates, 'Return','all');
M = table2array(allMasks);

gnames  = {gates.name};
parents = {gates.parent};

% Group index
[gidx, gkey] = findgroups(T(:, o.GroupBy));

parts = cell(numel(gates),1);
for k = 1:numel(gates)
    m = M(:,k);
    % parent mask for "% of parent"
    if ~isempty(parents{k}) && any(strcmp(gnames, parents{k}))
        mp = M(:, strcmp(gnames, parents{k}));
    else
        mp = true(size(m));
    end

    Nin  = splitapply(@sum, double(m),  gidx);
    Npar = splitapply(@sum, double(mp), gidx);

    tbl = gkey;
    tbl.Gate       = repmat(string(gnames{k}), height(tbl), 1);
    tbl.N          = Nin;
    tbl.N_parent   = Npar;
    tbl.Pct_parent = 100 * Nin ./ max(Npar, 1);
    tbl.Pct_total  = 100 * Nin ./ height(T);

    for c = 1:numel(o.Chans)
        ch = o.Chans{c};
        assert(ismember(ch, T.Properties.VariableNames), ...
            'cycif_flow_stats: channel "%s" not in table.', ch);
        med  = splitapply(@(a,b) median(a(logical(b)),'omitnan'), ...
                          T.(ch), double(m), gidx);
        iqrv = splitapply(@(a,b) local_iqr(a(logical(b))),       ...
                          T.(ch), double(m), gidx);
        tbl.([ch '_median']) = med;
        tbl.([ch '_iqr'])    = iqrv;
    end
    parts{k} = tbl;
end
S = vertcat(parts{:});
end

function v = local_iqr(x)
x = x(isfinite(x));
if numel(x) < 2; v = NaN; return; end
q = quantile(x, [0.25 0.75]);
v = q(2) - q(1);
end
