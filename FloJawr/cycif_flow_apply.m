function mask = cycif_flow_apply(T, gates, varargin)
%CYCIF_FLOW_APPLY  Evaluate gate(s) against a flow table (headless replay).
%
% mask = cycif_flow_apply(T, gates, Name, Value, ...)
%
% INPUT
%   T     : flow table from cycif_flow_build
%   gates : a single gate struct, a struct array, or a cell array of gate
%           structs (as produced by cycif_flow_gate).
%
% NAME-VALUE
%   'GateName' : name of the gate whose mask to return  (default: last gate)
%   'Return'   : 'mask' (logical N×1) | 'all' (N×G table, one column per gate)
%
% OUTPUT
%   mask : logical N×1 mask, OR an N×G table of all gate masks
%
% Notes
%   Gate coordinates are stored in TRANSFORMED space; this function applies
%   the recorded transform before testing inclusion. Parent gates are
%   evaluated first and AND-composed into their children automatically.

p = inputParser;
p.addParameter('GateName','');
p.addParameter('Return','mask', @(x) any(strcmpi(x,{'mask','all'})));
p.parse(varargin{:}); o = p.Results;

if iscell(gates); gates = [gates{:}]; end
if isempty(gates); mask = true(height(T),1); return; end

names = {gates.name};
byName = containers.Map(names, num2cell(1:numel(gates)));

N = height(T);
M = true(N, numel(gates));

% Topological order: parents before children
order = local_topo(gates, names);
for k = order(:).'
    G = gates(k);
    m = local_gate_mask(T, G);
    if ~isempty(G.parent) && byName.isKey(G.parent)
        m = m & M(:, byName(G.parent));
    end
    M(:,k) = m;
end

if strcmpi(o.Return,'all')
    mask = array2table(M, 'VariableNames', matlab.lang.makeValidName(names));
    return
end

if isempty(char(o.GateName))
    tgt = numel(gates);
else
    assert(byName.isKey(char(o.GateName)), ...
        'cycif_flow_apply: unknown gate "%s".', char(o.GateName));
    tgt = byName(char(o.GateName));
end
mask = M(:, tgt);
end

% -------------------------------------------------------------------------
function ord = local_topo(gates, names)
% Iterative topological sort: parents first. Cycles raise an error.
n   = numel(gates);
in  = zeros(n,1);                       % indegree
adj = cell(n,1);                        % adj{parent}=children
for i = 1:n
    pr = gates(i).parent;
    if ~isempty(pr)
        j = find(strcmp(names, pr), 1);
        if ~isempty(j)
            in(i) = in(i) + 1;
            adj{j}(end+1) = i; %#ok<AGROW>
        end
    end
end
ord = zeros(1,n); q = find(in==0).'; p = 0;
while ~isempty(q)
    i = q(1); q(1) = [];
    p = p + 1; ord(p) = i;
    for j = adj{i}
        in(j) = in(j) - 1;
        if in(j)==0; q(end+1) = j; end %#ok<AGROW>
    end
end
assert(p == n, 'cycif_flow_apply: gate hierarchy has a cycle.');
end

function m = local_gate_mask(T, G)
% Evaluate a single (non-hierarchical) gate: transform → shape test.
assert(ismember(G.xchan, T.Properties.VariableNames), ...
    'cycif_flow_apply: gate "%s" references missing column "%s".', G.name, G.xchan);
xt = local_trans(T.(G.xchan), G.xtrans);
yt = [];
if ~isempty(G.ychan)
    assert(ismember(G.ychan, T.Properties.VariableNames), ...
        'cycif_flow_apply: gate "%s" references missing column "%s".', G.name, G.ychan);
    yt = local_trans(T.(G.ychan), G.ytrans);
end

switch lower(G.shape)
    case 'rect'
        c = G.coords;
        m = xt >= c(1,1) & xt <= c(2,1) & yt >= c(1,2) & yt <= c(2,2);
    case 'poly'
        m = inpolygon(xt, yt, G.coords(:,1), G.coords(:,2));
    case 'thresh'
        m = xt >= G.coords(1);
    otherwise
        error('cycif_flow_apply: unknown shape ''%s''.', G.shape);
end

% Require finite values on whichever axes the gate uses
m = m & isfinite(xt);
if ~isempty(yt); m = m & isfinite(yt); end
end

function xt = local_trans(x, tr)
switch lower(tr)
    case 'linear'; xt = x;
    case 'log10';  xt = log10(max(x, eps));
    case 'asinh';  xt = asinh(x/150);
    otherwise;     error('cycif_flow_apply: unknown transform ''%s''.', tr);
end
end
