function [T, info] = cycif_filter_cells(T, group_by, exclude_spike, drop_ctrl)
% CYCIF_FILTER_CELLS
% Shared pre-plot cell filtering for the cycif_plot_* family, so spike/ctrl
% handling is defined ONCE and behaves identically across embedding, scatter,
% and distribution plots.
%
% INPUTS
%   T             table   cell table (from cycif_build_table)
%   group_by      char    the grouping column (used for the ctrl-name check)
%   exclude_spike logical drop ptx_is_spike wells. DEFAULT false -- for designs
%                         where Ptx is a core condition (e.g. glucose starve)
%                         the spike flag is not a disposable control. Guarded:
%                         if ALL cells are flagged spike, they are kept (with a
%                         warning) rather than emptying the table.
%   drop_ctrl     logical drop rows whose group_by value contains 'ctrl'
%                         (case-insensitive). DEFAULT false.
%
% OUTPUTS
%   T     table   filtered table
%   info  struct  .n_in .n_spike_dropped .n_ctrl_dropped .n_out
%
% Both filters are no-ops when the relevant column is absent, so callers can
% pass any group_by without guarding.

if nargin < 3 || isempty(exclude_spike); exclude_spike = false; end
if nargin < 4 || isempty(drop_ctrl);     drop_ctrl     = false; end

info = struct('n_in', height(T), 'n_spike_dropped', 0, ...
              'n_ctrl_dropped', 0, 'n_out', height(T));

% --- Spike wells ---
if exclude_spike && ismember('ptx_is_spike', T.Properties.VariableNames)
    is_spike = logical(T.ptx_is_spike);
    n_spike  = sum(is_spike);
    if n_spike == height(T)
        warning(['cycif_filter_cells: ALL %d cells are flagged ptx_is_spike; ' ...
                 'exclude_spike would drop everything. Keeping all cells. If Ptx ' ...
                 'is a core condition here, this flag is not meaningful for you.'], n_spike);
    elseif n_spike > 0
        T = T(~is_spike, :);
        info.n_spike_dropped = n_spike;
        fprintf('cycif_filter_cells: excluded %d spike-treatment cells.\n', n_spike);
    end
end

% --- Control treatments (group_by contains 'ctrl') ---
if drop_ctrl && ~isempty(group_by) && ismember(group_by, T.Properties.VariableNames)
    gv      = string(T.(group_by));
    is_ctrl = contains(lower(gv), 'ctrl');
    n_ctrl  = sum(is_ctrl);
    if n_ctrl == height(T)
        warning(['cycif_filter_cells: ALL cells have ''ctrl'' in %s; drop_ctrl ' ...
                 'would empty the table. Keeping all cells.'], group_by);
    elseif n_ctrl > 0
        dropped_labels = unique(gv(is_ctrl));
        T = T(~is_ctrl, :);
        info.n_ctrl_dropped = n_ctrl;
        fprintf('cycif_filter_cells: dropped %d cells in ctrl groups: %s\n', ...
                n_ctrl, strjoin(cellstr(dropped_labels), ', '));
    end
end

info.n_out = height(T);
end