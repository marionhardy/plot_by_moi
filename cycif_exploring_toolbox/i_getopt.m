% =========================================================================
% Shared private helpers (paste at bottom of each file, or into a +cycif pkg)
% =========================================================================

function v = i_getopt(opts, field, default)
    if isfield(opts, field), v = opts.(field); else, v = default; end
end

function markers = i_auto_markers(T)
    meta_cols = {'xy','ptx_name','ptx_conc','ptx_units','ptx_time','ptx_timeunit', ...
                 'ptx_label','tx1_name','tx1_conc','tx1_units','tx1_time', ...
                 'tx1_timeunit','tx2_name','tx2_conc','tx2_units','tx2_time', ...
                 'tx2_timeunit','tx1_label'};
    markers = setdiff(T.Properties.VariableNames, meta_cols, 'stable');
end