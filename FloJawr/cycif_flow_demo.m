%% CYCIF_FLOW_DEMO  End-to-end walkthrough of the CycIF flow toolkit
% -------------------------------------------------------------------------
% Runnable standalone: generates a synthetic dataloc that mirrors the
% project schema (THP1 polarization with CycIF endpoints + live sensors),
% then exercises every function in the toolkit. Swap
% `dataloc = makeSyntheticDataloc(...)` for your real load step and the
% rest of the script carries over unchanged.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ---- 0. PATH ----------------------------------------------------------
% On your machine, this is the folder holding cycif_flow_*.m
addpath('\\albecklab.mcb.ucdavis.edu\data\Code\Users\mhardy\MATLAB\plot_by_moi\FloJawr');

%% ---- 1. DATA ---------------------------------------------------------
% Replace with your real load:
%
%   dataloc  = DatalocHandler('initialize', true);
%   S        = load(fullfile(dataloc.fold.proc, dataloc.file.proc),'dataloc');
%   dataloc  = S.dataloc;
%
% Synthetic for demo purposes:
rng(42);
dataloc = makeSyntheticDataloc('nXY',24, 'nT',50, 'nCellsPerWell',400);

%% ---- 2. BUILD THE EVENT TABLE ----------------------------------------
T = cycif_flow_build(dataloc, ...
        'SensorSummaries', {'end','mean','max','auc'}, ...
        'IncludeCoords',   false);

fprintf('\nTable: %d cells × %d columns\n', height(T), width(T));
fprintf('Conditions: %s\n', strjoin(string(categories(T.pTx_Name)), ', '));
disp(head(T,3));

% Dynamic channel discovery — no hardcoded names
stainCols  = T.Properties.VariableNames(endsWith(T.Properties.VariableNames, ...
             {'_nuc','_cyt','_cell'}));
sensorCols = T.Properties.VariableNames(endsWith(T.Properties.VariableNames, ...
             {'_end','_mean','_max','_auc'}));
fprintf('Stains discovered : %s\n', strjoin(stainCols,  ', '));
fprintf('Sensors discovered: %s\n', strjoin(sensorCols, ', '));

%% ---- 3. INTERACTIVE GATE MANAGER (optional) --------------------------
% Uncomment to draw gates by hand. Gates persist to the .mat file so
% downstream sections remain reproducible.
%
gateFile = fullfile(tempdir,'demo_gates.mat');
gates    = cycif_flow_manager(T, 'GateFile', gateFile);

%% ---- 4. SCRIPTED GATES (for a runnable demo) -------------------------
% These exercise every gate shape and a two-level hierarchy. Coords are
% in TRANSFORMED axis space — same convention the manager writes.
gates = buildDemoGates();
fprintf('\nDemo gates (%d):\n', numel(gates));
for k = 1:numel(gates)
    par = gates(k).parent; if isempty(par); par = '(root)'; end
    fprintf('  %-12s  shape=%-6s  parent=%s\n', gates(k).name, gates(k).shape, par);
end

%% ---- 5. APPLY GATES --------------------------------------------------
allMasks = cycif_flow_apply(T, gates, 'Return','all');
disp(head(allMasks, 5));

%% ---- 6. STATS ---------------------------------------------------------
% Condition-level summary (pool wells with identical treatments)
S_cond = cycif_flow_stats(T, gates, ...
    'GroupBy', {'Cell','pTx_Name'}, ...
    'Chans',   {'HYLIGHT_end','HYLIGHT_mean'});
disp(S_cond);

% Per-well replicates (for error bars / statistics downstream)
S_well = cycif_flow_stats(T, gates, 'GroupBy','xy');
fprintf('Per-well stats: %d rows (= nGates × nWells)\n', height(S_well));

%% ---- 7. FIGURES -------------------------------------------------------
% 7a: biaxial with density, no grouping
figure('Name','Density scatter','Position',[100 100 640 520]);
cycif_flow_plot(T, 'iNOS_cell', 'Y','CD86_cell', ...
    'XTrans','log10','YTrans','log10', 'Density',true, 'Parent',gcf);
title('All cells — iNOS vs CD86 (density)');

% 7b: faceted by treatment, gate overlaid
figure('Name','Faceted by condition','Position',[120 120 900 600]);
cycif_flow_plot(T, 'iNOS_cell', 'Y','CD86_cell', ...
    'GroupBy','pTx_Name', 'Facet','tile', ...
    'XTrans','log10','YTrans','log10', ...
    'Gate', gates(strcmp({gates.name},'M1-like')), 'Parent',gcf);

% 7c: sensor histogram, grouped overlay
figure('Name','F16BP','Position',[140 140 640 520]);
cycif_flow_plot(T, 'HYLIGHT_end', ...
    'GroupBy','pTx_Name', 'Facet','overlay', ...
    'XTrans','linear', 'NBins',60, 'Parent',gcf);
title('ATP:ADP (HYLIGHT endpoint) by condition');

% 7d: one gate's population across conditions — bar chart from stats
figure('Name','M1-like % by condition','Position',[160 160 640 420]);
m1 = S_cond(S_cond.Gate == "M1-like", :);
bar(categorical(string(m1.pTx_Name)), m1.Pct_parent);
ylabel('% of parent (viable)'); title('M1-like % across treatments');
grid on;

%% ---- 8. HEADLESS REPLAY ----------------------------------------------
% Show that gates work detached from their original drawing context —
% save, re-load, apply to a different cell subset.
tmpGateFile = fullfile(tempdir,'demo_gates_replay.mat');
save(tmpGateFile,'gates');

L = load(tmpGateFile,'gates'); gates_reloaded = L.gates;
assert(isequal(gates, gates_reloaded), 'Gate round-trip mismatch');

% Apply only to THP1 wells (trivial here since all are THP1, but the point
% is that the same gate struct works on any event table with the same
% channel names).
mask_thp1_m1 = cycif_flow_apply(T(T.Cell=="THP1",:), gates_reloaded, ...
                                'GateName','M1-like');
fprintf('\nHeadless replay: %d M1-like cells in THP1 subset\n', nnz(mask_thp1_m1));

delete(tmpGateFile);

% =========================================================================
% LOCAL HELPERS — synthetic data + scripted gates for the demo
% =========================================================================
function dataloc = makeSyntheticDataloc(varargin)
% Minimal dataloc with the fields cycif_flow_build needs. Six conditions
% (M0, M1, M2) × two cell lines, ~400 cells each, with bimodal stain
% mixtures so the gates have something to separate.
p = inputParser;
p.addParameter('nXY', 24);
p.addParameter('nT',  50);
p.addParameter('nCellsPerWell', 400);
p.parse(varargin{:}); o = p.Results;

cellLines = ["THP1","U937"];
treatments = struct( ...
    'Name', {'Vehicle','IFNg_LPS','IL4','IFNg_LPS','IL4','Vehicle'}, ...
    'Fmix', {[0 0 1], [0.7 0.2 0.1], [0.1 0.7 0.2], ...
             [0.75 0.15 0.1], [0.05 0.75 0.2], [0 0 1]});        % M1 M2 M0

% Platemap cells — one row of wells per (cell × treatment)
pmd.xy   = cell(8,12);    pmd.xy(:) = {NaN};
pmd.Cell = cell(8,12,3);  pmd.Cell(:) = {''};
pmd.Gene = cell(8,12,3);  pmd.Gene(:) = {''};
pmd.pTx  = cell(8,12,5);  pmd.pTx(:)  = {''};
pmd.Tx1  = cell(8,12,10); pmd.Tx1(:)  = {''};

nTx = numel(treatments);
d   = cell(1, o.nXY);
xy  = 1;
for c = 1:numel(cellLines)
    for t = 1:nTx
        if xy > o.nXY; break; end
        for rep = 1:2                                 % 2 wells per condition
            if xy > o.nXY; break; end
            row = mod(xy-1,8)+1;
            col = floor((xy-1)/8)+1;
            pmd.xy{row,col}     = xy;
            pmd.Cell{row,col,1} = char(cellLines(c));
            pmd.Gene{row,col,1} = 'HYLIGHT';
            pmd.pTx{row,col,1}  = treatments(t).Name;

            d{xy} = synthesizeWell(o.nCellsPerWell, o.nT, treatments(t).Fmix);
            xy = xy + 1;
        end
    end
end

dataloc.d          = d;
dataloc.platemapd.pmd = pmd;
dataloc.movieinfo.tsamp = 6;
end

function W = synthesizeWell(N, nT, fmix)
% fmix = [f_M1 f_M2 f_M0], sums to 1. Draws each cell from a state and
% gives it CycIF + sensor values consistent with that state.
states = randsample(1:3, N, true, fmix);    % 1=M1, 2=M2, 3=M0

mu = struct('iNOS',  [3.2 2.0 2.0], ...
            'CD86',  [3.1 2.1 2.1], ...
            'CD206', [2.0 3.0 2.0], ...
            'CD163', [2.0 3.0 2.0]);

W.data = struct();
W.data.iNOS_cell  = nanpos(mu.iNOS(states)'  + 0.3*randn(N,1));
W.data.iNOS_nuc   = nanpos(W.data.iNOS_cell  - 0.3 + 0.1*randn(N,1));
W.data.CD86_cell  = nanpos(mu.CD86(states)'  + 0.3*randn(N,1));
W.data.CD206_cell = nanpos(mu.CD206(states)' + 0.3*randn(N,1));
W.data.CD163_cell = nanpos(mu.CD163(states)' + 0.3*randn(N,1));

% Drop 10% of stain measurements at random (mimic registration misses)
for f = ["iNOS_cell","iNOS_nuc","CD86_cell","CD206_cell","CD163_cell"]
    drop = rand(N,1) < 0.10;
    W.data.(f)(drop) = NaN;
end

% Live-cell sensor: ATP:ADP trajectory modulated by state
baseByState = [0.8 1.2 1.0];         % M1 lower, M2 higher
t = linspace(0, 1, nT);
drift = 0.2*sin(2*pi*t);             % shared oscillation
W.data.HYLIGHT = baseByState(states)' + drift + 0.05*randn(N,nT);

% XCoord/YCoord/nArea — optional sensors (excluded by default in build)
W.data.XCoord = randi(1600, N, nT);
W.data.YCoord = randi(1600, N, nT);
W.data.nArea  = 80 + 10*randn(N,nT);

W.cellindex = (1:N)';
end

function x = nanpos(x); x(x<0) = 0.01; end

function gates = buildDemoGates()
% Scripted gates in TRANSFORMED space (log10), mirroring what the manager
% would produce. One root rectangle, two polygon siblings, one threshold
% child — exercises every shape + one hierarchy level.

gates = struct('name',{},'shape',{},'xchan',{},'ychan',{}, ...
               'xtrans',{},'ytrans',{},'coords',{},'parent',{},'created',{});

gates(end+1) = mk('viable','rect','iNOS_cell','CD86_cell','log10','log10', ...
                  [-2 -2; 6 6], '');

gates(end+1) = mk('M1-like','poly','iNOS_cell','CD86_cell','log10','log10', ...
                  [2.4 2.4; 2.4 5; 5 5; 5 2.4], 'viable');

gates(end+1) = mk('M2-like','poly','CD206_cell','CD163_cell','log10','log10', ...
                  [2.4 2.4; 2.4 5; 5 5; 5 2.4], 'viable');

gates(end+1) = mk('M1_lowGlyc','thresh','HYLIGHT_end','','linear','linear', ...
                  1.0, 'M1-like');                  % cells with end Glycolysis < 1.0
% Note: thresh is >= by convention in cycif_flow_apply; flip sign or
% invert the column if you want the complementary side.
end

function G = mk(name, shape, xc, yc, xt, yt, coords, parent)
G = struct('name',name,'shape',shape,'xchan',xc,'ychan',yc, ...
           'xtrans',xt,'ytrans',yt,'coords',coords, ...
           'parent',parent,'created',datetime('now'));
end