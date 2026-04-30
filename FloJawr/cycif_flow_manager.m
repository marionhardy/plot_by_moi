function gates = cycif_flow_manager(T, varargin)
%CYCIF_FLOW_MANAGER  Interactive gate manager with per-gate hierarchy windows.
%
% gates = cycif_flow_manager(T, Name, Value, ...)
%
% INPUT
%   T        : flow table from cycif_flow_build
%
% NAME-VALUE
%   'GateFile'     : .mat path to load/save gates from/to   (default '')
%   'InitialGates' : struct array to seed the session       (default [])
%
% OUTPUT
%   gates : final gate struct array when the window is closed
%
% EACH GATE HAS ITS OWN PERSISTENT FIGURE showing the parent-subset
% population with the gate boundary overlaid. Figures are created on Plot
% and refreshed on Save. Editing a gate refocuses its figure; deleting it
% closes it.
%
% DRAW WORKFLOW (fixes the zoom/draw conflict):
%   1. Fill the editor form.
%   2. Click [Plot]          → figure opens (or focuses); use the axes
%                              toolbar (hover over plot) to zoom / pan to
%                              the region of interest.
%   3. Click [Draw shape]    → polygon/rect/threshold draw starts on that
%                              same axes; double-click to finish.
%   4. Click [Save gate]     → commits coords; figure stays open as the
%                              hierarchy view for that gate.

% ---------- parse ---------------------------------------------------------
p = inputParser;
p.addParameter('GateFile','', @(x) ischar(x) || isstring(x));
p.addParameter('InitialGates', [], @(x) isempty(x) || isstruct(x));
p.parse(varargin{:}); o = p.Results;

% ---------- state ---------------------------------------------------------
state.T        = T;
state.gateFile = char(o.GateFile);
state.dirty    = false;

if ~isempty(o.InitialGates)
    state.gates = o.InitialGates;
elseif ~isempty(state.gateFile) && exist(state.gateFile,'file')
    L = load(state.gateFile, 'gates');
    state.gates = L.gates;
else
    state.gates = local_empty_gate_array();
end

% Per-gate persistent windows (keyed by gate name; '__draft__' until save).
state.gateWindows   = containers.Map('KeyType','char','ValueType','any');
state.windowCounter = 0;            % for cascading positions
state.draftKey      = '__draft__';

% numeric channel inventory (exclude metadata)
vn = T.Properties.VariableNames;
isNum = false(1,numel(vn));
for i = 1:numel(vn); isNum(i) = isnumeric(T.(vn{i})); end
metaCols = {'cellID','xy','trackID','col','pTx_Conc'};
numCols  = setdiff(vn(isNum), metaCols, 'stable');

% session scratch
editingIdx    = [];   % index in state.gates; [] means 'new'
pendingCoords = [];   % coords from the most recent Draw (transformed space)

% ---------- layout --------------------------------------------------------
fig = uifigure('Name','CycIF Flow Manager', 'Position',[100 100 980 580]);
fig.CloseRequestFcn = @onClose;

root = uigridlayout(fig, [2 2]);
root.RowHeight    = {'1x', 40};
root.ColumnWidth  = {'1.2x','2.2x'};

% Left: tree + tree-level buttons
leftPanel = uipanel(root, 'Title','Gate hierarchy');
leftPanel.Layout.Row = 1; leftPanel.Layout.Column = 1;
lg = uigridlayout(leftPanel,[2 1]); lg.RowHeight = {'1x', 34};
tree = uitree(lg);
treeBtns = uigridlayout(lg,[1 2]); treeBtns.Padding = [0 0 0 0];
bNew = uibutton(treeBtns, 'Text','+ New gate');
bDel = uibutton(treeBtns, 'Text','− Delete');

% Right: editor
rightPanel = uipanel(root, 'Title','Gate editor');
rightPanel.Layout.Row = 1; rightPanel.Layout.Column = 2;
eg = uigridlayout(rightPanel, [11 2]);
eg.RowHeight   = repmat({28},1,11);
eg.ColumnWidth = {110, '1x'};

uilabel(eg,'Text','Name');          nameBox  = uieditfield(eg,'text');
uilabel(eg,'Text','Parent');        parentDD = uidropdown(eg,'Items',{'(root)'});
uilabel(eg,'Text','Shape');         shapeDD  = uidropdown(eg,'Items',{'poly','rect','thresh'});
uilabel(eg,'Text','X channel');     xDD      = uidropdown(eg,'Items',numCols);
uilabel(eg,'Text','X transform');   xtDD     = uidropdown(eg,'Items',{'log10','linear','asinh'});
uilabel(eg,'Text','Y channel');     yDD      = uidropdown(eg,'Items',[{'(none)'}, numCols(:)']);
uilabel(eg,'Text','Y transform');   ytDD     = uidropdown(eg,'Items',{'log10','linear','asinh'});

% action row: Plot / Draw / Save gate
actRow = uigridlayout(eg,[1 3]); actRow.Padding = [0 0 0 0];
actRow.Layout.Row = 8; actRow.Layout.Column = [1 2];
bPlot   = uibutton(actRow, 'Text','Plot');
bDraw   = uibutton(actRow, 'Text','Draw shape');
bSaveG  = uibutton(actRow, 'Text','Save gate');

% hint + stats
uilabel(eg,'Text','Hint'); ...
    uilabel(eg,'Text','Plot → zoom → Draw shape (right-click to place, Enter to finish)', ...
               'FontColor',[.35 .35 .35]);
uilabel(eg,'Text','Stats');
statsLbl = uilabel(eg,'Text','—','FontColor',[.2 .2 .2]);

% Bottom row: file actions
bottom = uigridlayout(root,[1 3]);
bottom.Layout.Row = 2; bottom.Layout.Column = [1 2];
bottom.ColumnWidth = {'1x','1x','1x'};
bLoad  = uibutton(bottom,'Text','Load .mat...');
bSaveF = uibutton(bottom,'Text','Save to file');
bClose = uibutton(bottom,'Text','Close');

% ---------- callbacks -----------------------------------------------------
tree.SelectionChangedFcn = @onTreeSelect;
bNew.ButtonPushedFcn     = @onNew;
bDel.ButtonPushedFcn     = @onDelete;
bPlot.ButtonPushedFcn    = @onPlot;
bDraw.ButtonPushedFcn    = @onDrawShape;
bSaveG.ButtonPushedFcn   = @onSaveGate;
bLoad.ButtonPushedFcn    = @onLoadFile;
bSaveF.ButtonPushedFcn   = @onSaveFile;
bClose.ButtonPushedFcn   = @onClose;

% ---------- initial paint -------------------------------------------------
refreshTree();
refreshParentDD();

% block until the window closes
uiwait(fig);

gates = state.gates;
closeAllGateWindows();
if isvalid(fig); delete(fig); end

% =========================================================================
% NESTED CALLBACKS
% =========================================================================
    function refreshTree()
        delete(tree.Children);
        rootNode = uitreenode(tree,'Text','(root)','NodeData',[]);
        if isempty(state.gates); expand(tree,'all'); return; end
        nodeMap = containers.Map('KeyType','char','ValueType','any');
        placed  = false(1, numel(state.gates));
        names   = {state.gates.name};
        while any(~placed)
            progressed = false;
            for k = find(~placed)
                par = state.gates(k).parent;
                if isempty(par) || ~any(strcmp(names, par))
                    anc = rootNode;
                elseif isKey(nodeMap, par)
                    anc = nodeMap(par);
                else
                    continue
                end
                nd = uitreenode(anc, ...
                    'Text', sprintf('%s  [%s]', state.gates(k).name, state.gates(k).shape), ...
                    'NodeData', k);
                nodeMap(state.gates(k).name) = nd;
                placed(k) = true;
                progressed = true;
            end
            if ~progressed; break; end
        end
        expand(tree,'all');
    end

    function refreshParentDD()
        items = {'(root)'};
        if ~isempty(state.gates)
            items = [items, {state.gates.name}];
        end
        if ~isempty(editingIdx)
            bad = descendants_of(editingIdx);
            bad = [bad, editingIdx];
            items = setdiff(items, {state.gates(bad).name}, 'stable');
        end
        prev = parentDD.Value;
        parentDD.Items = items;
        if ismember(prev, items); parentDD.Value = prev;
        else;                     parentDD.Value = '(root)'; end
    end

    function onTreeSelect(~,~)
        nd = tree.SelectedNodes;
        if isempty(nd) || isempty(nd.NodeData); clearEditor(); return; end
        k = nd.NodeData;
        editingIdx    = k;
        pendingCoords = state.gates(k).coords;
        loadGateIntoEditor(state.gates(k));
        refreshParentDD();
        updateStats();
        % focus the existing window if it's open; don't open a new one here
        nm = state.gates(k).name;
        if isKey(state.gateWindows,nm) && isvalid(state.gateWindows(nm))
            figure(state.gateWindows(nm));
        end
    end

    function loadGateIntoEditor(G)
        nameBox.Value  = G.name;
        shapeDD.Value  = G.shape;
        xDD.Value      = G.xchan;
        yDD.Value      = tern(isempty(G.ychan), '(none)', G.ychan);
        xtDD.Value     = G.xtrans;
        ytDD.Value     = G.ytrans;
        parentDD.Value = tern(isempty(G.parent), '(root)', G.parent);
    end

    function clearEditor()
        editingIdx    = [];
        pendingCoords = [];
        nameBox.Value = '';
        statsLbl.Text = '—';
        refreshParentDD();
    end

    function onNew(~,~)
        tree.SelectedNodes = [];
        clearEditor();
    end

    function onDelete(~,~)
        if isempty(editingIdx); return; end
        G = state.gates(editingIdx);
        kids = find(strcmp({state.gates.parent}, G.name));
        if ~isempty(kids)
            parentLabel = tern(isempty(G.parent),'(root)',G.parent);
            choice = uiconfirm(fig, ...
                sprintf('"%s" has %d descendant gate(s). Choose:', G.name, numel(kids)), ...
                'Delete gate', ...
                'Options',{['Reparent to ' parentLabel],'Delete descendants','Cancel'}, ...
                'DefaultOption',1, 'CancelOption',3);
            if strcmp(choice,'Cancel'); return; end
            if startsWith(choice,'Reparent')
                for c = kids; state.gates(c).parent = G.parent; end
                closeGateWindow(G.name);
            else
                todo = [editingIdx, descendants_of(editingIdx)];
                for k = todo; closeGateWindow(state.gates(k).name); end
                state.gates(todo) = [];
                postDelete(); return;
            end
        else
            closeGateWindow(G.name);
        end
        state.gates(editingIdx) = [];
        postDelete();
    end

    function postDelete()
        state.dirty = true;
        editingIdx  = [];
        refreshTree(); clearEditor();
    end

    function idxs = descendants_of(k)
        idxs = [];
        if isempty(state.gates) || k > numel(state.gates); return; end
        stack = k;
        while ~isempty(stack)
            cur = stack(end); stack(end) = [];
            kids = find(strcmp({state.gates.parent}, state.gates(cur).name));
            idxs  = [idxs, kids];                                 %#ok<AGROW>
            stack = [stack, kids];                                %#ok<AGROW>
        end
    end

    % -- PLOT: open/refresh the gate window, no drawing yet ---------------
    function onPlot(~,~)
        [xc,yc,pn] = readEditorChans();
        if isempty(xc); uialert(fig,'Pick an X channel first.','X channel'); return; end
        sub   = parentMask(pn);
        wkey  = currentGateKey();
        wname = nonempty(nameBox.Value, '(draft)');

        [gf, ax] = ensureGateWindow(wkey, wname);
        clf(gf); delete(gf.Children);
        if isempty(yc)
            cycif_flow_plot(state.T, xc, 'XTrans',xtDD.Value, 'Sub',sub, 'Parent',gf);
        else
            cycif_flow_plot(state.T, xc, 'Y',yc, ...
                'XTrans',xtDD.Value, 'YTrans',ytDD.Value, 'Sub',sub, 'Parent',gf);
        end
        ax = gca;                           % cycif_flow_plot created/cleared axes
        title(ax, sprintf('%s  (zoom via axes toolbar, then click Draw shape)', wname), ...
              'Interpreter','none');
        % overlay existing shape if we have one (e.g., editing an existing gate)
        if ~isempty(pendingCoords)
            G = previewGateStruct(xc,yc,pn);
            overlayShape(ax, G);
        end
        figure(gf);
    end

    % -- DRAW: right-click driven ROI on the current gate window ---------
    function onDrawShape(~,~)
        wkey = currentGateKey();
        if ~isKey(state.gateWindows,wkey) || ~isvalid(state.gateWindows(wkey))
            onPlot([],[]);   % implicit Plot if user skipped it
        end
        gf = state.gateWindows(wkey);
        figure(gf);
        ax = findobj(gf,'Type','axes');
        if isempty(ax); uialert(fig,'No axes in the plot window.','?'); return; end
        ax = ax(1);
        delete(findall(ax,'Tag','cycif_gate_overlay'));
        delete(findall(ax,'Tag','cycif_draw_preview'));

        coords = local_custom_draw(ax, shapeDD.Value);
        if isempty(coords)
            title(ax, sprintf('%s  (draw canceled)', nonempty(nameBox.Value,'(draft)')), ...
                  'Interpreter','none');
            return
        end
        pendingCoords = coords;

        G = previewGateStruct(readX(), readY(), readParent());
        overlayShape(ax, G);
        title(ax, sprintf('%s  (shape captured — click Save gate in manager)', ...
                          nonempty(nameBox.Value,'(draft)')), 'Interpreter','none');
        figure(fig);
    end

    function onSaveGate(~,~)
        nm = strtrim(nameBox.Value);
        if isempty(nm); uialert(fig,'Gate name required.','Missing name'); return; end
        if isempty(pendingCoords)
            uialert(fig,'No shape yet — click Plot, then Draw shape.','No shape'); return
        end
        [xc,yc,pn] = readEditorChans();

        existing = find(strcmp({state.gates.name}, nm));
        if ~isempty(existing) && (isempty(editingIdx) || existing ~= editingIdx)
            uialert(fig,'Name already exists.','Duplicate'); return
        end

        G = struct('name',nm,'shape',shapeDD.Value, ...
                   'xchan',xc,'ychan',yc, ...
                   'xtrans',xtDD.Value,'ytrans',ytDD.Value, ...
                   'coords',pendingCoords,'parent',pn, ...
                   'created',datetime('now'));

        if isempty(editingIdx)
            if isempty(state.gates); state.gates = G;
            else;                    state.gates(end+1) = G; end
            editingIdx = numel(state.gates);
        else
            oldName = state.gates(editingIdx).name;
            G.created = state.gates(editingIdx).created;
            state.gates(editingIdx) = G;
            if ~strcmp(oldName, nm)
                for i = 1:numel(state.gates)
                    if strcmp(state.gates(i).parent, oldName)
                        state.gates(i).parent = nm;
                    end
                end
                renameGateWindow(oldName, nm);
            end
        end
        % if this was a draft window, hand its key over to the saved gate
        if isKey(state.gateWindows, state.draftKey)
            renameGateWindow(state.draftKey, nm);
        end

        state.dirty = true;
        refreshTree(); refreshParentDD(); updateStats();
        refreshGateWindow(nm);
    end

    % -- Channel / parent read helpers ------------------------------------
    function [xc,yc,pn] = readEditorChans(); xc=readX(); yc=readY(); pn=readParent(); end
    function xc = readX();      xc = xDD.Value; end
    function yc = readY();      yc = yDD.Value; if strcmp(yc,'(none)'); yc=''; end; end
    function pn = readParent(); pn = parentDD.Value; if strcmp(pn,'(root)'); pn=''; end; end

    function G = previewGateStruct(xc, yc, pn)
        G = struct('name',nameBox.Value,'shape',shapeDD.Value, ...
                   'xchan',xc,'ychan',yc, ...
                   'xtrans',xtDD.Value,'ytrans',ytDD.Value, ...
                   'coords',pendingCoords,'parent',pn, ...
                   'created',datetime('now'));
    end

    function m = parentMask(pname)
        m = [];
        if isempty(pname); return; end
        try
            m = cycif_flow_apply(state.T, state.gates, 'GateName', pname);
        catch
            m = [];
        end
    end

    function updateStats()
        if isempty(editingIdx) || isempty(state.gates)
            statsLbl.Text = '—'; return
        end
        try
            G  = state.gates(editingIdx);
            m  = cycif_flow_apply(state.T, state.gates, 'GateName', G.name);
            N  = nnz(m);
            if ~isempty(G.parent)
                mp   = cycif_flow_apply(state.T, state.gates, 'GateName', G.parent);
                Npar = nnz(mp);
            else
                Npar = height(state.T);
            end
            statsLbl.Text = sprintf('N = %d   |   %% parent = %.1f   |   %% total = %.2f', ...
                N, 100*N/max(Npar,1), 100*N/height(state.T));
        catch ME
            statsLbl.Text = ['error: ' ME.message];
        end
    end

    % -- Gate window registry --------------------------------------------
    function key = currentGateKey()
        if isempty(editingIdx); key = state.draftKey;
        else;                   key = state.gates(editingIdx).name;
        end
    end

    function [gf, ax] = ensureGateWindow(key, name)
        if isKey(state.gateWindows,key) && isvalid(state.gateWindows(key))
            gf = state.gateWindows(key);
            set(gf,'Name',sprintf('Gate: %s', name));
        else
            state.windowCounter = state.windowCounter + 1;
            off = mod(state.windowCounter-1,8) * 28;
            gf  = figure('Name',sprintf('Gate: %s', name), ...
                         'Position',[260+off, 420-off, 640, 520], ...
                         'NumberTitle','off');
            state.gateWindows(key) = gf;
            gf.CloseRequestFcn = @(src,~) onGateFigClose(src, key);
        end
        ax = [];   % caller will fetch gca after cycif_flow_plot renders
    end

    function renameGateWindow(oldKey, newKey)
        if ~isKey(state.gateWindows, oldKey); return; end
        w = state.gateWindows(oldKey);
        remove(state.gateWindows, oldKey);
        if isvalid(w)
            state.gateWindows(newKey) = w;
            set(w, 'Name', sprintf('Gate: %s', newKey));
            w.CloseRequestFcn = @(src,~) onGateFigClose(src, newKey);
        end
    end

    function closeGateWindow(key)
        if isKey(state.gateWindows,key)
            w = state.gateWindows(key);
            if isvalid(w); delete(w); end
            remove(state.gateWindows, key);
        end
    end

    function closeAllGateWindows()
        ks = state.gateWindows.keys;
        for i = 1:numel(ks)
            w = state.gateWindows(ks{i});
            if isvalid(w); delete(w); end
        end
        state.gateWindows = containers.Map('KeyType','char','ValueType','any');
    end

    function onGateFigClose(src, key)
        if isKey(state.gateWindows, key); remove(state.gateWindows, key); end
        if isvalid(src); delete(src); end
    end

    % Re-plot a saved gate's window with its final shape overlay
    function refreshGateWindow(name)
        if ~isKey(state.gateWindows, name); return; end
        w = state.gateWindows(name);
        if ~isvalid(w); remove(state.gateWindows, name); return; end
        idx = find(strcmp({state.gates.name}, name), 1);
        if isempty(idx); return; end
        G   = state.gates(idx);
        sub = parentMask(G.parent);
        clf(w); delete(w.Children);
        if isempty(G.ychan)
            cycif_flow_plot(state.T, G.xchan, 'XTrans',G.xtrans, 'Sub',sub, 'Parent',w);
        else
            cycif_flow_plot(state.T, G.xchan, 'Y',G.ychan, ...
                'XTrans',G.xtrans, 'YTrans',G.ytrans, 'Sub',sub, 'Parent',w);
        end
        ax = gca;
        overlayShape(ax, G);
        % overlay immediate children in a muted color for context
        kids = find(strcmp({state.gates.parent}, name));
        for k = kids(:).'
            if strcmp(state.gates(k).xchan, G.xchan) && ...
               strcmp(nonempty(state.gates(k).ychan,''), nonempty(G.ychan,''))
                overlayShape(ax, state.gates(k), [0 0.45 0.9]);
            end
        end
        parentStr = tern(isempty(G.parent),'(root)', G.parent);
        title(ax, sprintf('%s    ← %s', G.name, parentStr), 'Interpreter','none');
    end

    % -- File actions ------------------------------------------------------
    function onLoadFile(~,~)
        start = tern(isempty(state.gateFile), pwd, state.gateFile);
        [f,p] = uigetfile({'*.mat','MAT files (*.mat)'},'Load gates', start);
        if isequal(f,0); return; end
        L = load(fullfile(p,f),'gates');
        closeAllGateWindows();
        state.gates    = L.gates;
        state.gateFile = fullfile(p,f);
        state.dirty    = false;
        editingIdx = []; pendingCoords = [];
        refreshTree(); refreshParentDD(); clearEditor();
    end

    function onSaveFile(~,~)
        if isempty(state.gateFile)
            [f,p] = uiputfile({'*.mat','MAT files (*.mat)'},'Save gates','gates.mat');
            if isequal(f,0); return; end
            state.gateFile = fullfile(p,f);
        end
        gates = state.gates;                                                %#ok<NASGU>
        save(state.gateFile,'gates');
        state.dirty = false;
        uialert(fig, sprintf('Saved %d gate(s) to %s', ...
                numel(state.gates), state.gateFile), 'Saved','Icon','success');
    end

    function onClose(~,~)
        if state.dirty
            sel = uiconfirm(fig, ...
                'You have unsaved changes. Save to file before closing?', ...
                'Unsaved changes', ...
                'Options',{'Save & close','Discard & close','Cancel'}, ...
                'DefaultOption',1, 'CancelOption',3);
            switch sel
                case 'Cancel';       return
                case 'Save & close'; onSaveFile();
            end
        end
        uiresume(fig);
    end
end

% =========================================================================
% LOCAL HELPERS
% =========================================================================
function g = local_empty_gate_array()
g = struct('name',{},'shape',{},'xchan',{},'ychan',{}, ...
           'xtrans',{},'ytrans',{},'coords',{}, ...
           'parent',{},'created',{});
end

function overlayShape(ax, G, col)
% Persistent shape overlay on an axes. Tagged so we can find+remove later.
if nargin < 3; col = [0.85 0.1 0.1]; end
hold(ax,'on');
switch lower(G.shape)
    case 'rect'
        c = G.coords;
        plot(ax, [c(1,1) c(2,1) c(2,1) c(1,1) c(1,1)], ...
                 [c(1,2) c(1,2) c(2,2) c(2,2) c(1,2)], ...
                 '-', 'Color',col, 'LineWidth',1.5, 'Tag','cycif_gate_overlay');
    case 'poly'
        V = G.coords;
        plot(ax, [V(:,1); V(1,1)], [V(:,2); V(1,2)], ...
                 '-', 'Color',col, 'LineWidth',1.5, 'Tag','cycif_gate_overlay');
    case 'thresh'
        yl = ylim(ax);
        plot(ax, [G.coords(1) G.coords(1)], yl, ...
                 '--', 'Color',col, 'LineWidth',1.5, 'Tag','cycif_gate_overlay');
end
end

function v = tern(c, a, b); if c; v = a; else; v = b; end; end
function s = nonempty(s, fallback); if isempty(char(s)); s = fallback; end; end

% =========================================================================
% CUSTOM RIGHT-CLICK ROI DRAWER
% =========================================================================
function coords = local_custom_draw(ax, shape)
%LOCAL_CUSTOM_DRAW  Right-click driven ROI drawing on an axes.
%
%   coords = local_custom_draw(AX, SHAPE)
%   SHAPE in {'poly','rect','thresh'}
%
% Interaction
%   right-click        — place a vertex
%   Enter / double-click — finish (if enough vertices)
%   Backspace          — remove the last vertex
%   Escape             — cancel (returns [])
%
% Output is in the axes' data space (same space cycif_flow stores).

fig = ancestor(ax,'figure');

% disable any interaction modes that would eat mouse clicks
try; zoom(fig,'off'); end        %#ok<TRYNC>
try; pan(fig,'off');  end        %#ok<TRYNC>
try; rotate3d(fig,'off'); end    %#ok<TRYNC>
try; datacursormode(fig,'off'); end %#ok<TRYNC>

% save state to restore when we're done
oldWBDF    = fig.WindowButtonDownFcn;
oldWKPF    = fig.WindowKeyPressFcn;
oldPointer = fig.Pointer;
oldAxMenu  = ax.UIContextMenu;

fig.Pointer      = 'crosshair';
ax.UIContextMenu = [];                          % suppress right-click menu

verts    = zeros(0,2);
preview  = gobjects(0);
canceled = false;

fig.WindowButtonDownFcn = @onClick;
fig.WindowKeyPressFcn   = @onKey;

showTitleHint();
uiwait(fig);

% restore
fig.WindowButtonDownFcn = oldWBDF;
fig.WindowKeyPressFcn   = oldWKPF;
fig.Pointer             = oldPointer;
ax.UIContextMenu        = oldAxMenu;
delete(preview(arrayfun(@isvalid, preview)));

if canceled || isempty(verts); coords = []; return; end
switch lower(shape)
    case 'poly';   coords = verts;
    case 'rect'
        c = verts(1:2,:);
        coords = [min(c(:,1)) min(c(:,2)); max(c(:,1)) max(c(:,2))];
    case 'thresh'; coords = verts(1,1);
end

% ---- nested callbacks ------------------------------------------------
    function onClick(~,~)
        cp = ax.CurrentPoint(1,1:2);
        if ~local_in_axes(ax, cp); return; end
        switch fig.SelectionType
            case 'alt'                     % right-click OR Ctrl+click
                verts(end+1,:) = cp;       %#ok<AGROW>
                redraw();
                autoFinish();
            case 'open'                    % double-click
                finishIfValid();
        end
    end

    function onKey(~, evt)
        switch evt.Key
            case 'return';    finishIfValid();
            case 'backspace'
                if ~isempty(verts); verts(end,:) = []; redraw(); end
            case 'escape'
                canceled = true; uiresume(fig);
        end
    end

    function autoFinish()
        % rect and thresh know exactly when they're complete
        switch lower(shape)
            case 'rect';   if size(verts,1) >= 2; uiresume(fig); end
            case 'thresh'; if size(verts,1) >= 1; uiresume(fig); end
        end
    end

    function finishIfValid()
        minV = struct('poly',3,'rect',2,'thresh',1);
        if size(verts,1) >= minV.(lower(shape)); uiresume(fig); end
    end

    function redraw()
        delete(preview(arrayfun(@isvalid, preview)));
        preview = gobjects(0);
        showTitleHint();
        if isempty(verts); return; end
        hold(ax,'on');
        preview(end+1) = plot(ax, verts(:,1), verts(:,2), 'ro', ...
            'MarkerFaceColor','r','MarkerSize',6, ...
            'Tag','cycif_draw_preview');                          %#ok<AGROW>
        switch lower(shape)
            case 'poly'
                if size(verts,1) >= 2
                    preview(end+1) = plot(ax, ...
                        [verts(:,1); verts(1,1)], [verts(:,2); verts(1,2)], ...
                        'r-','LineWidth',1.5,'Tag','cycif_draw_preview'); %#ok<AGROW>
                end
            case 'rect'
                if size(verts,1) >= 2
                    c = verts(1:2,:);
                    preview(end+1) = plot(ax, ...
                        [c(1,1) c(2,1) c(2,1) c(1,1) c(1,1)], ...
                        [c(1,2) c(1,2) c(2,2) c(2,2) c(1,2)], ...
                        'r-','LineWidth',1.5,'Tag','cycif_draw_preview'); %#ok<AGROW>
                end
            case 'thresh'
                yl = ylim(ax);
                preview(end+1) = plot(ax, [verts(1,1) verts(1,1)], yl, ...
                    'r--','LineWidth',1.5,'Tag','cycif_draw_preview');    %#ok<AGROW>
        end
    end

    function showTitleHint()
        switch lower(shape)
            case 'poly'
                msg = sprintf(['Polygon: right-click vertices (%d placed)   |   ' ...
                               'Enter = finish   Backspace = undo   Esc = cancel'], ...
                               size(verts,1));
            case 'rect'
                msg = sprintf('Rectangle: right-click 2 opposite corners (%d of 2)   |   Esc = cancel', ...
                              size(verts,1));
            case 'thresh'
                msg = 'Threshold: right-click once to place   |   Esc = cancel';
        end
        title(ax, msg, 'Interpreter','none');
    end
end

function tf = local_in_axes(ax, cp)
xl = xlim(ax); yl = ylim(ax);
tf = cp(1) >= xl(1) && cp(1) <= xl(2) && cp(2) >= yl(1) && cp(2) <= yl(2);
end