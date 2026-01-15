% plot_by_ND('treatment', dataloc, varargin)
% This code plots imaging experiment data by celltype or by treatment type.
%
% INPUTS
%   plotbyND - 'celltype', 'treatment', or 'pulseplot'. If set to 'celltype' function
%            will create a subplot for each cell type entered in data sheet
%            and plot all treatments on each subplot. 'treatment' creates a
%            subplot for each treatment.
%
%   dataloc - Nick's data structure
%        which is a stucture that includes.....
%        * ARE REQUIRED!! all else are optional (but some plots will not work)
%       *dataloc.d = the d cell array that ct_dataproc spits out
%       *dataloc.platemapd.pmd = the pmd output from putting your _platemap into iman_platemapreader
%       *dataloc.fold = structure with the following in them (only.fig is required)
%                .im
%                .proc
%       *            .fig = PATH TO WHAT FOLDER YOU WANT THE FIGURES TO BE SAVED TO
%                .if
%                .platemap
%       *dataloc.file = structure with the following in them (only .base is required)
%                .im
%                .celltracer
%       *            .base = the name of the folder where your imaging experiment is housed (or whatever title you want on the plot)
%                .global
%                .proc
%                .platemap
%       dataloc.z = the z cell array that ct_pulseanalysis spits out
%       dataloc.movieinfo = structure with the following movie metadata in it
%                .PixSizeX = converstion from pixels to um?
%                .PixSizeY = converstion from pixels to um?
%                .PixNumX = how many pixels wide is a movie? 
%                .PixNumY = how many pixels tall is a movie? 
%                .tsamp = how many minutes for each loop (e.g. 6)
%       dataloc.IFd = the d cell array that ct_dataproc spits out... if you ran it on IF data (that has been aligned already)
%       dataloc.normalizetype.(channel) = a string of how data was normalized (if it was) for each channel
%   
% OPTIONAL INPUTS
%   
%   channel - name of the channel(s) you wish to plot as a cell with
%   strings/chars. Default = {'ERKTR'}
%   If a cell array is given, it will do a plot for every channel in the
%   cell array. 
%
%   plottype - for all plotting, what type of plot (example: mean)
%   Options: mean, meanslope, meansemble, heatmap, stacks, spatialheatmap
%
%   for plotbyND = pulse analysis or arcos data
%   'analysischan' - 'frequency', 'dur', etc
%
%   looptime - number of minutes per loop (default = 6), for x (hour / time) scaling
%
%   numcells - number of cells to plot (default = all cells)
%
%   aftertreatment - which treatment do you want to consider when plotting (default = 1 aka the first treatment after the movie starts)
%   (NOTE: If you accidentally pick a tx higher than the max txs you have in that xy, it will default to the last tx for that xy)
%   
%   tstartaftertx - how many hours (IN HOURS) AFTER treatment (put in as aftertreatment) to START to consider for responses
%   tmaxaftertx - number of hours (IN HOURS) to consider after the treatment you specify (default = full length)
%
%   combinexys - combine xys that got the same treatments/celltype (default = true)
%
%   subset - only plot a subset of treatments or celltypes. For example, if
%   you want to plot only wells that were treated with EGF, subset = 'EGF'.
%   specsubset - only keep subset that matches ALL subset pars (default =
%   false)
%   
%   exclude - plot all treatments or celltypes excluding those specified.
%   For example, if you want to exclude treatments that got EGF, exclude =
%   'EGF'
%
%   tx_order - order the legened by a certain treatment. Defualt is Tx1.
%
%   font_size - size of the font for the titles of the figures
%
%   ymn = []; ymx = []; %if you want to fix the axis for the plots   
%
%   closefigs - close the figures the plot makes? (default = true)
%
%   ONLY CAN BE USED WITH plottype = 'stacks' & 'spatialstacks'
%   ntracks - number of tracks per stack (default = 5)
%   nstacks - number of stacks per condition (default = 1)
%   plotcolor - cell array of the plot colors (default = {'red', 'cyan'})
%   ntracks = 5; nstacks = 1;  normmax = 98; normmin = 2;
%   xys = [];
%
% looptime = 6; %in mins
% divtime = []
% p.specsubset = false;
% p.aftertreatment = 1; 
% p.ntracks = 5; p.nstacks = 1; p.plotcolor = {};
% p.plottype = {};
% p.analysischan = {'freq','durs'}; p.thisbythat = 'freq by dur';
% p.plotfromzero = false; p.withoutzeros = false;
% p.ncells = [];
% p.tstartaftertx = []; %is how many hours after treatment to START to consider for responses
% p.tmaxaftertx = []; %in hours,  max time for cells to be condidered after treatment
% p.zerohrtx = []; % which treatment should be the "zero" hour on the plot (default sets to same as p.aftertreatment)
% p.overlaptx = []; %tx where this tx is the second grouping for data (will be overlapped with the others)
% p.xpulse = 1;
% p.pulsethresh = 0.5; % min pulse frequency (per hour) to be considered pulsatile 
% p.imagesize = [1280, 1080]; p.nbins = [20,20];
% p.closefigs = true;
% p.overlapmeans = true;
% p.responsethresh = 0.2; % response threshold for pulse plots
% p.lengththresh = 1; % threshold for long versus short pulses (pulse plot)
% p.deltaprctl = [2, 98]; % for finding deltas in data to sort by
% p.addthreshline = true;
% p.ifchan = 'CY5_Nuc'; % IF channel to use for ifplot 
% p.iftrackwindow = 3; % number of tps to consider before end of plot for meaning
% p.showrawdata = true; %show raw data on pulse analysis and arcos box plots?
% p.errortype = 'sem'; %what type of stats for error bars in pulse and arcos analysis
% p.addstats = true; % do stats on the data it will ask which is the control condition (only for pulse and arcos analysis)
% p.addmedian = true; % show the median data on box plots (only for pulse and arcos analysis)
% p.squaresize = 150 % how big (in um) each square should be when subdividing an image
% p.nancolor = []; what color should nans be in heatmaps? (set to true for black, also if nonexistant color is given. Leave blank to leave to matlab's default)
% EXAMPLE
% plot_by('treatment',dataloc,'channel',{'ERKTR','AMPKAR'},'subset',{'EGF'},'plottype',{'mean','iBins'})
%
% -------------------------------------------------------------------------
% VERSION 4.0

function plot_by_MHY(plotby,dataloc,varargin)

%% Check paths
Path_Adder()

%% Input parsing
ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;

% Helper functions
isRGB = @(x) isnumeric(x) && numel(x) == 3 && all(x >= 0) && all(x <= 1);
isStringOrCellstr = @(x) ischar(x) || isstring(x) || iscellstr(x);
isNonNegScalar = @(x) isnumeric(x) && isscalar(x) && x >= 0;
isCellStr = @(x) iscell(x) && all(cellfun(@ischar, x));

% Add parameters (alphabetical)
addParameter(ip, 'addif', false, @islogical);
addParameter(ip, 'addmedian', true, @islogical);
addParameter(ip, 'addstats', true, @islogical);
addParameter(ip, 'addthreshline', true, @islogical);
addParameter(ip, 'aftertreatment', 1, isNonNegScalar);
addParameter(ip, 'analysischan', {'freq','durs'}, isCellStr);
addParameter(ip, 'channel', {'ERKTR'});
addParameter(ip, 'closefigs', true, @islogical);
addParameter(ip, 'combinexys', true, @islogical);
addParameter(ip, 'control', [], @(x) isempty(x) || ischar(x) || iscell(x));
addParameter(ip, 'deltaprctl', [2, 98], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'divtime', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'errortype', 'sem', @(x) any(validatestring(x, {'sem','std','ci'})));
addParameter(ip, 'exclude', [], @(x) isnumeric(x) || iscell(x));
addParameter(ip, 'facetby', {}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'font_size', 8, isNonNegScalar);
addParameter(ip, 'groupby', {}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'ifchan', 'CY5_Nuc', @ischar);
addParameter(ip, 'iftrackwindow', 3, isNonNegScalar);
addParameter(ip, 'imagesize', [1280, 1080], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'lengththresh', 1, isNonNegScalar);
addParameter(ip, 'looptime', 6, isNonNegScalar);
addParameter(ip, 'nancolor', [], @(x) isempty(x) || isRGB(x) || (islogical(x) && x));
addParameter(ip, 'nbins', [20, 20], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'ncells', [], @(x) isempty(x) || isscalar(x));
addParameter(ip, 'nogene', 0, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'nstacks', 1, isNonNegScalar);
addParameter(ip, 'ntracks', 5, isNonNegScalar);
addParameter(ip, 'overlapmeans', true, @islogical);
addParameter(ip, 'overlaptx', [], @(x) isempty(x) || ischar(x) || iscell(x));
addParameter(ip, 'plotcolor', {}, @(x) iscell(x) || isstring(x));
addParameter(ip, 'plotfromzero', false, @islogical);
addParameter(ip, 'plottype', {'mean'}, isCellStr);
addParameter(ip, 'printstyle', 'svg', @(x) any(validatestring(x, {'svg','png','pdf','eps'})));
addParameter(ip, 'pulsethresh', 0.5, isNonNegScalar);
addParameter(ip, 'qtiles', true, @islogical);
addParameter(ip, 'responsethresh', 0.2, isNonNegScalar);
addParameter(ip, 'showrawdata', true, @islogical);
addParameter(ip, 'smooth', [], @(x) isempty(x) || isvector(x));
addParameter(ip, 'specsubset', false, @islogical);
addParameter(ip, 'squaresize', 50, isNonNegScalar);
addParameter(ip, 'standardizeplots', true, @islogical);
addParameter(ip, 'subset', [], @(x) isnumeric(x) || iscell(x));
addParameter(ip, 'tbeforetx', [], @(x) isempty(x) || isscalar(x));
addParameter(ip, 'thisbythat', 'freq by dur', @ischar);
addParameter(ip, 'tmaxaftertx', [], @(x) isempty(x) || isscalar(x));
addParameter(ip, 'tmaxback', 3, isNonNegScalar);
addParameter(ip, 'tstartaftertx', [], @(x) isempty(x) || isscalar(x));
addParameter(ip, 'tx_order', 1, isNonNegScalar);
addParameter(ip, 'withoutzeros', false, @islogical);
addParameter(ip, 'xys', [], @(x) isempty(x) || isnumeric(x));
addParameter(ip, 'xpulse', 1, isNonNegScalar);
addParameter(ip, 'ymn', [], @(x) isnumeric(x) || iscell(x));
addParameter(ip, 'ymx', [], @(x) isnumeric(x) || iscell(x));
addParameter(ip, 'zerohrtx', [], @(x) isempty(x) || ischar(x));

% parse
parse(ip, varargin{:});
p = ip.Results;

% make sure p.ymn and p.ymx are cell arrays
if ~isempty(p.ymn) && ~iscell(p.ymn), p.ymn = {p.ymn}; end
if ~isempty(p.ymx) && ~iscell(p.ymx), p.ymx = {p.ymx}; end
if islogical(p.nancolor) && p.nancolor, p.nancolor = [0 0 0]; end

% If groupby/facetby are specified, override plotby behavior
use_custom_grouping = ~isempty(p.groupby) || ~isempty(p.facetby); 
if use_custom_grouping
    plotby = 'custom';  % special mode to trigger downstream logic
end

switch lower(plotby)
    case {'treatment','treatments','tx','txs'}
        plotby = 'treatment';
    case {'cell','celltype','celltypes'}
        plotby = 'celltype';
    case {'pulse','pulseanal','pulses','pulseplot'}
        plotby = 'pulseplot';
    case {'treatmentoverlay','overlay'}
        plotby = 'treatmentoverlay';
    case {'custom','group','grouped'}
        plotby = 'custom';
        % Don't raise an error, we are using groupby/facetby instead
    otherwise
        error('Accepted plotby types: treatment, celltype, pulseplot, or overlay. plotby must be a string.');
end

%% New dataloc parameters
if isfield(dataloc,'movieinfo')
    p.looptime = dataloc.movieinfo.tsamp;
    p.imagesize = [dataloc.movieinfo.PixNumY,dataloc.movieinfo.PixNumX];
    p.imagscaling = [dataloc.movieinfo.PixSizeY,dataloc.movieinfo.PixSizeX];
end 

p.tktm = 60/p.looptime; %set up how many tps per hour
if ~isempty(p.tstartaftertx);  p.tstartaftertx = p.tstartaftertx * p.tktm; end 
if ~isempty(p.tmaxaftertx); p.tmaxaftertx = p.tmaxaftertx * p.tktm; p.tmaxaftertx = floor(p.tmaxaftertx); end %convert respond time from hours to tps
if ~isempty(p.tbeforetx); p.tbeforetx = p.tbeforetx * p.tktm; p.tbeforetx = floor(p.tbeforetx); end
if ~isempty(p.divtime)
    p.divtime = p.divtime * p.tktm; % convert division time allowed into tps (from hours)
end

if ~iscell(p.channel); p.channel = {p.channel}; end %put the channel in a cell
if ~iscell(p.plottype); p.plottype = {p.plottype}; end %put the trackplot in a cell
if ~isempty(p.subset) && ~iscell(p.subset); p.subset = {p.subset}; end %put the subset in a cell
if ~isempty(p.exclude) && ~iscell(p.exclude); p.exclude = {p.exclude}; end %put the exclude in a cell
if ~iscell(p.analysischan); p.analysischan = {p.analysischan}; end


%% get cell type names
if isempty(dataloc.platemapd.pmd); warning('Your dataloc does not have platemap data, rerun fullproc and ensure it finds the platemap. Then try again.'); return; end
nct = ceil(size(dataloc.platemapd.pmd.Cell,3)/3);  ci = 3*((1:nct)-1) + 1;

if ~p.nogene % cat celltype to gene unless indicated
    genecat = cell(size(dataloc.platemapd.pmd.Gene(:,:,1)));
    for s = 1: size(dataloc.platemapd.pmd.Gene,3)
        genecat = cellfun(@(x,y)[x,y], genecat, dataloc.platemapd.pmd.Gene(:,:,s), ...
            'UniformOutput', false);
    end

%   Cat genes to cells
cnames = cellfun(@(x,y)[x,'_',y], dataloc.platemapd.pmd.Cell(:,:,ci), ...
            repmat(genecat,[1,1,numel(ci)]), ...
            'UniformOutput', false);
else
    cnames = cell(size(dataloc.platemapd.pmd.Cell(:,:,1)));
    for s = 1:numel(ci)
    cnames = cellfun(@(x,y)[x,y],cnames,dataloc.platemapd.pmd.Cell(:,:,s),'Un',0);
    end
    
    nancell = find(~cell2mat(cellfun(@ischar,cnames,'Un',0)));
    for s = 1:numel(nancell)
       cnames{nancell(s)} = ''; 
    end
end
%   Remove any non-word characters to make usable as names
cnames = regexprep(cnames, {'\W','^[\d_]*(\w)'}, {'','$1'});
%   Make a name for each unique catted string
[celltypes] = unique(cnames(cellfun(@ischar,cnames)));
%   Disregard invalid names and warn
gn = cellfun(@isvarname,celltypes);
if any(~gn)
    if any(~cellfun(@(x)strcmp(x,'_'),celltypes(~gn))); warning(['Invalid cell name found for idx: ', celltypes{~gn}]); end
    celltypes = celltypes(gn);
end

celltypes = celltypes(~cellfun(@isempty,celltypes)); % discard empty name fields
% get rid of @ density if present for fieldname
cellfn = cellfun(@(x)x{1},regexp(celltypes,'@','split'),'un',0);
% remove any spaces in name
cellfn = arrayfun(@(x)regexprep(cellfn{x},'[^a-zA-Z0-9]','_'),1:numel(cellfn),'un',0);
% add x to beginning of name starts with number
cellfn = arrayfun(@(x)regexprep(cellfn{x},'(^[\d_]+\w)','x$1'),1:numel(cellfn),'un',0);
% make list of xys for each celltype
cellfn=cellfn';
%   Assign matching xy positions to each compiled cell name
for s = 1:numel(celltypes)
    %   Assemble list of xys with the current cname
    idx.(celltypes{s}) = [dataloc.platemapd.pmd.xy{strcmp(celltypes{s},cnames)}];  
end


%% Find good xys that actually have data
goodxy = false(1,max([dataloc.platemapd.pmd.xy{:}])); 
for ii=1:size(dataloc.d,2)
    if ~isempty(dataloc.d{ii}) && isfield(dataloc.d{ii},'cellindex') && ~isempty(dataloc.d{ii}.cellindex)
        goodxy(ii)=1;
    end
end

%% get treatments and find unique combos
catTxnames = {'pTx','Tx1','Tx2','Tx3','Tx4'}; % possible treatment fields
[~,I]= intersect(catTxnames, fieldnames(dataloc.platemapd.pmd)); % tx fields in experiment
I = sort(I); catTxnames = catTxnames(I);

Txcat = cell(size(dataloc.platemapd.pmd.Cell,1),size(dataloc.platemapd.pmd.Cell,2)); % initialize cell array
linetp = cell(size(dataloc.platemapd.pmd.Cell,1),size(dataloc.platemapd.pmd.Cell,2)); % initialize cell array

for s = 1:numel(catTxnames) % cat each treatment in a loop
    numtreat = size(dataloc.platemapd.pmd.(catTxnames{s}),3)/5;
    tid = logical(repmat([1,1,1,0,0],1,numtreat));
    Txcat = cat(3,Txcat, cellfun(@(x)cat(2,x{:}), cat(3,num2cell(dataloc.platemapd.pmd.(catTxnames{s})(:,:,tid),3)),'Un',0));
    linetp = cat(3,linetp, cellfun(@(x)cat(2,x{:}),...
    cat(3,num2cell(dataloc.platemapd.pmd.(catTxnames{s})(:,:,4),3)),'Un',0));    
end
linetp = linetp(:,:,2:end); % remove empty layer from initialization
linetp = num2cell(linetp,3);
% color sequence for treatment lines
linecolor = parula(numel(catTxnames)+1);

for s=1:numel(Txcat)
   if strcmp(Txcat{s},'')
       Txcat{s} = NaN;
   end
end

stridx = cellfun(@ischar,Txcat);
spacer =  repmat({'+'},size(Txcat)); % add + between txs for readability
spacer(~stridx) = {[]};

% add + between treatment names
Txcat2 = Txcat(:,:,2);
for s = 2:size(Txcat,3)-1
    Txcat2 = cat(3,Txcat2,spacer(:,:,s+1),Txcat(:,:,s+1));
end
% cat Tx name and space together into one long string
Txcat2 = cellfun(@(x)cat(2,x{:}),cat(3,num2cell(Txcat2,3)),'Un',0);
stridx = cellfun(@ischar,Txcat2);
treatments = unique(Txcat2(cellfun(@ischar,Txcat2)));

% find the unique treatment combos
txs = unique(Txcat2(cellfun(@ischar,Txcat2)));
% replace spaces and / so txs can be used as fieldnames
for s = 1:numel(txs)
    txs{s} = regexprep(txs{s},char(0),'');
    txs{s} = regexprep(txs{s},'(^|\s|+)','');
    txs{s} = regexprep(txs{s},'[^a-zA-Z0-9]','_');
    
end
for s = 1:numel(Txcat2)
    if ischar(Txcat2{s})
    Txcat2{s} = regexprep(Txcat2{s},char(0),'');
    Txcat2{s} = regexprep(Txcat2{s},'(^|\s|+)','');
    Txcat2{s} = regexprep(Txcat2{s},'[^a-zA-Z0-9]','_');
    end
end

% restrict to subset
if ~iscell(p.subset) && ~isempty(p.subset); p.subset = {p.subset}; end
if ~isempty(p.subset)
    txI = []; ctI = [];
    for iSubset = 1:numel(p.subset)
        txI2 = contains(txs,p.subset{iSubset},'IgnoreCase',true);
        ctI2 = contains(cellfn,p.subset{iSubset},'IgnoreCase',true);
        txI = [txI, txI2]; ctI =[ctI, ctI2];
    end
    
    if p.specsubset
        txI = all(txI,2);
        ctI = all(ctI,2);
    else
        txI = any(txI,2);
        ctI = any(ctI,2); 
    end
    
    if any(txI)
       txs = txs(txI);
       treatments = treatments(txI);
    end
    if any(ctI) 
       cellfn = cellfn(ctI);
    end
end

if ~isempty(p.exclude)
    txI = []; ctI = [];
    for iSubset = 1:numel(p.exclude)
        txI2 = contains(txs,p.exclude{iSubset},'IgnoreCase',true);
        ctI2 = contains(cellfn,p.exclude{iSubset},'IgnoreCase',true);
        txI = [txI, txI2]; ctI =[ctI, ctI2];
    end  
    txI = any(txI,2);
    ctI = any(ctI,2); 
    if any(txI)
       txs = txs(~txI);
       treatments = treatments(~txI);
    end
    if any(ctI) 
       cellfn = cellfn(~ctI);
    end
end
% sort the treatments into descending order
splitTxs = cellfun(@(x)strsplit(x,'+'), treatments, 'UniformOutput', false);
if ~isempty(splitTxs)
    % check for inconsistent sizes
    sz = cellfun(@(x)numel(x),splitTxs);
    sz = max(sz);
    splitTxss = strings(numel(treatments),sz);
    for kk = 1:numel(treatments)
        for jj = 1:numel(splitTxs{kk})
            splitTxss(kk,jj) = splitTxs{kk}{jj};
        end
    end
    if numel(splitTxs)>2
        [~,I] = sort(splitTxss);
        txs = txs(I(:,p.tx_order)); treatments = treatments(I(:,p.tx_order));
    end
end

% beautify treatment names for labels
% remove extra spaces in name caused by NaNs to make pretty names
treatments = arrayfun(@(x)regexprep(treatments{x},char(0),''),1:size(treatments),'un',0)';
treatments = arrayfun(@(x)regexprep(treatments{x},'(^|\s)',''),1:size(treatments),'un',0)';
% remove double +s
treatments = arrayfun(@(x)regexprep(treatments{x},'++','+'),1:size(treatments),'un',0)';
% replace _ with ' ' 
treatments = arrayfun(@(x)regexprep(treatments{x},'_',' '),1:size(treatments),'un',0)';

% FIX determine if end point controls were used, if so, remove that from the
% titles and plotted data, also self normalized for the cell

% vector of xys for packing txs into idx structure.
xymat = dataloc.platemapd.pmd.xy(stridx);
% now pack concatenated txs into the idx format 
p.subsetNums = 0;
for s = 1:size(txs,1)
idx.(txs{s}) = cat(2,xymat{strcmp(txs{s},Txcat2(stridx))}); % find xy 
p.subsetNums = p.subsetNums + size(idx.(txs{s}),2);
end

linetp = linetp(stridx);

if p.combinexys
   [dataloc,idx,linetp,xymat,goodxy,p] = Combine_XYs(dataloc,txs,cellfn,idx,xymat,linetp,goodxy,p);
end

%% Plotting

switch plotby
    % Plot by Cell type
    case 'treatment'
        titlevec = cellfn;
        if numel(txs) > 2
            cmap = rainbow(numel(txs));
        else
            cmap = [0,0.4,0.7;0.6,0.1,0.2];
        end
        legname = treatments;
        for iChan = 1:numel(p.channel)
            channel = p.channel{iChan};
            main_plotting_func(dataloc,channel,cellfn,txs,idx,...
                xymat,cmap,linetp,titlevec,legname,goodxy,p);
        end %loop over channels

    case 'celltype'
        titlevec = treatments; % use beautified tx names
        if numel(cellfn) > 2
            cmap = rainbow(numel(cellfn)+1);
        else
            cmap = [0,0.4,0.7;0.6,0.1,0.2];
        end
        legname = cellfn;
        for iChan = 1:numel(p.channel)
            channel = p.channel{iChan};

            main_plotting_func(dataloc,channel,txs,cellfn,idx,...
                xymat,cmap,linetp,titlevec,legname,goodxy,p);

        end %loop over channels

    case 'custom'
        titlevec = treatments; % use beautified tx names
        if numel(cellfn) > 2
            cmap = rainbow(numel(cellfn)+1);
        else
            cmap = [0,0.4,0.7;0.6,0.1,0.2];
        end
        legname = cellfn;
        for iChan = 1:numel(p.channel)
            channel = p.channel{iChan};

            main_plotting_func(dataloc,channel,p.facetby,p.groupby,idx,...
                xymat,cmap,linetp,titlevec,legname,goodxy,p);

        end %loop over channels

    case {'pulseplot','pulse','pa','pulse analysis','pulseanalysis'}
        for iChan = 1:numel(p.channel)
            channel = p.channel{iChan};
            titlevec = treatments;
            Pulse_Plotter(dataloc,channel,cellfn,txs,idx,xymat,linetp,titlevec,goodxy,p)
        end
end
end
%% Main plotting function 
function [TXX,nrow,ncol] = main_plotting_func(dataloc,channel,facetby,groupby,idx,xymat,cmap,linetp,titlevec,legname,goodxy,p)
if ~iscell(channel) && ischar(channel); channel = {channel}; end
% remove any underscores in legend text to avoid subscript shenanigans
legname = arrayfun(@(x)regexprep(legname{x},'_',' '),1:numel(legname),'un',0)';
titlevec = arrayfun(@(x)regexprep(titlevec{x},'_',' '),1:numel(titlevec),'un',0)';
TXX=[]; nrow = []; ncol=[];


%% actual plotting part

for iPlot = 1:numel(p.plottype)
        for s = 1:numel(facetby)
            figgy = tiledlayout('flow','TileSpacing','compact','Padding','compact');
            legvec = []; legtxt = {}; sc = 0;
            for sa = 1:numel(groupby) % plot line for each unique treatment/cell
                xy = intersect(idx.(facetby{s}),idx.(groupby{sa})); % find xy(s) for treatment/celltype
                xys4sub = length(xy);
                firstxy = 0;
                for sb = 1:xys4sub % plot line for each xy
                    if xy(sb) > numel(dataloc.d) || isempty(dataloc.d{xy(sb)}) || any(~isfield(dataloc.d{xy(sb)}.data,channel)); else
                    if goodxy(xy(sb))% make sure there is data to plot.
                        if p.plottype{iPlot} == "mean" && p.overlapmeans
                            firstxy = firstxy + 1; % store tx times for each xy for plotting lines later
                            if firstxy == 1; sc = sc+1; AH = nexttile(figgy,sc); hold on; end
                        else %if youre not combining xys, plot them one at a time
                            AH = nexttile(figgy); hold on
                            firstxy = 1;
                            sc = sc+1; % index for saving individual line objects for input to legend.
                        end
                        
                        treatment4XY = str2double(linetp{cell2mat(cellfun(@(x)any(x == xy(sb)),xymat,'un',0))});
                        anyunique = unique(treatment4XY);
                        txs = anyunique(~isnan(anyunique)); txs = txs(txs>0); %get the treatment times not including pretreatments
                        if ~isempty(txs)
                            if numel(txs) < p.aftertreatment; thisTX=txs(end); else; thisTX=txs(p.aftertreatment); end % get treatment number or last treatment
                        end
                        
                        [numcells, tMax] = size(dataloc.d{xy(sb)}.data.(channel{1})); % how many cells and tps are there?
                        
                        if ~isempty(p.ncells) && p.ncells <= numcells; numcells = p.ncells; end % how many cells worth of data? and check you have that many cells
                        if isempty(p.tmaxaftertx); tracklength = tMax; else; tracklength = p.tmaxaftertx + thisTX; end % get the end point of the data plotted
                        if tracklength > tMax; tracklength = size(dataloc.d{xy(sb)}.data.(channel{1}),2); end % if the track is longer than the data, fix that

                        firsttp = 1;

                        if p.plotfromzero
                            if ~isempty(txs); firsttp = thisTX;  end % offset all data to the first treatment % thisTx had -1? idk
                        elseif ~isempty(p.tbeforetx) % how far before the tx do you want to plot? if its before time starts set it to tp 1.
                            firsttp = thisTX - p.tbeforetx;
                            if (firsttp > 0); if ~isempty(txs); txs=txs-firsttp; end  
                            else; firsttp = 1;
                            end
                        end
                        

                        if ~isempty(txs)
                            txs = txs(txs<tracklength);
                        end

                        if numcells < (p.ntracks * p.nstacks); ntracks = floor(numcells / p.nstacks); else; ntracks = p.ntracks; end %if doing stacks, make sure you have enough cells for it, if not fix it

                    PlotRules = PlotFormatRules(channel, p.plottype{iPlot});
                    NumChans = numel(channel);
                        switch p.plottype{iPlot}
                            case {'mean','means'}
                                 for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right;  AH=gca; end
                                    if iscell(p.ymn); ymn = p.ymn{iChan}; ymx = p.ymx{iChan};
                                    else; ymn = p.ymn; ymx= p.ymx; end
                                    pData = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    if ~isempty(p.smooth); pData = movmean(pData,p.smooth,2); end % smooth data if requested
                                    h{sc} = ct_trackvis_MHY(AH,'mean',pData,'cmap',cmapcol,'ymn',ymn,'ymx',ymx,'PLOTQTILES',p.qtiles);
                                    clear pData;
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                    p.plottype{iPlot} = 'mean';
                                 end

                            case {'mean vs mean'}
                                xx = nanmean(dataloc.d{xy(sb)}.data.(channel{1})(:,firsttp:tracklength),1);
                                yy = nanmean(dataloc.d{xy(sb)}.data.(channel{2})(:,firsttp:tracklength),1);

                                if ~isempty(txs)
                                    colorlst = {'autumn','cool','winter','copper','pink'}; % 6 colors to pick from
                                    allowedTxs = (firsttp <= txs) & (txs <= tracklength);
                                    tTxs = txs(allowedTxs);
                                    tTxs = [firsttp-1, tTxs', tracklength];
                                    cc = [];
                                    for iTx = 1:numel(tTxs)-1
                                        cc2 = eval([colorlst{iTx},'(',num2str((tTxs(iTx+1)-tTxs(iTx))),')']);
                                        cc = [cc;cc2];
                                    end
                                else; cc = parula((tracklength-firsttp)+1); % get colors per time
                                end
                                h{sc} = scatter(xx,yy,25,cc,'filled','MarkerEdgeColor','k','LineWidth',0.25,'MarkerFaceAlpha',0.75);
                                xlabel(['mean ',channel{1}])
                                ylabel(['mean ',channel{2}])
                                colormap(cc);
                                c = colorbar;
                                clim([firsttp-1,tracklength])
                                tickz = linspace(firsttp-1,tracklength,5);
                                if ~isempty(txs)
                                    tickz = [tickz, txs(allowedTxs)'];
                                    tickz = sort(tickz)';
                                    tickz2 = tickz/p.tktm;
                                    tickzName = string(num2str(tickz2,3));
                                    whereTxs = ismember(tickz,txs(allowedTxs)');
                                    tickzName(whereTxs) = 'Treatment';
                                else; tickz2 = tickz/p.tktm;
                                      tickzName = string(num2str(tickz2,3));
                                end
                                c.Ticks = tickz;
                                c.TickLabels = tickzName;
                                c.Label.String = 'Time (hrs)';
                                labPos = c.Label.Position; labPos(1) = -1;
                                c.Label.Position = labPos;
                                h1 = lsline;
                                set(h1(1),'color','black')
                                p.plottype{iPlot} = 'mean vs mean';
                            
                            case {'meanslope', 'mean slope'}
                                 for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right;  AH=gca; end
                                    if iscell(p.ymn); ymn = p.ymn{iChan}; ymx = p.ymx{iChan};
                                    else; ymn = p.ymn; ymx= p.ymx; end
    
                                    datnanmean = nanmean(dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength),1);
                                    h{sc} = ct_trackvis_MHY(AH,'mean',dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength),'cmap',cmapcol,'ymn',ymn,'ymx',ymx);
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                    
                                    % get the slope of the channel from the tpof the treatment until the next one and so on
                                    if ~isempty(txs)
                                        numTxs = numel(txs);
                                        % loop through every treatment and fit lines btwn them or the end
                                        for iTreat = 1:numTxs
                                            if iTreat < numTxs
                                                xs = txs(iTreat):txs(iTreat+1);
                                                y1 = datnanmean(:,txs(iTreat):txs(iTreat+1));
                                            else
                                                xs = txs(iTreat):tracklength;
                                                y1 = datnanmean(:,txs(iTreat):tracklength);
                                            end
                                            P = polyfit(xs,y1,1); %fit a slope to the data
                                            yfit = P(1)*xs+P(2);  % P(1) is the slope and P(2) is the intercept
                                            plot(xs,yfit,'k-.','LineWidth',2);
                                            text(floor(mean(xs,'all')),(max(y1)*1.1),num2str(P(1)));
                                        end 
                                    end % empty tx check
                                 end

                            
                            case {'mean hill', 'meanhill'}
                                if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right; AH=gca; end
                                datnanmean = nanmean(dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength),1);
                                h{sc} = plot(datnanmean,'Color',cmapcol,'LineStyle','-','LineWidth',1.1);
                                set(gca,'YColor','k');
                                ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                
                                % get the slope of the channel from the tpof the treatment until the next one and so on
                                if isempty(txs); else
                                    numTxs = numel(txs);
                                    % loop through every treatment and fit lines btwn them or the end
                                    for iTreat = 1:numTxs
                                        if iTreat < numTxs
                                            xs = txs(iTreat):txs(iTreat+1);
                                            y1 = datnanmean(:,txs(iTreat):txs(iTreat+1));
                                        else
                                            xs = txs(iTreat):tracklength;
                                            y1 = datnanmean(:,txs(iTreat):tracklength);
                                        end
                                        % second cell array of HillOutput is [max, slope, half max (,and intercept)]
                                        HillOutput = HillFunctionFit(1:size(xs,2),y1,3,max(y1),((max(y1)-min(y1))/size(y1,2)),max(y1)/2,[]);
                                        plot(xs,HillOutput{1}(2,:),'k-.','LineWidth',1);
                                        %text(median(xs,'all'),(max(y1)),num2str(HillOutput{2}(1)));

                                    end 
                                end % empty tx check
                                p.plottype{iPlot} = 'mean hill';

                            case 'meansemble'
                                 for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    pData = dataloc.d{xy(sb)}.data.(channel{iChan})(1:numcells,firsttp:tracklength);
                                    if ~isempty(p.smooth); pData = movmean(pData,p.smooth,2); end % smooth data if requested
                                    h{sc} = ct_trackvis(AH, 'meansemble',pData,'cmap',cmapcol,'morelabels',false);
                                    xlim([0, tracklength]);
                                 end
                            case 'meanoverlaytreat'
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    h{sc} = plot(nanmean(dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength)),...
                                        'Color',cmapcol,'LineStyle','-','LineWidth',1.1);
                                    xlim([0, tracklength-firsttp]); 
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                end

                            case {'stacks','stack'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if NumChans > 1; PlotQs = false; else; PlotQs = true; end
                                    if ~isempty(p.ymn) && ~isempty(p.ymx); ymn = p.ymn(iChan); ymx = p.ymx(iChan); else; ymx = []; ymn=[]; end
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right; AH=gca; end
                                    pData = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    if ~isempty(p.smooth); pData = movmean(pData,p.smooth,2); end % smooth data if requested
                                    h{sc} = ct_trackvis_MHY(AH, 'stack', pData,'ntracks',ntracks,'nstacks',p.nstacks,'cmap',cmapcol,'PLOTQTILES',PlotQs,'morelabels',false,'ymn',ymn,'ymx',ymx,'treatment',txs); 
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                end
                                p.plottype{iPlot} = 'stacks';
                                
                            case {'stacks sorted','sorted stacks'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if NumChans > 1; PlotQs = false; else; PlotQs = true; end
                                    if ~isempty(p.ymn) && ~isempty(p.ymx); ymn = p.ymn{iChan}; ymx = p.ymx{iChan}; else; ymx = []; ymn=[]; end
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right; AH=gca; end
                                    pDat = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    if ~isempty(p.smooth); pDat = movmean(pDat,p.smooth,2); end % smooth data if requested
                                    [~,I]=sort(sum(isnan(pDat),2));
                                    h{sc} = ct_trackvis_MHY(AH, 'stack',pDat(I,:),'ntracks',ntracks,'nstacks',p.nstacks,'cmap',cmapcol,'PLOTQTILES',PlotQs,'morelabels',false,'ymn',ymn,'ymx',ymx,'treatment',txs); 
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                end
                                p.plottype{iPlot} = 'sorted stacks';
                                
                            case {'stacks sorted v2','sorted stacks v2'}
                                [~,I]=sort(sum(isnan(dataloc.d{xy(sb)}.data.(channel{1})(:,firsttp:tracklength)),2));
                                for iCell = 1:ntracks*p.nstacks
                                    hold on;
                                    for iChan = 1:NumChans
                                        if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right; end
                                        pData = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                        if ~isempty(p.smooth); pData = movmean(pData,p.smooth,2); end % smooth data if requested
                                        plot(pData(I(iCell),:),'Color',PlotRules(iChan).CmapColor,'LineStyle','-');
                                        set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                    end
                                    title([legname{sa}, ' Cell ',num2str(iCell)])
                                    if ~isempty(txs)
                                        xline(txs,'--','LineWidth', 1)
                                    end
                                    findZero(txs,thisTX,tracklength,firsttp,p)
                                    xlim([firsttp-1,tracklength])
                                    hold off;
                                    if ~(iCell == ntracks*p.nstacks); nexttile; end
                                end
                                p.plottype{iPlot} = 'sorted stacks v2'; firstxy = false;

                                    

                            case {'spatialstacks','spatialstack','spatial stacks'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if NumChans > 1; PlotQs = false; else; PlotQs = true; end
                                    if ~isempty(p.ymn) && ~isempty(p.ymx); ymn = p.ymn(iChan); ymx = p.ymx(iChan); else; ymx = []; ymn=[]; end
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right; AH=gca; end
                                    if iChan ==1
                                        averagex = nanmean(dataloc.d{xy(sb)}.data.XCoord,2);
                                        averagey = nanmean(dataloc.d{xy(sb)}.data.YCoord,2);
                                        coords = [averagex, averagey];
                                        DIST = pdist2(coords,coords);
                                        LINK = linkage(DIST,'average');
                                        leafOrder = optimalleaforder(LINK,DIST);
                                    end
                                    
                                    h{sc} = ct_trackvis_MHY(AH, 'stack', dataloc.d{xy(sb)}.data.(channel{iChan})(leafOrder',firsttp:tracklength),'ntracks',ntracks,'nstacks',p.nstacks,'cmap',cmapcol,'PLOTQTILES',PlotQs,'morelabels',false,'ymn',ymn,'ymx',ymx,'treatment',txs); 
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                end
                                p.plottype{iPlot} = 'spatial stacks';
                                %fix quantiles being weird in ct_trackvis

                            case {'selected stacks'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    tp2Check = input('Check which tp? ',"s");
                                    tp2Check = str2double(tp2Check);
    
                                    x = dataloc.d{xy(sb)}.data.XCoord(:,tp2Check);
                                    y = dataloc.d{xy(sb)}.data.YCoord(:,tp2Check);
                                    delMe = figure; scatter(x,y,'filled');
                                    a = (1:size(x,1))'; b = num2str(a); c = cellstr(b);
                                    text(x, y, c);
                                    pause(1)
                                    cellz = input('Pick your cells: ',"s");
                                    cellz = str2double(cellz)';
                                    close(delMe)
    
                                    if NumChans > 1; PlotQs = false; else; PlotQs = true; end
                                    if ~isempty(p.ymn) && ~isempty(p.ymx); ymn = p.ymn(iChan); ymx = p.ymx(iChan); else; ymx = []; ymn=[]; end
                                    if iChan == 1 && NumChans > 1; yyaxis left; elseif iChan > 1; yyaxis right;  AH=gca; end
                                    pDat = dataloc.d{xy(sb)}.data.(channel{iChan})(cellz,firsttp:tracklength);
                                    h{sc} = ct_trackvis_MHY(AH, 'stack',pDat,'ntracks',ntracks,'nstacks',p.nstacks,'cmap',cmapcol,'PLOTQTILES',PlotQs,'morelabels',false,'ymn',ymn,'ymx',ymx,'treatment',txs); 
                                    set(gca,'YColor','k');
                                    ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                end 
                                firstxy = false;

                            case {'intensity bin','intensity bins','ibins','ibin','i bins','i bin'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    if size(dataloc.d{xy(sb)}.data.(channel{iChan}),1) > 4
                                        h{sc} = ct_trackvis_MHY(AH, 'histogram', dataloc.d{xy(sb)}.data.(channel{iChan})(1:numcells,firsttp:tracklength),'cmap',parula,'nolabel',true);
                                        set(gca, 'YDir','reverse','YColor','k');
                                        ylabel(strrep(channel{iChan},'_',' '),'Color',PlotRules(iChan).CmapColor);
                                    end
                                end
                                p.plottype{iPlot} = 'intensity bin';

                            case {'heatmap','heatmaps'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    h{sc} = ct_trackvis_MHY(AH, 'heatmap', dataloc.d{xy(sb)}.data.(channel{iChan})(1:numcells,firsttp:tracklength),'cmap',cmapcol,'nolabel',true);
                                    colorbar(AH,'off');
                                end
                                p.plottype{iPlot} = 'heatmap';
                                axis tight
                                % set(AH, 'YLim', [0.5 40.5]); %set for
                                % limited number of cells
                                set(AH, 'YDir', 'reverse');  % or 'normal'
                            case {'spatialheatmap', 'spatial heatmap','spatialheatmaps', 'spatial heatmaps'} % A GIFT FROM RAM
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    h{sc} = spatialheatmap(gca,dataloc,xy(sb),channel{iChan},cmapcol,p,firsttp,tracklength);
                                    colorbar(AH,'off');
                                    axis tight
                                end
                                p.plottype{iPlot} = 'spatialheatmap';

                            case {'sorted heatmap','sorted heatmaps'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    pDat = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    [~,I]=sort(sum(isnan(pDat),2));
                                    pDat = pDat(I,:);
                                    h{sc} = ct_trackvis_MHY(AH, 'heatmap', pDat(1:numcells,:),'cmap',cmapcol,'nolabel',true);
                                    %AH.XLim = [0, tracklength]; AH.YLim = [0, numcells];
                                    colorbar(AH,'off');
                                    axis tight
                                end
                                    p.plottype{iPlot} = 'sorted heatmap';
                                    clear pDat;

                            case {'mean sorted heatmap','mean sorted heatmaps'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    pDat = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    mDat = mean(pDat,2,"omitnan");
                                    [~,I]=sort(mDat,'descend');
                                    pDat = pDat(I,:);
                                    h{sc} = ct_trackvis_MHY(AH, 'heatmap', pDat(1:numcells,:),'cmap',cmapcol,'nolabel',true);
                                    colorbar(AH,'off');
                                    axis tight
                                end
                                    p.plottype{iPlot} = 'mean sorted heatmap';
                                    clear pDat mDat;

                            case {'delta heatmap','delta heatmaps'}
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    pDat = dataloc.d{xy(sb)}.data.(channel{iChan})(:,firsttp:tracklength);
                                    
                                    datMin = min(pDat,[],2,"omitnan");
                                    datMax = max(pDat,[],2,"omitnan");

                                    [~,I]=sort(sum(isnan(pDat),2));
                                    pDat = pDat(I,:);
                                    h{sc} = ct_trackvis_MHY(AH, 'heatmap', pDat(1:numcells,:),'cmap',cmapcol,'nolabel',true);
                                    %AH.XLim = [0, tracklength]; AH.YLim = [0, numcells];
                                    colorbar(AH,'off');
                                    axis tight
                                end
                                    p.plottype{iPlot} = 'delta heatmap';
                                    clear pDat pDelt pLength I

                            case 'mean vs if'
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    trackhold = []; ifhold = []; keepem = [];
                                    trackhold = nanmean(dataloc.d{xy(sb)}.data.(channel{iChan})(:,tracklength-p.iftrackwindow:tracklength),2); 
                                    ifhold = dataloc.IFd{xy(sb)}.data.(p.ifchan);
                                    keepem = all([~isnan(trackhold), ~isnan(ifhold)],2);
                                    trackhold = trackhold(keepem); ifhold = ifhold(keepem);
                                    h{sc} = scatter(trackhold,ifhold,15,[0,0,1],'filled');
                                    P = polyfit(trackhold,ifhold,1); %fit a slope to the data
                                    yfit = P(1)*trackhold+P(2);  % P(1) is the slope and P(2) is the intercept
                                    hold on;
                                    plot(trackhold,yfit,'k-.','LineWidth',2);
                                    text(max(trackhold)*0.8,max(ifhold)*0.8,num2str(P(1)));
                                    ylabel([strrep(p.ifchan,'_',' '), ' signal'])
                                    xlabel([strrep(channel{iChan},'_',' '), ' signal'])
                                    hold off;
                                end
                            
                            case 'slope vs if'
                                for iChan = 1:NumChans
                                    cmapcol = PlotRules(iChan).CmapColor;
                                    trackhold = []; ifhold = [];
                                    trackhold = dataloc.d{xy(sb)}.data.(channel{iChan})(:,tracklength-p.iftrackwindow:tracklength); 
                                    ifhold = dataloc.IFd{xy(sb)}.data.(p.ifchan);
                                    hasIFd = ~isnan(ifhold);
                                    trackhold = trackhold(hasIFd,:); ifhold = ifhold(hasIFd);
                                    [numCells,numTps] = size(trackhold);
                                    linez = NaN(numCells,1);
                                    timz = 1:numTps;
                                    for ii = 1:numCells
                                        P=polyfit(timz,trackhold(ii,:),1);
                                        linez(ii,1) = abs(P(1));
                                    end
    
                                    h{sc} = scatter(ifhold,linez,15,[0,0,1],'filled');
                                    xlabel([strrep(p.ifchan,'_',' '), ' signal'])
                                    ylabel([strrep(channel{iChan},'_',' '), ' signal slope'])
                                    %hold off;
                                end
                        end %switch
                    
                    if firstxy
                        if ~contains(p.plottype{iPlot}, ["vs if","mean vs mean"])
                        legvec = [legvec,1]; % identify unique legend entries.
                        legtxt = cat(1,titlevec,legname{sa});
                        lineTxs = txs;
                        if p.plotfromzero
                            for ixline = 1:numel(lineTxs)
                                if contains(channel{1},"Corr",'IgnoreCase',true)
                                elseif contains(p.plottype{iPlot},"heatmap") 
                                    xline(AH,lineTxs(ixline)-firsttp+1,'--r','LineWidth', 1);
                                else
                                    xline(AH,lineTxs(ixline)-firsttp+1,'--','LineWidth', 1);
                                end
                            end %xline
                        else
                            for ixline = 1:numel(lineTxs)
                                if contains(channel{1},"Corr",'IgnoreCase',true)
                                elseif contains(p.plottype{iPlot},"heatmap") 
                                    xline(AH,lineTxs(ixline),'--r','LineWidth', 1);
                                else
                                    xline(AH,lineTxs(ixline),'--','LineWidth', 1);
                                end
                            end %xline
                        end %plot from zero
                        end %not an IF plot
                        % axis and title      
                        tit = legname{sa};
                        if p.combinexys || (p.plottype{iPlot} == "mean" && p.overlapmeans)
                        else
                            tit = [tit, ' xy', num2str(xy(sb))];
                        end
                        title(tit);
                        

                    else
                        legvec = [legvec,0];
                    end %first xy
                    
                    if ~contains(channel{1},"Corr",'IgnoreCase',true) && ~contains(p.plottype{iPlot},["vs if","mean vs mean"])
                        findZero(txs,thisTX,tracklength,firsttp,p)
                    elseif contains(channel{1},"Corr",'IgnoreCase',true); xlabel('');
                    end %don't shift the axis
                    else
                        warning(['XY ',num2str(xy(sb)),' contains no data.'])
                    end
                    end% field check
                end %sb 
            end %sa

        % Make the new file name
        Chanz = [];
        for iChan = 1:numel(channel); iChannel = channel{iChan}; 
            if startsWith(iChannel,'n') && ~endsWith(iChannel, 'Corr')
                if ischar(dataloc.normalizetype); iChannel = [iChannel, ' ', dataloc.normalizetype]; 
                else; iChannel = [iChannel, ' ', dataloc.normalizetype.(iChannel)]; 
                end
            end
            Chanz = [Chanz, ' ', iChannel];  
        end
        
        if p.plottype{iPlot} == "mean" && p.overlapmeans; ChanSaveName = ['plotbyMHY',Chanz,' overlaping means'];
            else; ChanSaveName = ['plotbyMHY',Chanz,' ', p.plottype{iPlot}]; end
        if ~isempty(p.ncells); ChanSaveName = [ChanSaveName, '_ncells', num2str(p.ncells)]; end
        if ~isempty(p.tmaxaftertx); ChanSaveName = [ChanSaveName, ' ' num2str(p.tmaxaftertx/p.tktm) ' hrs after tx']; end
        if ~isempty(p.subset); for iSubset = 1:numel(p.subset); ChanSaveName = [ChanSaveName, ' ', p.subset{iSubset}]; end; end
        if ~isempty(p.exclude); ChanSaveName = [ChanSaveName, ' ', 'ex']; for iSubset = 1:numel(p.exclude); ChanSaveName = [ChanSaveName, ' ', p.exclude{iSubset}]; end; end
        if p.combinexys; ChanSaveName = [ChanSaveName, ' combxys']; end
        if numel(facetby) > 1; ChanSaveName = [ChanSaveName,' ',facetby{s}]; end
        ChanSaveName = strrep(ChanSaveName,'_',' ');
        ChanSaveName = strrep(ChanSaveName,'.',' pt ');

        if endsWith(dataloc.fold.fig,'\'); FullSaveName = [dataloc.fold.fig,ChanSaveName]; 
        else; FullSaveName = [dataloc.fold.fig,'\',ChanSaveName]; 
        end

        % Give the figure a title
        tittxt = [dataloc.file.base, ' ', ChanSaveName];
        sgtitle(strrep(tittxt,'_',' '));

        if ~exist('h','var') && contains(p.plottype{iPlot}, "sorted stacks v2")
            h = get(gcf, 'Children');
            h = h.Children;
        end
        
        if (~contains(p.plottype{iPlot}, "intensity bin")) && exist('h','var') && p.standardizeplots && ~isempty(figgy) %&& ~isempty(legtxt)
            subplot_standardizer(figgy, p.plottype{iPlot},p,channel) % make axis limits the same for all plots
            if (~contains(p.plottype{iPlot}, "intensity bin")) && exist('h','var') && p.standardizeplots && ~isempty(figgy)
                subplot_standardizer(figgy, p.plottype{iPlot}, p, channel)
            
                if contains(p.plottype{iPlot},"heatmap")
                    % Use the last heatmap axes you actually plotted into (store it!)
                    if exist('AH','var') && ~isempty(AH) && isgraphics(AH,'axes')
                        extraleg(AH, channel{iChan});
                    end
                end
            end
        end

        fontsize(gcf,p.font_size,'points')
        if ~isempty(get(gcf, 'Children')) %dont save an empty figuer
        switch p.printstyle
            case 'pdf'
                set(gcf,'visible','off','Resize','off','Renderer', 'painter')
                set(gcf,'units','inches','PaperSize', [(3*figgy.GridSize(2)), (3*figgy.GridSize(1))],'position', [0.5,0.5,(3*figgy.GridSize(2))+0.5, (3*figgy.GridSize(1))+0.5],'OuterPosition',[0.5,0.5,(3*figgy.GridSize(2))+0.5, (3*figgy.GridSize(1))+0.5]) %[0,0,figgy.GridSize(2)*270,figgy.GridSize(1)*270])
                fontsize(gcf,p.font_size,"points")
                fontname(gcf,'Calibri')
                exportgraphics(gcf, [FullSaveName,'.pdf'],'Resolution', 300 ) % 300 dpi printing
            case 'svg'
                set(gcf,'Renderer', 'painter','units','inches') % ,'normalized','outerposition',[0 0 1 1])
                set(gcf,'Position',[0.5,0.5,2*figgy.GridSize(2)+0.5, (2*figgy.GridSize(1))+0.5],"Papersize",[2*figgy.GridSize(2)+1, (2*figgy.GridSize(1))+1],'Resize',false) %
                fontsize(gcf,p.font_size,"points"); fontname(gcf,'Calibri');
                if ~isempty(findall(gcf,'Type','axes'))
                    saveas(gcf, FullSaveName, 'svg')
                end
            case 'test'
                set(gcf,'Position',[0,50,360*figgy.GridSize(2), (360*figgy.GridSize(1))+50],'Renderer', 'painter') %
                saveas(gcf, FullSaveName, 'svg')
        end

            if exist('printTable','var') %print the table for the data
                writetable(printTable,[FullSaveName,'_data'],'FileType','spreadsheet')
                clear printTable
            end
        else; warning([dataloc.file.base, ChanSaveName, ' has no data, so no plots will appear.']) %warn people of empty 
        end
        if p.closefigs; close(gcf); end
        end %facetby loop
        clear h;
end %iplot loop
end     
%% function to set same x and y axis limits for all plots
function [] = subplot_standardizer(fig, plottype, p, channel)

sp = get(fig, 'Children');

if contains(plottype, "heatmap")
    plottype = 'heatmap';
elseif contains(plottype, ["stacks","mean vs mean"])
    sp = findall(fig,'Type','Axes');
end

if ~isempty(sp)
switch plottype
    case {'heatmap', 'spatialheatmap'}
        for i = 1:numel(sp)
            clim_min(i) = sp(i).CLim(1);
            clim_max(i) = sp(i).CLim(2);
        end
        if ~isempty(p.ymn) && ~isempty(p.ymx)
            cvals(1) = p.ymn; cvals(2) = p.ymx;
        else
            cvals(1) = min(clim_min); cvals(2) = max(clim_max);
            % scale it down a little
            cvals(1) = cvals(1)*1.05; cvals(2) = cvals(2)*0.95;
        end

        if ~isempty(p.nancolor); cmp = [p.nancolor; colormap]; % Take the colormap and make it have nans be black
        else; cmp = colormap; % Take the colormap
        end

        for i = 1:numel(sp)
            if iscell(cvals); cvals = cell2mat(cvals); end
           set(sp(i),'CLim', cvals);
           set(sp(i),'Colormap',cmp);
        end

        if any(contains(channel,"corr", "IgnoreCase",true),"all")
            xLims = sp(1).XLim;
            midPt = (((xLims(2)-xLims(1))-1)/2);
            tickz = [-midPt,0,midPt];
            for i = 1:numel(sp)
                set(sp(i),'XTick', [xLims(1),(((xLims(2)-xLims(1)))/2)+0.5,xLims(2)]);
                set(sp(i),'XTickLabel',tickz);
                set(sp(i),'XLim',xLims)
            end
        end

    case 'hist3'
        for i = 1:numel(fig)
            clim_min(i) = fig(i).CLim(1);
            clim_max(i) = fig(i).CLim(2);
        end
        
        cvals(1) = min(clim_min); cvals(2) = max(clim_max);
        
        for i = 1:numel(fig)
            set(fig(i),'CLim', cvals);
        end

    case 'sorted stacks v2'
            Y = get(sp,'YAxis');
            Y = [Y{:}];
            if numel(channel) > 1
                 if ~isempty(p.ymn) && ~isempty(p.ymx)
                    set(Y(1,:),'Limits', [p.ymn{1}, p.ymx{1}]);
                    set(Y(2,:),'Limits', [p.ymn{2}, p.ymx{2}]);
                 else
                    lRange = [Y(1,:).Limits];
                    rRange = [Y(2,:).Limits];
                    set(Y(1,:),'Limits', [min(lRange(:)), max(lRange(:))]);
                    set(Y(2,:),'Limits', [min(rRange(:)), max(rRange(:))]);
                 end
            else 
                if ~isempty(p.ymn) && ~isempty(p.ymx)
                    set(Y(1,:),'Limits', [p.ymn{1}, p.ymx{1}]);
                else
                    Range = [Y(:).Limits]; set(Y(:),'Limits', [min(Range(:)), max(Range(:))]);
                end
            end
            
    case 'mean'
            Y = get(sp,'YAxis');
            if iscell(Y) && numel(Y) > 1; Y = [Y{:}]; end
            if numel(channel) > 1
                lRange = [Y(1,:).Limits];
                rRange = [Y(2,:).Limits];
                set(Y(1,:),'Limits', [min(lRange(:)), max(lRange(:))]);
                set(Y(2,:),'Limits', [min(rRange(:)), max(rRange(:))]);
            else; Range = [Y(:).Limits]; set(Y(:),'Limits', [min(Range(:)), max(Range(:))]);
            end
            
    otherwise
        if ~contains(plottype,'stacks')
            for i = 1:numel(sp)
                xlim_min(i) = sp(i).XLim(1);
                xlim_max(i) = sp(i).XLim(2);
                ylim_min(i) = sp(i).YLim(1);
                ylim_max(i) = sp(i).YLim(2);
            end
        
            xvals(1) = min(xlim_min); xvals(2) = max(xlim_max);
            yvals(1) = min(ylim_min); yvals(2) = max(ylim_max);
        
            for i = 1:numel(sp)
                if ~isempty(p.ymn) && ~isempty(p.ymx) && ~contains(plottype,"stacks"|"mean vs mean")
                    set(sp(i),'YLim',[p.ymn,p.ymx])
                elseif ~(plottype == "stacks")
                    set(sp(i),'YLim', yvals);
                end
                set(sp(i),'XLim', xvals);
            end
        end

end

end %empty plot check
end
%% function for placing legend in extra plot
function extraleg(ax, chan)
% extraleg - attach a colorbar to a specific axes (tiledlayout-friendly)

    if nargin < 2 || isempty(chan); chan = ''; end
    if isempty(ax) || ~isgraphics(ax,'axes')
        return; % nothing plotted -> no legend/colorbar
    end

    % Create/update the colorbar for this axes
    cb = colorbar(ax);
    cb.Label.String = chan;

    % If we're in a tiledlayout, allow it to dock on the east
    try
        cb.Layout.Tile = 'east';
    catch
        % Older MATLAB / non-tiledlayout contexts: ignore
    end
end

%% Colormap
function [cmap] = rainbow(n)
values = [
213,62,79
244,109,67
253,174,97
254,224,139
171,221,164
102,194,165
50,136,189
94,79,162]./256;

cmap = interp1(1:size(values,1), values, linspace(1,size(values,1),n), 'linear');
end
%% PlotRulesDefiner
function PlotRules = PlotFormatRules(channels, plottype) %FIX TO ADD OPTION FOR NON-NORMALIZED DATA LABELS

% PlotRulesLines is set up like this: 'Channel to look for', 'Y label', 'cmap color', regexpi settings;

PlotRulesList = {...
    'Signal Cross Corr',    'Cross correlations',                           '#C30F0E',       '(Corr)$'; ...
    'RAMPKAR',              'AMPK Activity (Red-Arbitrary Units)',          '#C30F0E',       '^(RAMPKAR2)$'; ...
    'RAMPKAR_Low',          'AMPK Activity (Red-Arbitrary Units)',          'rgb(25 118 210)',      '^(RAMPKAR2_Low)$'; ...
    'Hylight',              'FBP levels',                                   '#1976D2',     '(n)?(HYLIGHT)'; ...
    'Perceval',             'ATP-to-ADP Ratio',                             '#168039',     '(n)?(PERCEVAL)'; ...
    'AMPK',                 'AMPK Activity (Arbitrary Units)',              '#04BFBF',      '^(AMPKAR)$'; ...
    'ERK',                  'ERK (EKAR) Activity (Arbitrary Units)',        'magenta',   '^(nEKAR)'; ... 
    'ERK',                  'ERK (EKAR) Activity (Arbitrary Units)',        'magenta',   '^(EKAR)'; ... 
    'ERK',                  'ERK (EKAR) Activity (Arbitrary Units)',        '#C30F0E',       '^(REKAR)'; ... 
    'ERK',                  'ERK (ERKTR) Activity (Arbitrary Units)',       '#C30F0E',       '(n)?ERKTR'; ...
    'JNK',                  'JNK (JNKTR) Activity (Arbitrary Units)',       '#D17300',       '(n)?JNKTR'; ...    
    'AKT',                  'AKT (EXRAI AKT) Activity (Arbitrary Units)',   '#168039',     '(n)?AKT(2)?KTR'; ...
    'Granularity',          'Granularity',                                  '#C30F0E',       '^(Granularity)'; ...
    'P65KTR',               'NFkB Activity (Arbitrary Units)',              '#168039',     '(n)?P65(KTR)?'; ...
    'mCard',                'mCard Intensity (Arbitrary Units)',            '#C30F0E',       'mCard_Nuc'; ...
    'Med_ERK',              'Median ERK Activity (Arbitrary Units)',        '#C30F0E',       '(n)?Med_E(R)?K(TR|AR)'; ... 
    'Med_P65KTR',           'Median NFkB Activity (Arbitrary Units)',       '#168039',     '(n)?Med_P65(KTR)?'; ...
    'Velocity',             'Instantaneous Speed',                          '#168039',     '^(Velocity)'...
    };
 
% Check to see if this plot has formating rules
if ischar(channels) && ~iscell(channels); channels = {channels}; end

PlotRulesIdxs = cellfun(@(x)find(~cellfun(@isempty,regexpi(x, PlotRulesList(:,4)))), channels, 'UniformOutput',false);
NumChans = numel(channels);
for ii = 1:NumChans
    if isempty(PlotRulesIdxs{ii})
        %fprintf('Could not find plot rules for %s, using defaults. \n', channels{ii})
        PlotRules(ii).Channel = channels{ii}; %#ok<*AGROW>
        PlotRules(ii).Ylabel = [channels{ii},' Activity (Arbitrary Units)'];
        PlotRules(ii).CmapColor = 'red';
    else
        PlotRules(ii).Channel = PlotRulesList{PlotRulesIdxs{ii},1};
        PlotRules(ii).Ylabel = PlotRulesList{PlotRulesIdxs{ii},2};
        PlotRules(ii).CmapColor = PlotRulesList{PlotRulesIdxs{ii},3};
    end
if contains(plottype,'heatmap'); PlotRules(ii).CmapColor = 'parula'; end
end

end
%% Plots the Pulse Analysis Data by Treatment (can combine xys per treatment)
function Pulse_Plotter(dataloc,channel,facetby,groupby,idx,xymat,linetp,titlevec,goodxy,p)
if ~iscell(channel) && ischar(channel); channel = {channel}; end
ThisChan = channel{1};
if numel(channel) > 1; ThisChan2 = channel{2}; end
titlevec = arrayfun(@(x)regexprep(titlevec{x},'/',' '),1:numel(titlevec),'un',0)';
% remove any underscores in legend text to avoid subscript shenanigans
TXX=[];

%% actual plotting part
for iPlot = 1:numel(p.analysischan)
     for s = 1:numel(facetby)
            legvec = []; legtxt = {};  groupings = struct('layer1',[],'layer2',[],'nums',[],'origin',[],'iszero',[],'respgrps',[],'data',[]); PlotSaveName = []; lowestValue = Inf;
            for sa = 1:numel(groupby) % plot line for each unique treatment/cell
                xy = intersect(idx.(facetby{s}),idx.(groupby{sa})); % find xy(s) for treatment/celltype
                xys4sub = length(xy); sc = 0;
                for sb = 1:xys4sub % collect each xys data
                    ThisXY = xy(sb);
                    if goodxy(xy(sb)) && isfield(dataloc.d{xy(sb)}.data,(ThisChan))% make sure there is data to plot.
                        sc = sc+1; % index for saving group2 into.
                        
                        treatment4XY = str2double(linetp{cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0))});
                        anyunique = unique(treatment4XY);
                        txs=anyunique(~isnan(anyunique)); txs = txs(txs>0); %get the treatment times not including pretreatments
                        if ~isempty(txs)
                        if numel(txs) < p.aftertreatment; thisTx=txs(end); else; thisTx=txs(p.aftertreatment); end %get treatment number or last treatment
                        else; thisTx = 1;
                        end
                        movieLength = size(dataloc.d{xy(sb)}.data.(ThisChan),2);
                        if ~isempty(p.tstartaftertx); tStart = p.tstartaftertx+thisTx; else; tStart = thisTx; end
                        if tStart > movieLength; tStart = movieLength; end % if the movie isn't that long... don't get the data
                        if isempty(p.tmaxaftertx); tEnd = movieLength; else; tEnd = p.tmaxaftertx+thisTx-1; end
                        if tEnd > movieLength; tEnd = movieLength; end
                        goodData = ~isnan(dataloc.d{ThisXY}.data.(ThisChan)(:,(tStart:tEnd)));

                        switch  p.analysischan{iPlot}
                            case {'frequency', 'freq'}
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                a = cellfun(@(x) sum(tEnd>x & x>tStart),thisData,'UniformOutput',false); 
                                a=cell2mat(a); %a = a(all(goodData,2));
                                Div = (sum(goodData,2)/ p.tktm); %Div = Div(all(goodData,2)); 
                                a = a ./ Div;
                                moreThan1 = unique(a); moreThan1 = sort(moreThan1); %janky way of getting the number for just 1 pulse
                                if any(moreThan1>0)
                                    moreThan1 = moreThan1(moreThan1 > 0); moreThan1 = moreThan1(1); 
                                    if moreThan1<lowestValue; lowestValue = moreThan1; end
                                end
                                catRange = [0,lowestValue,lowestValue+0.000001,inf]; catGroups = {'no response','one pulse','> 1 pulse'};
                                perXLabel = 'Responder distribution';
                                tite = [ThisChan,' frequency over time (pulses per hour)'];
                                PlotSaveName = 'Pulse frequency';
                                p.analysischan{iPlot} = 'freq';

                            case {'frequency v2', 'freq v2'}
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                a = cellfun(@(x) sum(tEnd>x & x>tStart),thisData,'UniformOutput',false); 
                                a=cell2mat(a); a = a(all(goodData,2));
                                Div = (sum(goodData,2)/ p.tktm); Div = Div(all(goodData,2)); a = a ./ Div;
                                catRange = [0,p.pulsethresh-0.000001,inf]; catGroups = {'inactive',[' active (> ',num2str(p.pulsethresh) ,' pulses per hour)']};
                                perXLabel = 'Responder distribution';
                                tite = [ThisChan,' frequency over time (pulses per hour)'];
                                PlotSaveName = 'Pulse frequency v2';
                                p.analysischan{iPlot} = 'freq v2';

                            case {'mean','means'}
                                a = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tStart:tEnd)),2,'omitnan');
                                a = a(~isnan(a));
                                tite = [ThisChan,' mean activity'];
                                catRange = [0,p.responsethresh,Inf]; catGroups = {'inactive',['active (signal > ',num2str(p.responsethresh),')']};
                                perXLabel = 'active cell distribution';
                                PlotSaveName = 'Mean activity';
                                p.analysischan{iPlot} = 'mean';

                            case {'delta mean','delta means'}
                                afterTx = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tStart+1:tEnd)),2,'omitnan');
                                tpBack = tStart-floor(p.tmaxback*p.tktm); % what tp to go back to after convert hours to tps
                                if tpBack < 1; tpBack = 1; end % if you try to go back before time exists, stop you
                                beforeTx = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tpBack:tStart-1)),2,'omitnan');
                                deltaMean = afterTx - beforeTx; % get the diff in means
                                a = deltaMean(~isnan(deltaMean));
                                tite = [ThisChan,' delta mean activity ', num2str(p.tmaxback), 'h before tx up to ', num2str(tEnd/p.tktm),' h after tx'];
                                tite = strrep(tite,'_', ' ');
                                catRange = [0,p.responsethresh,Inf]; catGroups = {'non-responder',['responder (delta change > ',num2str(p.responsethresh),')']};
                                perXLabel = 'percent responder distribution';
                                PlotSaveName = 'delta mean activity';
                                p.analysischan{iPlot} = 'delta mean';

                            case {'perc responder','perc responders'}
                                afterTx = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tStart+1:tEnd)),2,'omitnan');
                                tpBack = tStart-floor(p.tmaxback*p.tktm); % what tp to go back to after convert hours to tps
                                if tpBack < 1; tpBack = 1; end % if you try to go back before time exists, stop you
                                beforeTx = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tpBack:tStart-1)),2,'omitnan');
                                deltaMean = afterTx - beforeTx; % get the diff in means
                                a = deltaMean(~isnan(deltaMean));
                                tite = [ThisChan,' delta mean activity ', num2str(p.tmaxback), 'h before tx up to ', num2str(tEnd/p.tktm),' h after tx'];
                                tite = strrep(tite,'_', ' ');
                                catRange = [0,p.responsethresh,Inf]; catGroups = {'non-responder',['responder (delta change > ',num2str(p.responsethresh),')']};
                                perXLabel = 'percent responder distribution';
                                PlotSaveName = 'delta mean activity';
                                p.analysischan{iPlot} = 'perc responder';

                            case {'ifmeans'}
                                a = mean(dataloc.d{ThisXY}.data.(ThisChan)(:,(tStart:tEnd)),2);
                                a = a(~isnan(a));
                                tite = [ThisChan,'if mean activity'];
                                catRange = [0,p.responsethresh,inf]; catGroups = {'inactive',['active (signal > ',num2str(p.responsethresh),')']};
                                perXLabel = 'Responder distribution';
                                PlotSaveName = 'IF Mean activity';
                                p.analysischan{iPlot} = 'ifmeans';

                            case {'peakdurs', 'durs', 'duration'}
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                b = cellfun(@(x)(tEnd>x & x>tStart),thisData,'UniformOutput',false); 
                                b = cell2mat(b);
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).dur}';%get the durs
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                a = cell2mat(thisData); a = a(b) / p.tktm; % get the right durs get the durs over hours
                                catRange = [-inf,0.000001,p.lengththresh,inf]; catGroups = {'no response',['short (< ',num2str(p.lengththresh),' hr) pulses'],['long (>',num2str(p.lengththresh),'hr) pulses']};
                                perXLabel = 'Response distribution';
                                tite = ['Duration of ', ThisChan,' peaks per condition (in hours)'];
                                PlotSaveName = 'Duration of Peaks';
                                p.analysischan{iPlot} = 'durs';

                            case {'dur x response', 'durxresponse'}
                                theseCells = all(goodData,2); % get the cells that have tracks in the span of question
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}'; %get the peak data
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0}; %fill the empties
                                thisData = thisData(theseCells); % keep the data from cells in question
                                a = nan(numel(thisData),1); %make an 0 time reponse array (for % responders later)

                                responders = cellfun(@(x) sum(tEnd>x & x>(tStart-1)),thisData); % get the cells that respond in the given time
                                responders = responders >= p.xpulse; % get the response you care about (1st, 2nd, etc, response)
                                responseTimes = thisData(responders); % pull the times of the responses
                                properResponse = cellfun(@(x) x>(tStart-1), responseTimes,'UniformOutput',false); % establish which responses are in the right window
                                
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).dur}'; %get the dur data that we care about
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0}; %fill the empties
                                thisData = thisData(theseCells); % keep the data from cells in question
                                thisData = thisData(responders); %keep the responding cell data

                                properResponse = cellfun(@(x,y) x(y), thisData,properResponse,'UniformOutput',false); % pull durs of responses in the right window
                                durXResponse = cellfun(@(x) x(p.xpulse), properResponse); % get the xpulse dur of proper responders
                                durXResponse = durXResponse / p.tktm; % a is time to first response (make it into minutes)
                                a(responders) = durXResponse; %put the dur times in the 0 array
                                %a = durXResponse; %give the durations of x response (within time span)
                                %sum(b > 0)/size(b,1);
                                catRange = [-inf,0,inf]; catGroups = {'no response','responder'};
                                perXLabel = 'Response distribution';
                                tite = ['Duration of ', ThisChan,' ', num2str(p.xpulse)  ,' peak to treatment (in hours)'];
                                PlotSaveName = ['Duration of ', num2str(p.xpulse)  ,' response'];
                                p.analysischan{iPlot} = 'durxresponse';

                            case 'numshortpeaks'
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                a = cellfun(@(x)(tEnd>x & x>tStart),thisData,'UniformOutput',false); 
                                a = cell2mat(a);
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).dur}';%get the durs
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                b = cell2mat(thisData); b = b(a) / p.tktm; % get the right durs get the durs over hours
                                shortPs = (b / p.tktm) <= p.lengththresh;
                                a = a(shortPs);
                                tite = [' Number of short (< ', num2str(p.lengththresh),' hr) ', ThisChan,' peaks '];
                                PlotSaveName = 'Number of short peaks ';
                                
                            case {'time2pulse', 'time2xpulse'} % fix - add AND PERCENT RESPONDERS
                                thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                %a = zeros(numel(thisData),1);
                                responders = cellfun(@(x) sum(tEnd>x & x>tStart),thisData);
                                responders = responders >= p.xpulse;
                                responseTimes = thisData(responders);
                                properResponse = cellfun(@(x) x>tStart, responseTimes,'UniformOutput',false);
                                properResponse = cellfun(@(x,y) x(y), responseTimes,properResponse,'UniformOutput',false);
                                time2XResponse = cellfun(@(x) x(p.xpulse), properResponse);
                                responseTime = time2XResponse - thisTx;
                                responseTime = responseTime / p.tktm; % a is time to first response (make it into hours)
                                a = responseTime; % a(responders) = responseTime;
                                catRange = [0,0.0000001,p.lengththresh,inf]; catGroups = {'no response',['fast response (within ', num2str(p.lengththresh),' hours)'], 'slow response'};
                                perXLabel = 'Response distribution';
                                tite = [ThisChan,' time (hours) to ',num2str(p.xpulse),' pulse response to ligand'];
                                PlotSaveName = ['Time to ' num2str(p.xpulse),' response'];

                            case 'roughdivisions'
                                b=[];
                                for ii=1:NumXys
                                    ThisXY = XyNums(ii);
                                    treatment4XY = str2double(linetp{cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0))});
                                    anyunique = unique(treatment4XY); 
                                    lsttx=anyunique(~isnan(anyunique)); lsttx = lsttx(lsttx>0); %get the treatment times not including pretreatments
                                    if numel(lsttx) < p.aftertreatment; tStart=lsttx(end); else; tStart=lsttx(p.aftertreatment); end %get treatment number or last treatment
                                    nCellsStart = sum(~isnan(dataloc.d{ThisXY}.data.XCoord(:,tStart)));
                                    if size(dataloc.d{ThisXY}.data.XCoord,2)<(p.divtime+tStart); p.divtime = []; end
                                    if ~isempty(p.divtime)
                                        nCellsEnd = sum(~isnan(dataloc.d{ThisXY}.data.XCoord(:,(p.divtime+tStart))));
                                        movieLength = p.divtime / p.tktm;
                                    else
                                        nCellsEnd = sum(~isnan(dataloc.d{ThisXY}.data.XCoord(:,end)));
                                        movieLength = (size(dataloc.d{ThisXY}.data.XCoord,2) - tStart) / p.tktm;
                                    end
                                    b(ii) = ((nCellsEnd - nCellsStart) / nCellsStart)*100;
                                end
                                tite = ['Percent change in cell number over ', num2str(movieLength) ,' hours (~divisions)'];
                                %yLabe = 'Change in number cells';
                                % how many hours was the movie
                                % Make the new file name
                                ThisChan = 'Perc_Divisions_over';
                                PlotSaveName = [num2str(movieLength) ,'_hours'];
                      
                            case 'meandurspercell'
                                b=[];
                                for ii=1:NumXys
                                    try
                                        ThisXY = XyNums(ii);
                                        a=cell2mat({dataloc.z{ThisXY}.data.(ThisChan).dur}');
                                        b(ii)=numel(a);
                                    catch
                                        b(ii)=0;
                                    end
                                end
                                maxnumdurs = max(b);
                                clear b;
                                b = NaN(maxnumdurs,NumXys);
                                for ii=1:NumXys
                                    ThisXY = XyNums(ii);
                                    treatment4XY = str2double(linetp{cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0))});
                                    anyunique = unique(treatment4XY); 
                                    lsttx=anyunique(~isnan(anyunique)); lsttx = lsttx(lsttx>0); %get the treatment times not including pretreatments
                                    if numel(lsttx) < p.aftertreatment; lsttx=lsttx(end); else; lsttx=lsttx(p.aftertreatment); end %get treatment number or last treatment
                                    thisData = {dataloc.z{ThisXY}.data.(ThisChan).pkpos}';
                                    dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                    a = cellfun(@(x) x>lsttx,thisData,'UniformOutput',false); 
                                    thisData = {dataloc.z{ThisXY}.data.(ThisChan).dur}';%get the durs
                                    dEmpty = cellfun(@isempty, thisData); thisData(dEmpty) = {0};
                                    c = cellfun(@(x,y) x(y), thisData,a, 'UniformOutput',false);
                                    c = cellfun(@mean, c, 'UniformOutput',false);
                                    c = cell2mat(c);
                                    b(1:numel(c),ii)= c / p.tktm; %get the durs over hours
                                end
                                plotType = p.pptype;
                                tite = ['Mean Duration of ', PlotRules(iChan).Channel,' peaks per condition (in hours)'];
                                PlotSaveName = ['Mean Duration of Peaks ', plotType, ' '];
                        end% switch for plot type
                        
                        groupings.iszero = [groupings.iszero; ~(a > 0)];
                        groupings.data = [groupings.data; a]; 
                        if ~isempty(p.overlaptx)
                            %get matrix of overlaping treatments minus the treatment provided
                            SplitTxs = split(titlevec{sa},'+');
                            if (numel(SplitTxs) < p.overlaptx) %if not enough treatments, don't split them
                                groupings.layer1 = [groupings.layer1; repelem(string(titlevec{sa}),size(a,1))'];
                                groupings.layer2 = [groupings.layer2; repelem(" ",size(a,1))']; 
                            else
                                thisTx = false(numel(SplitTxs),1); thisTx(p.overlaptx) = 1;
                                otherTxs = {join(SplitTxs(~thisTx),'+')};
                                groupings.layer1 = [groupings.layer1; repelem(string(otherTxs),size(a,1))']; %give group by all other txs
                                groupings.layer2 = [groupings.layer2; repelem(string(SplitTxs{p.overlaptx}),size(a,1))']; %seperate out by picked tx                         
                            end
                        else
                            groupings.layer1 = [groupings.layer1; repelem(string(titlevec{sa}),size(a,1))'];
                            groupings.layer2 = [groupings.layer2; repelem(string(['rep ', num2str(sc)]),size(a,1))'];
                        end
                        groupings.nums = [groupings.nums; (ones(size(a,1),1) * sa)];
                        groupings.origin = [groupings.origin; repelem(string(titlevec{sa}),size(a,1))'];
                    end
                end %sb 
            end %sa
    if ~isempty(groupings.data)        
    groupings.respgrps = discretize(groupings.data,catRange,'categorical',catGroups);
    groupings.respgrps(isundefined(groupings.respgrps)) = catGroups{1};
    groupings.originwzeros = groupings.origin;
    if p.withoutzeros
        groupings.data = groupings.data(~groupings.iszero);
        groupings.layer1 = groupings.layer1(~groupings.iszero);
        groupings.layer2 = groupings.layer2(~groupings.iszero);
        groupings.nums = groupings.nums(~groupings.iszero);
        groupings.origin = groupings.origin(~groupings.iszero);
    end
    groups2Plot = unique(groupings.origin,'stable'); 

    % make the plot
    for iPlotType = 1:numel(p.plottype)
    plotType = p.plottype{iPlotType};
    switch plotType
        case 'violin'
            violinplot(groupings.data,groupings.origin,'GroupOrder',titlevec)
            xtickangle(90); ylabel(tite)
            depthNum = 1080;
        case {'compactboxplot'}
            boxplot(groupings.data,{groupings.layer1, groupings.layer2},'orientation','horizontal','PlotStyle','compact','ColorGroup',groupings.layer2)
            legend(findobj(gca,'Tag','Box'),flip(unique(groupings.layer2,'stable')))
            xlabel(tite)
            depthNum = numel(unique(groupings.layer1))*(90*numel(unique(groupings.layer2)));
        case {'box','boxplot'}
            boxplot_ND(groupings.data,{groupings.layer1, groupings.layer2},'orientation','horizontal','ColorGroup',groupings.layer2,'showalldata',p.c);
            xlabel(tite); xlim([0,inf])
            legend(findobj(gca,'Tag','Box'),flip(unique(groupings.layer2,'stable')),'location','best');
            depthNum = numel(unique(groupings.layer1))*(90*numel(unique(groupings.layer2)));
        case {'barh','bar'}
            barh(groupings.data);
            set(gca, 'YTick', [1:numel(titlevec)])
            set(gca, 'YTickLabel', titlevec)

        case {'scatter'}
            scatter(groupings.data);
            set(gca, 'YTick', [1:numel(titlevec)])
            set(gca, 'YTickLabel', titlevec)

        case {'bar_v2','barh_v2'}
            barHold = []; barHold.origin = groupings.origin; barHold.data = groupings.data;
            barHold = struct2table(barHold);
            barHold = grpstats(barHold,"origin",["mean","median",p.errortype]);
            barHold.origin = reordercats(categorical(barHold.origin),barHold.origin);
            barh(barHold.origin,barHold.mean_data,...
                'FaceColor',[0.75,0.75,0.75],...
                'EdgeColor','k', 'LineWidth',1.5);
            hold on;
            if size(barHold.([p.errortype,'_data']),2) > 1
                errorbar(barHold.mean_data,barHold.origin,barHold.([p.errortype,'_data'])(:,1),barHold.([p.errortype,'_data'])(:,2),'.','horizontal','Color','k','LineWidth', 1.5,'MarkerSize', 1);
            else
                errorbar(barHold.mean_data,barHold.origin,barHold.([p.errortype,'_data']),'.','horizontal','Color','k','LineWidth', 1.5,'MarkerSize', 1);

            end    
            if p.showrawdata
                swarmchart(groupings.data,categorical(groupings.origin),'XJitter','none','YJitter','density','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]) % ,'MarkerFaceColor','b'
            end

            if p.addmedian
                plot(barHold.median_data,barHold.origin,'d','MarkerFaceColor','r','MarkerEdgeColor','k') % ,'MarkerFaceColor','b'
            end

            if p.showrawdata; maxX = max(groupings.data); minX = min(groupings.data);
            else; maxX = (max(barHold.mean_data) + max(barHold.([p.errortype,'_data']))); minX = (min(barHold.mean_data) - max(barHold.([p.errortype,'_data'])));
            end

            if p.addstats && (size(barHold.origin,1) < 26)
                [~,~,stats] = anova1(groupings.data,groupings.origin,'off');
                if isempty(p.control)
                    nConds = size(stats.gnames,1);
                    for i = 1:nConds
                        printMe = strjoin([num2str(i) ' == ' string(stats.gnames{i})]);
                        fprintf('%s \n', printMe);
                    end
                    thisOne = input('Which condition is the control? \n');
                else; thisOne = p.control; 
                end
                [results,~,~,gnames] = multcompare(stats,"CriticalValueType","dunnett",'ControlGroup',thisOne,'Display','off'); 
%                 resultsTbl = array2table(results,"VariableNames", ["Group","Control Group","Lower Limit","Difference","Upper Limit","P-value"]);
%                 resultsTbl.("Group") = gnames(resultsTbl.("Group"));
%                 resultsTbl.("Control Group") = gnames(resultsTbl.("Control Group"));
                barHold.pValue = nan(size(barHold,1),1); % create a nan placeholder for the r squared 
                % swarmchart(resultsTbl,'P-Value','Group','*','MarkerEdgeColor','r','MarkerFaceColor',[1,0,0])
                for i = 1:size(results,1)
                    thisGroup = gnames{results(i,1)}; % get the group the stats are for
                    barHold{thisGroup,'pValue'} = results(i,6); % pull this data's pValue value and put it into the barHold data
                    barLoc = find(barHold.origin == thisGroup);
                    if results(i,6) < 0.05
                        text(maxX*1.05, barLoc, '\ast ','HorizontalAlignment','center','FontSize',8,'Color','red'); % , num2str(results(i,6))
                    else; text(maxX*1.05, barLoc, 'N.S.','HorizontalAlignment','center','FontSize',8,'Color','red')
                    end
                end
                text(maxX*1.05, results(1,2), 'Cntrl','HorizontalAlignment','center','FontSize',10,'Color','red')
            end
            xlim([(minX-((maxX*1.1)-maxX)),maxX*1.1])
            hold off;
            xlabel(tite);
            depthNum = numel(groups2Plot) * 270;

        case 'percbar'
            percResps = [];
            [~, histEdges] = histcounts(groupings.data,'BinMethod','sturges');
            tiledlayout(numel(groups2Plot), 2,'TileSpacing','compact','Padding','tight');
            for iGroup = 1:numel(groups2Plot)
                percResps(iGroup,:) = histcounts(groupings.respgrps(groupings.originwzeros == groups2Plot{iGroup}));
                subP(iGroup) = nexttile((iGroup*2));
                gHistData = histcounts(groupings.data(groupings.origin == groups2Plot{iGroup}),histEdges);
                gHistData = gHistData ./ sum(gHistData,2);
                histogram('BinEdges',histEdges,'BinCounts',gHistData)
                subP(iGroup).XTick=histEdges; subP(iGroup).YAxisLocation = 'right'; subP(iGroup).YAxis.FontSize = subP(iGroup).YAxis.FontSize*0.6;
                if iGroup < numel(groups2Plot); subP(iGroup).XTickLabel = []; subP(iGroup).YTickLabel = []; end
                if p.addthreshline && exist('catRange',"Var"); hold on; xline(catRange(2:end-1),'--r','LineWidth',2); hold off; end
            end
            xlabel(tite);
            percResps = flip(percResps);
            percResps = (percResps ./ sum(percResps,2))*100;            
            nexttile(1,[numel(groups2Plot),1])
            bb = barh(percResps,'stacked');
            for ibb = 1:size(bb,2)
                xtips = bb(ibb).XEndPoints;
                if ibb > 1; ytips = (bb(ibb-1).YEndPoints + (bb(ibb).YEndPoints))/2;
                else; ytips = (bb(ibb).YEndPoints)/2;
                end
                labelz = string(num2str(bb(ibb).YData','%.2f'))';
                alinz = 'center';
                text(ytips,xtips,labelz,'HorizontalAlignment',alinz,'VerticalAlignment','middle')
            end
            xlim([0, 100]); xlabel(perXLabel)
            ylim([0.5, numel(groups2Plot)+0.5])
            [subP(:).YLim] = deal([min([subP.YLim]),max([subP.YLim])]);
            yticks([1:numel(groups2Plot)])
            yticklabels(flip(groups2Plot))
            legend(catGroups,'location','best')
            depthNum = numel(groups2Plot) * 270;
            barHold = table;
            barHold.treatments = groups2Plot;
            for i = 1:numel(catGroups)
                barHold.(catGroups{i}) = percResps(:,i);
            end
            barHold.pValue = ones([size(groups2Plot,1),1]);
            clear subP;
            
        case 'roughhist'
            percResps = [];
            histEdges = unique(groupings.data)';
            tiledlayout(numel(groups2Plot), 2,'TileSpacing','compact','Padding','tight');
            for iGroup = 1:numel(groups2Plot)
                percResps(iGroup,:) = histcounts(groupings.respgrps(groupings.originwzeros == groups2Plot{iGroup}));
                subP(iGroup) = nexttile((iGroup*2));
                gHistData = [];
                thisData = groupings.data(groupings.origin == groups2Plot{iGroup});
                thisData = thisData(~isnan(thisData));
                histEdges = histEdges(~isnan(histEdges));
                for iHist = 1:numel(histEdges)
                    gHistData(iHist) = numel(thisData(thisData == histEdges(iHist)));
                end
                gHistData = gHistData ./ sum(gHistData,2);
                bar(histEdges,gHistData,0.95); xticks(histEdges); xticklabels(num2str(histEdges',2));
                if iGroup < numel(groups2Plot); subP(iGroup).XTickLabel = []; subP(iGroup).YTickLabel = []; end
                subP(iGroup).YAxisLocation = 'right'; subP(iGroup).YAxis.FontSize = subP(iGroup).YAxis.FontSize*0.6;
            end
            xlabel(tite);
            percResps = flip(percResps);
            percResps = (percResps ./ sum(percResps,2))*100;            
            nexttile(1,[numel(groups2Plot),1])
            bb = barh(percResps,'stacked');
            for ibb = 1:size(bb,2)
                xtips = bb(ibb).XEndPoints;
                if ibb > 1; ytips = (bb(ibb-1).YEndPoints + (bb(ibb).YEndPoints))/2;
                else; ytips = (bb(ibb).YEndPoints)/2;
                end
                labelz = string(num2str(bb(ibb).YData','%.2f'))';
                alinz = 'center';
                text(ytips,xtips,labelz,'HorizontalAlignment',alinz,'VerticalAlignment','middle')
            end            
            xlim([0, 100]); xlabel(perXLabel)
            ylim([0.5, numel(groups2Plot)+0.5])
            [subP(:).YLim] = deal([min([subP.YLim]),max([subP.YLim])]);
            yticks([1:numel(groups2Plot)])
            yticklabels(flip(groups2Plot))
            legend(catGroups,'location','best')
            depthNum = numel(groups2Plot) * 270;
            clear subP;       
        
        case 'hist'
            [~, histEdges] = histcounts(groupings.data);
            [fitData, gNames] = fitdist(groupings.data,'Normal','By',{groupings.origin});
            plottingData = [];
            for iFit = 1:numel(fitData)
                 pdfData = pdf(fitData{iFit},histEdges);
                 plottingData(:,iFit) = pdfData / sum(pdfData);
            end
            imagesc(plottingData');
            set(gca, 'YTick', [1:numel(gNames)]); set(gca, 'YTickLabel', gNames)            
            set(gca, 'XTick', [1:numel(histEdges)]); set(gca, 'XTickLabel', histEdges)

            
    end
    
    % hide and Resize the plot
    %figgy = gcf;
    %set(figgy,'Position',[0,0,1920,(depthNum)],'Renderer', 'painter')%,'visible','off'
    % Make the new file name
    ChanSaveName = ThisChan;
    ChanSaveName = [ChanSaveName, ' ', PlotSaveName];
    if p.withoutzeros; tite = [tite, ' without zeros']; ChanSaveName = [ChanSaveName, ' without zeros']; end
    if ~isempty(p.tstartaftertx); ChanSaveName = [ChanSaveName, ' after ' num2str(p.tstartaftertx/p.tktm) ' hours']; end
    if ~isempty(p.tmaxaftertx); ChanSaveName = [ChanSaveName, ' before ' num2str(p.tmaxaftertx/p.tktm) ' hours']; end
    if ~isempty(p.ncells); ChanSaveName = [ChanSaveName, '_ncells', num2str(p.ncells)]; end
    if ~isempty(p.subset); for iSubset = 1:numel(p.subset); ChanSaveName = [ChanSaveName, ' ', p.subset{iSubset}]; end; end
    if ~isempty(p.exclude); ChanSaveName = [ChanSaveName, ' ', 'excluding']; for iSubset = 1:numel(p.exclude); ChanSaveName = [ChanSaveName, ' ', p.exclude{iSubset}]; end; end
    if p.combinexys; ChanSaveName = [ChanSaveName, ' combinedxys']; end
    if numel(facetby) > 1; ChanSaveName = [ChanSaveName,' ',facetby{s}]; end
    if ~isempty(p.overlaptx); ChanSaveName = [ChanSaveName, ' overlaptx ', num2str(p.overlaptx)]; end
    ChanSaveName = [ChanSaveName, ' afterTX ', num2str(p.aftertreatment)];
    ChanSaveName = [ChanSaveName, ' ' , plotType];
    ChanSaveName = strrep(ChanSaveName, '.',' pt ');
    
    FullSaveName = [dataloc.fold.fig,'\',ChanSaveName];

    % Give the figure a title
    tittxt = [dataloc.file.base,' ', ChanSaveName];
    tittxt = strrep(tittxt,'_',' ');
    sgtitle(tittxt,'FontSize', p.font_size);
    fontsize(gcf,p.font_size,"points"); fontname(gcf,'Calibri');
    if ~isempty(findall(gcf,'Type','Axes')) %dont save an empty figuer 
        saveas(gcf, FullSaveName, 'svg')
    end
    if exist('barHold','var') && any("pValue" == string(barHold.Properties.VariableNames)) %print the table for the data
        writetable(barHold,[FullSaveName,'_data'],'FileType','spreadsheet')
        clear barHold
    end
    if p.closefigs; close(gcf); end
    end %plot type loop
    end %check for empty plots 
    end %facetby loop
end %iplot loop
    
    %save  arrayfun(@(x)cellfun(@cell2mat,{dataloc.z{x}.(ThisChan).pkpos},'UniformOutput',false),XyNums,'UniformOutput',false)
end %for PulsePlotter Function
%% Combine XYs Function
function [dataloc,idx2,linetp,xymat,goodxy,p] = Combine_XYs(dataloc,txs,cellfn,idx,xymat,linetp,goodxy,p)
fns = fieldnames(idx);
for ifn = 1:numel(fns)
    idx2.(fns{ifn}) = []; 
end
linetp2 = {};

%figure out how many xys there needs to be
totalXYs = 0;
for iCellline = 1:numel(cellfn)
    for iTxs = 1:numel(txs)
        XYs2Comb = intersect(idx.(txs{iTxs}),idx.(cellfn{iCellline}));
        if XYs2Comb>0; totalXYs = totalXYs + 1; end
    end
end

PlotD = struct('d',[],'z',[],'IFd',[],'txx',[]);
PlotD.d = cell([1,totalXYs]); PlotD.z = cell([1,totalXYs]); PlotD.IFd = cell([1,totalXYs]);

newXYNums = 1:totalXYs;
xyCounter = 0;
sizeD = size(dataloc.d,2); % what is the max d number

%find matching data to combine
for iCellline = 1:numel(cellfn) %cell line loop
    for iTxs = 1:size(txs,1)
        XYs2Comb = intersect(idx.(txs{iTxs}),idx.(cellfn{iCellline}));
        XYs2Comb = XYs2Comb((XYs2Comb <= sizeD)); % check that the xys requested are not outside the bounds of d
        goodxys = arrayfun(@(x)~isempty(dataloc.d{x}),XYs2Comb);
        XYs2Comb = XYs2Comb(goodxys);
        numXYs = length(XYs2Comb);
        if ~isempty(XYs2Comb)
        fieldNames = fieldnames(dataloc.d{XYs2Comb(1)}.data);% get the channels in the xys to combine
        if any(contains(fieldNames, 'Corr')); fieldNames = fieldNames(~contains(fieldNames, 'Corr')); end
        NumChans = length(fieldNames);
        firstxy = true;
        for iXY = 1:numXYs
            ThisXY = XYs2Comb(iXY);
            %xyIndex = cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0));
            if goodxy(ThisXY) && ThisXY <= length(dataloc.d)
                if firstxy
                    xyCounter = xyCounter + 1; newXYNum = newXYNums(xyCounter); idx2.(txs{iTxs}) = [idx2.(txs{iTxs}), newXYNum]; 
                end
                if ~isempty(dataloc.d{ThisXY}) %check the XY isn't empty
                if isfield(dataloc.d{ThisXY}, 'data') %check the XY has data
                if ~isfield(PlotD.d{newXYNum},'data'); PlotD.d{newXYNum}.data = struct(); end
                if ~isfield(PlotD.z{newXYNum},'data'); PlotD.z{newXYNum}.data = struct(); end
                if ~isfield(PlotD.IFd{newXYNum},'data'); PlotD.IFd{newXYNum}.data = struct(); end

                for iChan = 1:NumChans
                    ThisChan = fieldNames{iChan};        % get the channel name
                    if isfield(dataloc.d{ThisXY}.data, ThisChan) %check if that XY has the channel in it
                        HasZ = false;
                        if isfield(dataloc, 'z') && ~isempty(dataloc.z) && ThisXY <= length(dataloc.d)
                            if isfield(dataloc.z{ThisXY}, 'data') && ~isempty(dataloc.z{ThisXY}.data)
                                HasZ = isfield(dataloc.z{ThisXY}.data, ThisChan); %check if that XY has the z data in it
                            end
                        end %check for z in dataloc
                        
                        HasIFd = false;
                        if isfield(dataloc, 'IFd')
                            if ~isempty(dataloc.IFd)
                                HasIFd = isfield(dataloc.IFd{ThisXY}.data, p.ifchan); %check if that XY has the IF data in it
                            end
                        end %check for IFd in dataloc
                        
                    if ~isfield(PlotD.d{newXYNum}.data, ThisChan)
                        PlotD.d{newXYNum}.data.(ThisChan) = dataloc.d{ThisXY}.data.(ThisChan); %hold the data
                        if firstxy
                            linetp2{newXYNum} = linetp{cell2mat(cellfun(@(x)any(x == ThisXY),xymat,'un',0))}; % get the txx tps for xlines
                            firstxy = false;
                        end
                        
                        if HasZ
                            PlotD.z{newXYNum}.data.(ThisChan) = dataloc.z{ThisXY}.data.(ThisChan);
                        end %has Z data
                        
                        if (HasIFd && iChan < 2)
                            PlotD.IFd{newXYNum}.data.(p.ifchan) = dataloc.IFd{ThisXY}.data.(p.ifchan);
                        end %has IF data
                        
                    else %add the xy to the other data that already exists
                        PlotD.d{newXYNum}.data.(ThisChan) = [PlotD.d{newXYNum}.data.(ThisChan); dataloc.d{ThisXY}.data.(ThisChan)];
                        
                        if HasZ
                            PlotD.z{newXYNum}.data.(ThisChan) = [PlotD.z{newXYNum}.data.(ThisChan); dataloc.z{ThisXY}.data.(ThisChan)];
                        end %has Z data
                        
                        if (HasIFd && iChan < 2)
                            PlotD.IFd{newXYNum}.data.(p.ifchan) = [PlotD.IFd{newXYNum}.data.(p.ifchan); dataloc.IFd{ThisXY}.data.(p.ifchan)];
                        end %has IF data
                        
                    end %channel is field check
                    end %xy having channel check   
                end %channel loop
                end %xy has data field
                end %check if xy is empty
            end %good XY check
        end %each xy loop
        end %xys to comb
    if  numel(XYs2Comb) > 0 && exist('newXYNum', 'var'); idx2.(cellfn{iCellline}) = [idx2.(cellfn{iCellline}), newXYNum]; end 
    end %treatment loop  
    %now add the treatments that celline got to it's treatment list
    
end %cell line loop

% Now replace dataloc and tx with the new stuff

dataloc.d = cell(size(PlotD.d,2)); dataloc.d = PlotD.d;
dataloc.IFd = cell(size(PlotD.IFd,2)); dataloc.IFd = PlotD.IFd;
dataloc.z = cell(size(PlotD.z,2)); dataloc.z = PlotD.z;
goodxy = ~cellfun(@isempty, dataloc.d)'; 
linetp = linetp2';

xymat = cell(numel(newXYNums),1);
for ixymat = 1:numel(newXYNums)
    xymat{ixymat} = newXYNums(ixymat);
end
end %end Combine_XYs function

%% Path re-adder function
function Path_Adder()
    CurrPaths = path;

    % first one is the path, second is if it needs genpath too (1 is yes)
    PathsToCheck = {...
    '\\albecklab.mcb.ucdavis.edu\data\Code\DatalocHandler', 0;...
    '\\albecklab.mcb.ucdavis.edu\data\Code\Cell Trace', 0;...
    '\\albecklab.mcb.ucdavis.edu\data\Code\Image Analysis', 0;...
    '\\albecklab.mcb.ucdavis.edu\data\Code\PLS', 0};

    pathCell = regexp(CurrPaths, pathsep, 'split');
    pathCell = cellfun(@lower, pathCell, 'UniformOutput',false);
    % FIX
%     [~,b] = system('net use');
%     tScan = textscan(b,'%s'); tScan = tScan{1};
%     hasLetters = cellfun(@(x)regexpi(x,'\S:'),tScan,'UniformOutput',false);
%     nothingHere = cellfun(@isempty,hasLetters);
%     hasLetters(nothingHere)=deal({0}); hasLetters = logical(cell2mat(hasLetters));
%     theseLetters = tScan(hasLetters);
%     theseLetters(:,2) = tScan(find(hasLetters)+1);
%     itHasThis = cell2mat(cellfun(@(x)contains(pathCell,x),theseLetters(:,1),'UniformOutput',false));
%     if any(itHasThis)
%         replaceThisLetter = find(itHasThis);
%         icalledyou = strrep(icalledyou,theseLetters(replaceThisLetter,1),theseLetters{replaceThisLetter,2});
%     end
%     if iscell(icalledyou); icalledyou = icalledyou{1}; end

    NotOnPath = ~cell2mat(cellfun(@(x)ismember(x, pathCell), lower(PathsToCheck(:,1)), 'UniformOutput',false));
    AddThese = PathsToCheck(NotOnPath,:);
    if ~isempty(AddThese)
        for iAdd = 1:size(AddThese,1)
            if AddThese{iAdd,2} == 1
                addpath(genpath(AddThese{iAdd,1}))
            else
                addpath(AddThese{iAdd,1})
            end
        end
    end
    %% Now check if the subfolders have been added for each folder with subfolders
    PathsToCheck = PathsToCheck(1:3,1); % fix to be a cell fun
    for iAdd = 1:size(PathsToCheck,1)
        CurFold = lower(PathsToCheck{iAdd});
        CurDir = dir(CurFold);
        DirFlags = [CurDir.isdir];
        SubFolders = {CurDir(DirFlags).name};
        SubFolders2 = cellfun(@(x)strcat([CurFold,'\'],x), lower(SubFolders), 'UniformOutput', false);
        SubFolders = SubFolders2(3:numel(SubFolders2));
        NotOnPath2 = ~cell2mat(cellfun(@(x)ismember(x, pathCell), SubFolders, 'UniformOutput', false));
        AddThese2 = SubFolders(NotOnPath2);
        if ~isempty(AddThese2)
            for iAdd2 = 1:size(AddThese2,2)
                addpath((AddThese2{1,iAdd2}))
            end
        end
    end

    %% Now remove the evil utrack histogram thing
    hasEvilUtrack = contains(pathCell,"u-track");
    if any(hasEvilUtrack)
        PathsToRemove = pathCell(hasEvilUtrack);
        cellfun(@(x)rmpath(x),PathsToRemove)
    end


end
%% findZero
% Makes the plot offset to zero for a treatment
function findZero(txs,thisTX,tracklength,firsttp,p)
    %fix the axes so time 0 is the given tx
    %if any(arrayfun(@strlength,xlabs)>3); arrayfun(@(x) x(1:3),xlabs(arrayfun(@strlength,xlabs)> 3)) %keep the x labels a reasonable length
    if ((tracklength-firsttp) > (12*p.tktm)); ticklengthz = 2; % if the tracks are more than 12 hrs do a tick every other hour
    elseif ((tracklength-firsttp) > (4*p.tktm)); ticklengthz = 1; % if the tracks are more than 4 hrs do a tick every hour
    else; ticklengthz = 0.5; % if the tracks are less than 4 hrs do a tick every half hour
    end
    if isempty(p.zerohrtx) % see if you want a different zero hour tx if not use aftertreatment
        zeroHrTx = thisTX;
    else
        if numel(txs) < p.zerohrtx; zeroHrTx = txs(end); 
        else; zeroHrTx=txs(p.zerohrtx); end % get treatment number or last treatment
    end
    if ~isempty(zeroHrTx) && ~p.plotfromzero %offset all data to the given treatment
        xBack = zeroHrTx:-p.tktm*ticklengthz:0; xBack = sort(xBack(2:end));
        xtix = [xBack, zeroHrTx:p.tktm*ticklengthz:tracklength]; %make a tick for every given interval
        xlabs = (xtix-zeroHrTx)/p.tktm; % off set the time and make the labels
    
    elseif ~isempty(zeroHrTx) && p.plotfromzero %offset all data to the given treatment
        xtix = [1:p.tktm*ticklengthz:(tracklength-zeroHrTx)]; %make a tick for every given interval
        xlabs = (xtix-1)/p.tktm; % off set the time and make the labels                        
    else
        xtix = [0:p.tktm*ticklengthz:tracklength]; %make a tick for every given interval
        xlabs = xtix/p.tktm;
    end
    xticks(xtix); xlabs = string(xlabs);
    xticklabels(xlabs);
    xtickangle(0);
    xlabel('');

end % findZero
%% Spatialheatmap function
function ah = spatialheatmap(ah,dataloc,iwell,channel,color_map,p,firsttp,tracklength)
%%create a spatial heatmap of data from iwell in d where d is the output of
%%ct_dataproc
%%% This works best if you have full length tracks and can see spatial
%%% activtiy in the movie.

%iwell = well that you want to visualize
%channel = which channel you want to visualize ex: 'ampkar'
%OPTIONAL INPUTS
%color_map = default is parula. Dont input as string
%xname and yname = if you x and y coord is named something other than
%XCoord and YCoord, inupt them here as strings

% if you want to change heatmap xticklabels and etc. use sh2. ex:
%xlabel(sh2, 'time')

if p.addif
    a1 = subplot(1,2,1,'Parent',gca);
else
    axes(ah)
end
if ~exist('colormap','var')
    %  parameter does not exist, so default it to something
    color_map = parula;
end

%get time average x and y coordinates and put them into a matrix
averagex = nanmean(dataloc.d{iwell}.data.XCoord(:,firsttp:tracklength),2);
averagey = nanmean(dataloc.d{iwell}.data.YCoord(:,firsttp:tracklength),2);
coords = [averagex, averagey];

keepData = ~any(isnan(coords),2); % find the non-nan data
coords = coords(keepData,:);  % keep real data

%% get distance of each cell from others and perform hierarchical clustering based on distance
DIST = pdist2(coords,coords);
LINK = linkage(DIST,'average');
leafOrder = optimalleaforder(LINK,DIST);
if sum(keepData) < 3 || numel(leafOrder) < 3 ; return; end
dHold = dataloc.d{iwell}.data.(channel)(keepData,:);
dHold = dHold(leafOrder',firsttp:tracklength);
if ~isempty(p.ncells) && (size(dHold,1) >= p.ncells)
    dHold = dHold(1:p.ncells,:);
end
imagesc(dHold)
if ~isempty(p.ymn) && ~isempty(p.ymx); clim([p.ymn, p.ymx]); end
colormap(color_map)
if p.addif
    a2 = subplot(1,2,2,'Parent',ah);
    ifHold = dataloc.IFd{iwell}.data.(p.ifchan)(keepData,1);
    ifHold = ifHold(leafOrder',1);
    if ~isempty(p.ncells) && (size(dHold,1) >= p.ncells)
        ifHold = ifHold(1:p.ncells,:);
    end
    imagesc(ifHold)
    clear ifHold;
end
clear dHold;
yticks([])
end
