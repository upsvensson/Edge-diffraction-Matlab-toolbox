function plothandles = EDplotmodel(eddatafile,varargin);
% EDplotmodel - Plots a model which is given in an eddatafile.
%
% Input parameters:
%   eddatafile (optional)   If an input file is not specified, a file
%                           opening window will be presented.
%   plotoptions (optional)
%                           The text-string 'plotoptions',with an integer can give extra options:
%                           if bit0 = 1 (= 1) => plot sources
%                           if bit1 = 1 (= 2) => plot Rdata.receivers
%                           if bit2 = 1 (= 4) => plot plane normal vectors
%                           if bit3 = 1 (= 8) => print plane numbers
%                           if bit4 = 1 (=16) => print edge numbers (and
%                           indicate start end each edge with a little
%                           circle)
%                           if bit5 = 1 (=32) => print corner numbers
%                           if bit6 = 1 (=64) => print plane numbers using the
%                                                CAD file numbering
%                           if bit7 = 1 (=128)=> print corner numbers using the
%                                                CAD file numbering
%                           Example: the integer 11 = 1011 binary,
%                           so bits 0,1, and 3 are set.
%   'sounumbers',vector of source numbers (optional) 
%                           If the text-string 'sounumbers' is given,
%                           followed by a vector of integers, only those
%                           source numbers are plotted.
%   'recnumbers',vector of source numbers (optional) 
%                           If the text-string 'recnumbers' is given,
%                           followed by a vector of integers, only those
%                           source numbers are plotted.
%   'edgenumbers',vector of edge numbers (optional)
%                           If the text-string 'edgenumbers' is given,
%                           followed by a vector of integers, only those
%                           edge numbers are plotted with a solid line.
%                           The others are plotted with dashed lines. 
%   'figurewindow' a desired figure number. If not specified, figure number
%   1 will be used.
% 
% Output parameters:
%   plothandles     A struct with fields that give handles to the various
%                   parts of the figure:
%                   .parent   To set the size of the axis values, change
%                             the parameter 'FontSize' of this figure
%                             handle etc.
%                   .sources
%                   .receivers
%                   .edges        A list, one value for each edge
%                   .edgecircles  A list, one value for each edge
%                   .edgenumbers  A list, one value for each edge
%                   .conumbers    A list, one value for each corner
%                   .nvecs        A list, one value for each nvec
%                   .nveccircles  A list, one value for each nvec
%                   .planenumbers A list, one value for each plane
% 
% Sources and Rdata.receivers are taken from an sdatafile and an rdatafile, the file name of
% which is assumed to be similar to the eddatafile.
%
% Uses functions EDstrpend, EDstrpblnk
%
% Peter Svensson (peter.svensson@ntnu.no) 14 March 2021
% 
% plothandles = EDplotmodel(eddatafile,varargin);

% 29 April 2016 Completed the use of structs planedata and edgedata
% 1 Nov 2017 Added EDstrpblnk for the sdata etc file names. Also skipped the
% question about view. Always chooses [1 1 1].
% 17 Nov. 2017 Added a little marker at the start of each edge.
% 17 Nov. 2017 Implemented the inputParser functionality
% 30 Nov. 2017 Removed a remaining "plotoptions2"
% 18 Jan 2018 Copied from ESIE2toolbox
% 26 Jan 2018 Small fix: Sdata and Rdata filenames were wrong (but still
% worked??). Also removed the little circles on edges if edgenumbers are
% not plotted.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 14 March 2021 Added plothandles as output parameter, and figurewindow as
% optional input parameter

p = inputParser;

p.addOptional('plotoptions',0,@isscalar);
p.addParamValue('sounumbers',[],@isnumeric);
p.addParamValue('recnumbers',[],@isnumeric);
p.addParamValue('edgenumbers',[],@isnumeric);
p.addParamValue('figurewindow',1,@isscalar);

p.parse(varargin{:})
inputs = p.Results;

plotoptions = inputs.plotoptions;
sounumbers = inputs.sounumbers;
recnumbers = inputs.recnumbers;
edgenumbers = inputs.edgenumbers;
figurewindow = inputs.figurewindow;

%--------------------------------------------------------------
% Output struct

plothandles = struct('parent',[],...
    'sources',[],'receivers',[],'edges',[],'edgecircles',[],...
    'edgenumbers',[],'conumbers',[],'nvecs',[],'nveccircles',[],...
    'planenumbers',[]);
                     
%--------------------------------------------------------------
% Read the eddatafile

if nargin == 0 || (nargin >= 1 && isstr(eddatafile) ~= 1) || (nargin >=1 && isstr(eddatafile) && length(eddatafile)==0)
	[eddatafile,eddatafilepath] = uigetfile('*eddata.mat','Please select the eddatafile');
    [eddatafilepath,~,~] = fileparts(eddatafilepath);
	if ~isstr(eddatafile) || isempty(eddatafile)
		return
	end
else
	[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
    eddatafile = [eddatafile,fileext];
end
Filestem = EDstrpend(eddatafile,'_eddata');

disp(['eddatafile is ',eddatafilepath,filesep,eddatafile])

%--------------------------------------------------------------
% Read and decode the second input parameter

plotsources = bitget(plotoptions,1);
plotreceivers = bitget(plotoptions,2);
plotnvecs = bitget(plotoptions,3);
plotplnumbers = bitget(plotoptions,4);
plotednumbers = bitget(plotoptions,5);
plotconumbers = bitget(plotoptions,6);
plotplCADnumbers = bitget(plotoptions,7);
plotcoCADnumbers = bitget(plotoptions,8);

if plotplCADnumbers, plotplnumbers = 1; end
if plotcoCADnumbers, plotconumbers = 1; end

%--------------------------------------------------------------------------
% Load the needed input files

if ~isempty(eddatafilepath)
    eval(['load ''',eddatafilepath,filesep,eddatafile,''''])
else
    eval(['load ''',eddatafile,''''])    
end
    
ncornersperplanevec = double(planedata.ncornersperplanevec);
if plotsources
    sdatafile = EDstrpblnk([eddatafilepath,Filestem,'_Sdata.mat']);
    if exist(sdatafile) == 2
        eval(['load ''',sdatafile,''''])
    else
        error(['ERROR: The sdata file named ',sdatafile,' could not be opened'])    
    end
end
if plotreceivers
    rdatafile = EDstrpblnk([eddatafilepath,Filestem,'_Rdata.mat']);
    if exist(rdatafile) == 2
        eval(['load ''',rdatafile,''''])
    else
        error(['ERROR: The rdata file named ',rdatafile,' could not be opened'])    
    end
end

if plotplCADnumbers || plotcoCADnumbers
    cadgeofile = EDstrpblnk([eddatafilepath,Filestem,'_cadgeo.mat']);
    if exist(cadgeofile) == 2
        eval(['load ''',cadgeofile,''''])
    else
        error(['ERROR: The cadgeo file named ',cadgeofile,' could not be opened'])            
    end
end
%--------------------------------------------------------------

ncorners = size(planedata.corners,1);
nedges = size(edgedata.edgecorners,1);
nplanes = size(planedata.planeeqs,1);

planelist = 1:nplanes;

viewpos = [1 1 1];

%--------------------------------------------------------------
% Plot the edges. Check if some should have dashed lines

if isempty(edgenumbers)
    linemarker = ones(nedges,1);
else
    linemarker = 2*ones(nedges,1);    
    linemarker(edgenumbers) = 1;
end
linemarkertype = ['k- ';'k--'];
linewidthvalues = [2 1];

figure(figurewindow)
clf(figurewindow)
for ii = 1:nedges
	co1 = edgedata.edgecorners(ii,1);
	co2 = edgedata.edgecorners(ii,2);
	iv = [co1;co2];
	plothandles.edges(ii) = plot3(planedata.corners(iv,1),planedata.corners(iv,2),planedata.corners(iv,3),linemarkertype(linemarker(ii),:));
    set(plothandles.edges(ii),'LineWidth',linewidthvalues(linemarker(ii)));
    if ii ==1
        view(viewpos)
        hold on
        plothandles.parent = get(plothandles.edges(1),'Parent');
    end
    if plotednumbers
        costart = planedata.corners(co1,:) + 0.1*(planedata.corners(co2,:)-planedata.corners(co1,:));
        plothandles.edgecircles(ii) = plot3(costart(1),costart(2),costart(3),'ko');    
    end
end

if plotnvecs
	for ii = 1:nplanes
        midpoint = mean(planedata.corners(planedata.planecorners(ii,1:(planedata.ncornersperplanevec(ii))),:));
        
        endpoint = midpoint + planedata.planeeqs(ii,1:3);
        bothpoints = [midpoint;endpoint];
        plothandles.nvecs(ii) = plot3(bothpoints(:,1),bothpoints(:,2),bothpoints(:,3));
        plothandles.nveccircles(ii) = plot3(midpoint(1),midpoint(2),midpoint(3),'ro');
	end
end

if plotplnumbers
	for ii = 1:nplanes
        midpoint = mean(planedata.corners(planedata.planecorners(ii,1:planedata.ncornersperplanevec(ii)),:));
        endpoint = midpoint + planedata.planeeqs(ii,1:3)*0.1;
        if plotplCADnumbers == 1
            plothandles.planenumbers(ii) = text(endpoint(1),endpoint(2),endpoint(3),['p',int2str(planenumbers(ii))]);
            set(plothandles.planenumbers(ii),'FontSize',14);
        else
            plothandles.planenumbers(ii) = text(endpoint(1),endpoint(2),endpoint(3),['p',int2str(ii)]);
            set(plothandles.planenumbers(ii),'FontSize',14);
        end
    
	end
end

if plotednumbers
    for ii = 1:nedges
        midpoint = mean(planedata.corners(edgedata.edgecorners(ii,1:2),:));
        plothandles.edgenumbers(ii) = text(midpoint(1),midpoint(2),midpoint(3),['e',int2str(ii)]);
            set(plothandles.edgenumbers(ii),'FontSize',18);
    end
end

if plotconumbers
    for ii = 1:ncorners
        if plotcoCADnumbers == 1
            plothandles.conumbers(ii) = text(planedata.corners(ii,1),planedata.corners(ii,2),planedata.corners(ii,3),['c',int2str(cornernumbers(ii))]); 
            set(plothandles.conumbers(ii),'FontSize',18);
        else
            plothandles.conumbers(ii) = text(planedata.corners(ii,1),planedata.corners(ii,2),planedata.corners(ii,3),['c',int2str(ii)]);    
            set(plothandles.conumbers(ii),'FontSize',18);
        end
    end
    
end

if plotsources == 1
    if isempty(sounumbers)
        plothandles.sources = plot3(Sdata.sources(:,1),Sdata.sources(:,2),Sdata.sources(:,3),'*');
    else
        plothandles.sources = plot3(sources(sounumbers,1),sources(sounumbers,2),sources(sounumbers,3),'*');       
    end
    set(plothandles.sources,'LineWidth',3)
end

if plotreceivers == 1
    if isempty(recnumbers)
        plothandles.receivers = plot3(Rdata.receivers(:,1),Rdata.receivers(:,2),Rdata.receivers(:,3),'ro');
    else
        plothandles.receivers = plot3(Rdata.receivers(recnumbers,1),Rdata.receivers(recnumbers,2),Rdata.receivers(recnumbers,3),'ro');
    end
    set(plothandles.receivers,'LineWidth',3)
end

xlabel('x')
ylabel('y')
zlabel('z')
axis equal
rotate3d on


