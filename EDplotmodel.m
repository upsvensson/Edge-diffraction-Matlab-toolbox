function EDplotmodel(eddatafile,varargin)
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
%
% Sources and Rdata.receivers are taken from an sdatafile and an rdatafile, the file name of
% which is assumed to be similar to the eddatafile.
%
% Uses functions EDstrpend, EDstrpblnk
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 26 Jan 2018

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

p = inputParser;

p.addOptional('plotoptions',0,@isscalar);
p.addParamValue('sounumbers',[],@isnumeric);
p.addParamValue('recnumbers',[],@isnumeric);
p.addParamValue('edgenumbers',[],@isnumeric);

p.parse(varargin{:})
inputs = p.Results;

plotoptions = inputs.plotoptions;
sounumbers = inputs.sounumbers;
recnumbers = inputs.recnumbers;
edgenumbers = inputs.edgenumbers;
 
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

% if nargin == 1
%     if isstr(eddatafile) ~= 1
%         plotoptions = eddatafile;
%     else
%         plotoptions = 0;    
%     end
% elseif nargin == 0
%     plotoptions = 0;    
% end
% 
% if nargin < 3
%     plotoptions2 = [];
% else
%     nplotoptions2 = size(plotoptions2,1);
% end

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
    eval(['load ',eddatafilepath,filesep,eddatafile])
else
    eval(['load ',eddatafile])    
end
    
ncornersperplanevec = double(planedata.ncornersperplanevec);
if plotsources
    sdatafile = EDstrpblnk([eddatafilepath,Filestem,'_Sdata.mat']);
    if exist(sdatafile) == 2
        eval(['load ',sdatafile])
    else
        error(['ERROR: The sdata file named ',sdatafile,' could not be opened'])    
    end
end
if plotreceivers
    rdatafile = EDstrpblnk([eddatafilepath,Filestem,'_Rdata.mat']);
    if exist(rdatafile) == 2
        eval(['load ',rdatafile])
    else
        error(['ERROR: The rdata file named ',rdatafile,' could not be opened'])    
    end
end

if plotplCADnumbers || plotcoCADnumbers
    cadgeofile = EDstrpblnk([eddatafilepath,Filestem,'_cadgeo.mat']);
    if exist(cadgeofile) == 2
        eval(['load ',cadgeofile])
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

% figure(1)
hold off
for ii = 1:nedges
	co1 = edgedata.edgecorners(ii,1);
	co2 = edgedata.edgecorners(ii,2);
	iv = [co1;co2];
	h = plot3(planedata.corners(iv,1),planedata.corners(iv,2),planedata.corners(iv,3),linemarkertype(linemarker(ii),:));
    set(h,'LineWidth',linewidthvalues(linemarker(ii)));
    if ii ==1
        view(viewpos)
        hold on
    end
    if plotednumbers
        costart = planedata.corners(co1,:) + 0.1*(planedata.corners(co2,:)-planedata.corners(co1,:));
        h = plot3(costart(1),costart(2),costart(3),'ko');    
    end
end

if plotnvecs
	for ii = 1:nplanes
        midpoint = mean(planedata.corners(planedata.planecorners(ii,1:(planedata.ncornersperplanevec(ii))),:));
        
        endpoint = midpoint + planedata.planeeqs(ii,1:3);
        bothpoints = [midpoint;endpoint];
        plot3(bothpoints(:,1),bothpoints(:,2),bothpoints(:,3));
        plot3(midpoint(1),midpoint(2),midpoint(3),'ro');
	end
end

if plotplnumbers
	for ii = 1:nplanes
        midpoint = mean(planedata.corners(planedata.planecorners(ii,1:planedata.ncornersperplanevec(ii)),:));
        endpoint = midpoint + planedata.planeeqs(ii,1:3)*0.1;
        if plotplCADnumbers == 1
            h = text(endpoint(1),endpoint(2),endpoint(3),['p',int2str(planenumbers(ii))]);
            set(h,'FontSize',14);
        else
            h = text(endpoint(1),endpoint(2),endpoint(3),['p',int2str(ii)]);
            set(h,'FontSize',14);
        end
    
	end
end

if plotednumbers
    for ii = 1:nedges
        midpoint = mean(planedata.corners(edgedata.edgecorners(ii,1:2),:));
        h = text(midpoint(1),midpoint(2),midpoint(3),['e',int2str(ii)]);
            set(h,'FontSize',18);
    end
end

if plotconumbers
    for ii = 1:ncorners
        if plotcoCADnumbers == 1
            h = text(planedata.corners(ii,1),planedata.corners(ii,2),planedata.corners(ii,3),['c',int2str(cornernumbers(ii))]); 
            set(h,'FontSize',18);
        else
            h = text(planedata.corners(ii,1),planedata.corners(ii,2),planedata.corners(ii,3),['c',int2str(ii)]);    
            set(h,'FontSize',18);
        end
    end
    
end

if plotsources == 1
    if isempty(sounumbers)
        plot3(Sdata.sources(:,1),Sdata.sources(:,2),Sdata.sources(:,3),'*');
    else
        plot3(sources(sounumbers,1),sources(sounumbers,2),sources(sounumbers,3),'*')
    end
end

if plotreceivers == 1
    if isempty(recnumbers)
        plot3(Rdata.receivers(:,1),Rdata.receivers(:,2),Rdata.receivers(:,3),'ro');
    else
        plot3(Rdata.receivers(recnumbers,1),Rdata.receivers(recnumbers,2),Rdata.receivers(recnumbers,3),'ro');
    end
end

xlabel('x')
ylabel('y')
zlabel('z')
axis equal
rotate3d on


