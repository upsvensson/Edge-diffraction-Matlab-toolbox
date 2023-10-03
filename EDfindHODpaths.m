function [hodpaths,hodpathsalongplane,outpar,existingfilename] = EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
difforder,EDversionnumber,filehandlingparameters)
% EDfindHODpaths finds the higher order diffraction paths for a convex
% object, or a non-convex object where specular reflections are ignored.
% 
% A version 1 of this function was used up to v 0.221 of EDtoolbox
% and version 2 after that.
%
% Input parameters:
%   edgeseesedge        Matrix, [nedges,nedges], of -1, 0 and 1
%   visedgesfroms       Matrix, [nedges,nsources], of 0 and 1
%   visedgesfromr       Matrix, [nedges,nreceivers], of 0 and 1
%   difforder
%   EDversionnumber
%   filehandlingparameters  Ignored in version 1 of this function and 
%                       obligatory for version 2.
%
% Output parameters:
%   hodpaths            Cell variable, hodpaths{2} has size [npaths,2] etc
%                       Each row in each matrix has a sequence of edge
%                       numbers.
%   hodpathsalongplane  Cell variable, hodpaths{2} has size [npaths,1] etc
%                       Each row in each matrix has values 0 or 1, which tells
%                       for each hod sequence "leg" whether it passes along
%                       a plane or not.
%   outpar                  This parameter is different for version 1 and
%                           version 2 of this function. 
%       v1 outpar = EDinputdatahash       
%                           This is a string of characters which is
%                           uniquely representing the input data. 
%                           An existing result file with the same value of
%                           this EDinputdatahash will be reused.
%       v2 outpar = elapsedtimehodpaths
%                           This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimefindHODpaths has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        For v2 of this function, if an existing file was
%                           found that was reused, then the reused file
%                           name is given here. If no existing file could be 
%                           reused then this variable is empty. For v1 of
%                           this function, this variable is also empty.
% 
% Uses the function EDrecycleresultfiles from EDtoolbox
% Uses the function Datahash from Matlab Central
% 
% Peter Svensson (peter.svensson@ntnu.no) 29 Sep. 2023
% 
% [hodpaths,hodpathsalongplane,outpar,existingfilename] = EDfindHODpaths...
% (edgeseesedge,visedgesfroms,visedgesfromr,...
% difforder,EDversionnumber,filehandlingparameters)

% 15 Mar 2018 First version
% 16 Mar 2018 Expanded to several sources and receivers
% 21 May 2019 Made some changes so that negative values of edgeseesedge can
% be handled. Introduced the output variable hodpathsalongplane.
% 1 June 2019 Fixed bug; hadn't finished the handling of negative values of
% edgeseesedge.
% 29 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.

t00 = clock;

if nargin < 6  % Must be the old version
	functionversion = 1;
else   % nargin = 6 -> must be the new version
    functionversion = 2;
end

EDinputdatastruct = struct('difforder',difforder,...
        'edgeseespartialedge',edgeseesedge,...
        'vispartedgesfroms',visedgesfroms,'vispartedgesfromr',visedgesfromr,...
        'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

if functionversion == 1
	outpar = EDinputdatahash;
	existingfilename = '';
elseif functionversion == 2
	%---------------------------------------------------------------
	% Sort out the file business: can an existing file be used?
	% Then copy the existing file to a new copy. Should the data be saved in a file? 

	if filehandlingparameters.suppressresultrecycling == 1
		foundmatch = 0;
		existingfilename = '';
	else
		[foundmatch,existingfilename] = ... 
			EDrecycleresultfiles(filehandlingparameters.outputdirectory,...
			'_hodpaths',EDinputdatahash);
	end

	desiredname = [filehandlingparameters.outputdirectory,filesep,...
		filehandlingparameters.filestem,'_hodpaths.mat'];

	if foundmatch == 1
		eval(['load(''',existingfilename,''')'])
		if ~strcmp(existingfilename,desiredname)
			copyfile(existingfilename,desiredname);
		end
		elapsedtimehodpaths_new = etime(clock,t00);
		elapsedtimehodpaths = [elapsedtimehodpaths_new elapsedtimehodpaths];
        outpar = elapsedtimehodpaths;
		return
	end
end

%----------------------------------------------------------------
% No existing file can be used

nedges = size(edgeseesedge,2);
nsources = size(visedgesfroms,2);
nreceivers = size(visedgesfromr,2);

hodpaths = cell(difforder,nreceivers,nsources);
hodpathsalongplane = cell(difforder,nreceivers,nsources);

classtocheck = class(edgeseesedge);
if strcmp(classtocheck,'int8')
    visedgesfromr = int8(visedgesfromr);
else
    if strcmp(classtocheck,'int16')
        visedgesfromr = int16(visedgesfromr);
    end
end
edgeseesedge_lastround = zeros(nedges,nedges,nreceivers);
for ii = 1:nreceivers
    edgeseesedge_lastround(:,:,ii) = edgeseesedge.*visedgesfromr(:,ones(1,nedges)*ii);        
end

for isou = 1:nsources
    startedges = find(visedgesfroms(:,isou));

    pathspattern_to_propagate = startedges;
    
    for norder = 2:difforder

        npaths = size(pathspattern_to_propagate,1);

        for ii = 1:nreceivers
            edgeseesedgevalues = edgeseesedge_lastround(:,pathspattern_to_propagate(:,norder-1),ii); 
            ivec = find(edgeseesedgevalues);
            [paths1,paths2] = ind2sub([nedges,npaths],ivec);
            hodpaths{norder,ii,isou}           = [pathspattern_to_propagate(paths2,:) paths1];
            indexlist = sub2ind([nedges nedges],pathspattern_to_propagate(paths2,end),paths1);
            hodpathsalongplane_coltoadd = (edgeseesedge(indexlist)+1)/2;
            if norder > 2
                hodpathsalongplane{norder,ii,isou} = [pathsalong_to_propagate(paths2,:) hodpathsalongplane_coltoadd];
            else
                hodpathsalongplane{norder,ii,isou} = hodpathsalongplane_coltoadd;
            end
            
%             negivec = find(edgeseesedgevalues(ivec) == -1)
%             hodpathsalongplane{norder,ii,isou}(negivec,norder-1) = 0; 
%             if norder >= 3
%                 
%                 hodpathsalongplane{norder,ii,isou}(:,1:norder-2) = ...
%                     hodpathsalongplane{norder-1,ii,isou}(paths2,:);
%             end
        end

        ivec = find(edgeseesedge(:,pathspattern_to_propagate(:,norder-1)));
        [paths1,paths2] = ind2sub([nedges,npaths],ivec);
        pathspattern_to_propagate = [pathspattern_to_propagate(paths2,:) paths1];
        if norder == 2
           indexlist = sub2ind([nedges nedges],pathspattern_to_propagate(:,1),pathspattern_to_propagate(:,2)) ;
           pathsalong_to_propagate = (edgeseesedge(indexlist)+1)/2;
        else            
            indexlist = sub2ind([nedges nedges],pathspattern_to_propagate(:,end-1),pathspattern_to_propagate(:,end));
            pathsalong_to_propagate_coltoadd = (edgeseesedge(indexlist)+1)/2;
            pathsalong_to_propagate = [pathsalong_to_propagate(paths2,:) pathsalong_to_propagate_coltoadd];
        end
    end

end

if functionversion == 2
	elapsedtimehodpaths = etime(clock,t00);
    outpar = elapsedtimehodpaths;

	if filehandlingparameters.savehodpaths == 1
    	eval(['save(''',desiredname,''',''hodpaths'',''hodpathsalongplane'',''EDinputdatahash'',''elapsedtimehodpaths'');'])
	end
end
