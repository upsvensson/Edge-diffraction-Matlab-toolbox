function [hodpaths,hodpathsalongplane,elapsedtimehodpaths,existingfilename] = ...
    EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
    difforder,EDversionnumber,filehandlingparameters)
% EDfindHODpaths - Finds the higher-order diffraction paths, for the IR calculation.
% Finds all the edge sequences, for diffraction orders 2 to difforder, from
% each source to each receiver, for a convex object (or a non-convex
% object where specular reflections are ignored). The results are saved
% in the stage file <filestem>_hodpaths.mat (if
% filehandlingparameters.savehodpaths = 1), which is reused when the input
% data are the same.
%
% Input parameters:
%   edgeseesedge    [nedges,nedges] with the values -1, 0 or 1: the sign of
%                   edgetoedgedata.edgeseespartialedge (1: the path runs
%                   along a plane, -1: it does not)
%   visedgesfroms   [nedges,nsources] 0 or 1, the sign of
%                   Sdata.vispartedgesfroms
%   visedgesfromr   [nedges,nreceivers] 0 or 1, the sign of
%                   Rdata.vispartedgesfromr
%   difforder       The highest diffraction order
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem,
%                   .suppressresultrecycling and .savehodpaths are used.
%
% Output parameters:
%   hodpaths        Cell array {difforder,nreceivers,nsources}. The cell
%                   hodpaths{N,ireceiver,isource} has one row per path of
%                   diffraction order N, with the N edge numbers.
%   hodpathsalongplane  Cell array of the same size. Each row has, for each
%                   of the N-1 legs of the path, 1 if the leg runs along a
%                   plane, otherwise 0.
%   elapsedtimehodpaths
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [hodpaths,hodpathsalongplane,elapsedtimehodpaths,existingfilename] = ...
%    EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
%    difforder,EDversionnumber,filehandlingparameters);

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
% 29 Oct. 2023 CHanged the input parameters; only version 2.
% 14 Oct. 2026 The EDinputdatahash contains a stageversion from
% EDstageversions instead of the overall EDversionnumber. load/save
% without eval.
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).
% 15 Oct. 2026 The list of called functions was replaced by a See also line.

t00 = clock;

% The version of this function, and of the subroutines it uses, are
% defined in EDstageversions. They are part of the hash, so that saved
% stage files are not reused after a change in the computations.
stageversion = EDstageversions('EDfindHODpaths');

EDinputdatastruct = struct('stage','EDfindHODpaths','stageversion',stageversion,...
        'difforder',difforder,...
        'edgeseespartialedge',edgeseesedge,...
        'vispartedgesfroms',visedgesfroms,'vispartedgesfromr',visedgesfromr);
EDinputdatahash = DataHash(EDinputdatastruct);

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
	filedata = load(existingfilename,'hodpaths','hodpathsalongplane','elapsedtimehodpaths');
	hodpaths = filedata.hodpaths;
	hodpathsalongplane = filedata.hodpathsalongplane;
	elapsedtimehodpaths = filedata.elapsedtimehodpaths;
	if ~strcmp(existingfilename,desiredname)
		copyfile(existingfilename,desiredname);
	end
	elapsedtimehodpaths_new = etime(clock,t00);
	elapsedtimehodpaths = [elapsedtimehodpaths_new elapsedtimehodpaths];
	return
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

elapsedtimehodpaths = etime(clock,t00);

if filehandlingparameters.savehodpaths == 1
	save(desiredname,'hodpaths','hodpathsalongplane','EDinputdatahash','elapsedtimehodpaths','stageversion','EDversionnumber');
end
