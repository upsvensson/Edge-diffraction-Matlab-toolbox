function [hodpaths,hodpathsalongplane,EDinputdatahash] = EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
difforder,EDversionnumber)
% EDfindHODpaths finds the higher order diffraction paths for a convex
% object, or a non-convex object where specular reflections are ignored.
%
% Input parameters:
%   edgeseesedge        Matrix, [nedges,nedges], of -1, 0 and 1
%   visedgesfroms       Matrix, [nedges,nsources], of 0 and 1
%   visedgesfromr       Matrix, [nedges,nreceivers], of 0 and 1
%   difforder
%   EDversionnumber
%
% Output parameters:
%   hodpaths            Cell variable, hodpaths{2} has size [npaths,2] etc
%                       Each row in each matrix has a sequence of edge
%                       numbers.
%   hodpathsalongplane  Cell variable, hodpaths{2} has size [npaths,1] etc
%                       Each row in each matrix has values 0 or 1, which tells
%                       for each hod sequence "leg" whether it passes along
%                       a plane or not.
%   EDinputdatahash
% 
% Uses the function Datahash from Matlab Central
% 
% Peter Svensson (peter.svensson@ntnu.no) 1 June 2019
% 
% [hodpaths,hodpathsalongplane,EDinputdatahash] = EDfindHODpaths(edgeseesedge,...
% visedgesfroms,visedgesfromr,difforder,EDversionnumber);

% 15 Mar 2018 First version
% 16 Mar 2018 Expanded to several sources and receivers
% 21 May 2019 Made some changes so that negative values of edgeseesedge can
% be handled. Introduced the output variable hodpathsalongplane.
% 1 June 2019 Fixed bug; hadn't finished the handling of negative values of
% edgeseesedge.

EDinputdatastruct = struct('difforder',difforder,...
        'edgeseespartialedge',edgeseesedge,...
        'vispartedgesfroms',visedgesfroms,'vispartedgesfromr',visedgesfromr,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

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
           pathsalong_to_propagate = edgeseesedge(indexlist);
        else            
            indexlist = sub2ind([nedges nedges],pathspattern_to_propagate(:,end-1),pathspattern_to_propagate(:,end));
            pathsalong_to_propagate_coltoadd = edgeseesedge(indexlist);
            pathsalong_to_propagate = [pathsalong_to_propagate(paths2,:) pathsalong_to_propagate_coltoadd];
        end
    end

end

