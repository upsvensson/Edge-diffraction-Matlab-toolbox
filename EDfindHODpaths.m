function [hodpaths,EDinputdatahash] = EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
difforder,EDversionnumber)
% EDfindHODpaths finds the higher order diffraction paths for a convex
% object.
%
% Input parameters:
%   edgeseesedge        Matrix, [nedges,nedges], of 0 and 1
%   visedgesfroms       Matrix, [nedges,nsources], of 0 and 1
%   visedgesfromr       Matrix, [nedges,nreceivers], of 0 and 1
%   difforder
%   EDversionnumber
%
% Output parameters:
%   hodpaths            Cell variable, hodpaths{2} has size [npaths,2] etc
%   EDinputdatahash
% 
% Uses the function Datahash from Matlab Central
% 
% Peter Svensson (peter.svensson@ntnu.no) 15 Mar 2018
% 
% hodpaths = EDfindHODpaths(edgeseesedge,visedgesfroms,visedgesfromr,...
% difforder,EDversionnumber);

% 15 Mar 2018 First version
% 16 Mar 2018 Expanded to several sources and receivers

EDinputdatastruct = struct('difforder',difforder,...
        'edgeseespartialedge',edgeseesedge,...
        'vispartedgesfroms',visedgesfroms,'vispartedgesfromr',visedgesfromr,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

nedges = size(edgeseesedge,2);
nsources = size(visedgesfroms,2);
nreceivers = size(visedgesfromr,2);

hodpaths = cell(difforder,nreceivers,nsources);

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
            ivec = find(edgeseesedge_lastround(:,pathspattern_to_propagate(:,norder-1),ii));
            [paths1,paths2] = ind2sub([nedges,npaths],ivec);
            hodpaths{norder,ii,isou} = [pathspattern_to_propagate(paths2,:) paths1];
        end

        ivec = find(edgeseesedge(:,pathspattern_to_propagate(:,norder-1)));
        [paths1,paths2] = ind2sub([nedges,npaths],ivec);
        pathspattern_to_propagate = [pathspattern_to_propagate(paths2,:) paths1];
    end

end
