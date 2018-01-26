function firstorderpathdata = EDfindconvexGApaths(planedata,edgedata,...
    sources,visplanesfromS,vispartedgesfromS,receivers,visplanesfromR,...
    vispartedgesfromR,difforder,directsound,doallSRcombinations,showtext)
% EDfindconvexGApaths - Finds all the first-order specular and (possibly)
% first-order diffraction paths for a convex object.
%
% Input parameters:
%   planedata, edgedata,edgetoedgedata   structs
% 	sources
% 	visplanesfromS, vispartedgesfromS       From the Sdata struct
%   receivers
% 	visplanesfromR, vispartedgesfromR       From the Rdata struct
%   difforder
%   directsound
%   doallSRcombinations
%   showtext                (optional) Default: 0
%
% Output parameters:
%   firstorderpathdata      struct with the fields:
%   .specrefllist           matrix, [nIS,2] with all the valid specular
%                           reflections. First column has the source
%                           numbers and the second column has the receiver
%                           numbers.
%   .specreflIScoords       matrix, [nIS,3] with IS coordinates
%   .diffpaths              matrix, [nreceivers,nsources,nedges] with
%                           logical 0 or 1
%   .edgeisactive           vector, [nedges,1] with logical 0 or 1
%   .directsoundlist        matrix, [ncomponents,2] with all the direct
%                           sound components that are visible. First column
%                           has the source number and the second column has
%                           the receiver number.
%   .ncomponents            vector, [1,3], with the number of components
%                           for the direct sound, specular reflections, diffraction.
%
% Uses functions  EDfindis EDchkISvisible 
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
% Peter Svensson (peter.svensson@ntnu.no) 26 Jan. 2018
%
% firstorderpathdata = EDfindconvexGApaths(planedata,edgedata,edgetoedgedata,...
% sources,visplanesfromS,vispartedgesfromS,receivers,visplanesfromR,vispartedgesfromR,...
% difforder,directsound,doallSRcombinations,showtext)

% 27 Dec. 2017 First start
% 28 Dec. 2017 Functioning version for diff
% 11 Jan. 2018 First complete version
% 12 Jan. 2018 Small bug.fixes
% 15 Jan. 2018 Added the direct sound amplitude: 1, 0.5 or 0.25 for edge
% and corner hits. Also for the specular reflections.
% 17 Jan 2018 Introduced the input parameter difforder, so one can skip the
% diffraction, if wanted.
% 18 Jan 2018 Fixed an important bug for the direct sound, which happened
% when there were many sources and many receivers.
% 19 Jan 2018 Fixed a bug: error occured if there were no specular
% reflections.
% 23 Jan 2018 Fixed a bug with the direct sound amplitude. See explanation
% in source code.
% 25 Jan 2018 Fixed a bug: preciously, the direct sound was calculated even if
% .directsound = 0.
% 26 Jan 2018: V 0.107: introduced the doallSRcombinations parameter

if nargin < 13
   showtext = 0; 
end

% planedata.corners = size(planedata.corners,1);
nplanes = length(planedata.planeisthin);
nedges = length(edgedata.closwedangvec);
nsources = size(visplanesfromS,2);
nreceivers = size(visplanesfromR,2);

numberofcomponents = zeros(1,3);

%--------------------------------------------------------------------------
% First the first-order diffraction paths. Easy for a convex scatterer: if
% both the source and the receiver can see an edge, then it is fully
% visible. Otherwise it is not at all visible.

if difforder > 0
    if showtext >= 2
       disp(['      Finding first-order diffraction paths']) 
    end
    diffpaths = zeros(nreceivers,nsources,nedges);
    edgeisactive = zeros(nedges,1,'uint8');

    for ii = 1:nreceivers
        if doallSRcombinations == 1
            for jj = 1:nsources
                tempvec = vispartedgesfromR(:,ii).*vispartedgesfromS(:,jj);
                edgeisactive = edgeisactive + tempvec;        
                diffpaths(ii,jj,:)=( tempvec>0 );
            end
        else
           jj = ii; 
           tempvec = vispartedgesfromR(:,ii).*vispartedgesfromS(:,jj);
           edgeisactive = edgeisactive + tempvec;        
           diffpaths(ii,jj,:)=( tempvec>0 );
        end
    end

    edgeisactive = (edgeisactive>0);
else
    diffpaths = [];
    edgeisactive = [];
end

%--------------------------------------------------------------------------
% Then the first-order specular paths. 
%
% We find all potentially possible S-P-R combos.
%
% visplanesfromS has size [nplanes,nsources]
% visplanesfromR has size [nplanes,nreceivers]

if showtext >= 2
   disp(['      Finding first-order specular reflection paths']) 
end

min_number_elements = min([nsources nreceivers nplanes]);

if nsources == min_number_elements
    possibleSPR = [];
    for ii = 1:nsources
        visplanesfromoneS = visplanesfromS(:,ii);
        tempmatrix = visplanesfromR.*visplanesfromoneS(:,ones(1,nreceivers));
        ivpotential = find(tempmatrix);
        npotentialIS = length(ivpotential);
        if npotentialIS > 0
            [Pnumber_potentialIS,Rnumber_potentialIS] = ind2sub([nplanes,nreceivers], ivpotential);
            possibleSPR = [possibleSPR;[ii*ones(npotentialIS,1) Pnumber_potentialIS Rnumber_potentialIS]];
        end
    end
else 
    possibleSPR = [];
    for ii = 1:nreceivers
        visplanesfromoneR = visplanesfromR(:,ii);
        tempmatrix = visplanesfromS.*visplanesfromoneR(:,ones(1,nsources));
        ivpotential = find(tempmatrix);
        npotentialIS = length(ivpotential);
        if npotentialIS > 0
            [Pnumber_potentialIS, Snumber_potentialIS ] = ind2sub([nplanes,nsources], ivpotential);
            possibleSPR = [possibleSPR;[ Snumber_potentialIS  Pnumber_potentialIS ii*ones(npotentialIS,1)]];
        end
    end
end

if doallSRcombinations == 0
    iv = find(possibleSPR(:,1) ~= possibleSPR(:,3));
    possibleSPR(iv,:) = [];
    npotentialIS = size(possibleSPR,1);
end

if npotentialIS > 0
    coords_potentialIS = EDfindis(sources(possibleSPR(:,1),:),possibleSPR(:,2),planedata.planeeqs);

    [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ...
        EDchkISvisible(coords_potentialIS,receivers(possibleSPR(:,3),:),...
        planedata.planeeqs(possibleSPR(:,2),4),planedata.planeeqs(possibleSPR(:,2),1:3),...
        planedata.minvals(possibleSPR(:,2),:),planedata.maxvals(possibleSPR(:,2),:),...
        planedata.planecorners(possibleSPR(:,2),:),planedata.corners,planedata.ncornersperplanevec(possibleSPR(:,2)));
    
    specreflamp = ones(size(hitplanes));    
    validIScoords = coords_potentialIS(hitplanes,:);
    validsounumber = possibleSPR(hitplanes,1);
    validrecnumber = possibleSPR(hitplanes,3);
    
    if ~isempty(edgehits)
        validIScoords  = [validIScoords;coords_potentialIS(edgehits,:)];
        validsounumber = [validsounumber;possibleSPR(edgehits,1)];
        validrecnumber = [validrecnumber;possibleSPR(edgehits,3)];
        specreflamp = [specreflamp;0.5*ones(size(edgehits))];
    end
    if ~isempty(cornerhits)
        validIScoords  = [validIScoords;coords_potentialIS(cornerhits,:)];
        validsounumber = [validsounumber;possibleSPR(cornerhits,1)];
        validrecnumber = [validrecnumber;possibleSPR(cornerhits,3)];
        specreflamp = [specreflamp;0.25*ones(size(cornerhits))];        
    end
    
    numberofcomponents(2) = length(specreflamp);
        
else
    validIScoords = [];
    validsounumber = [];
    validrecnumber = [];  
    specreflamp = [];
end
    
    
%--------------------------------------------------------------------------
% Then the direct sound 
%
% We find all potentially obstructing S-P-R combos
%
% visplanesfromS has size [nplanes,nsources]
% visplanesfromR has size [nplanes,nreceivers]

if directsound ~= 0
    if showtext >= 2
       disp(['      Finding direct sound paths']) 
    end

    if nsources == min_number_elements
        possibleSPR_obstruct = [];
        for ii = 1:nsources
            visplanesfromoneS = visplanesfromS(:,ii);
            ivpotential = find(visplanesfromR ~= visplanesfromoneS(:,ones(1,nreceivers)));
            npotentialobstruct = length(ivpotential);
            if npotentialobstruct > 0
                [Pnumber_potentialobstruct, Rnumber_potentialobstruct] = ind2sub([nplanes,nreceivers], ivpotential);
                possibleSPR_obstruct = [possibleSPR_obstruct;[ii*ones(npotentialobstruct,1) Pnumber_potentialobstruct Rnumber_potentialobstruct]];
            end
        end
    else 
        possibleSPR_obstruct = [];
        for ii = 1:nreceivers
            visplanesfromoneR = visplanesfromR(:,ii);
            ivpotential = find(visplanesfromS ~= visplanesfromoneR(:,ones(1,nsources)) );
            npotentialobstruct = length(ivpotential);
            if npotentialobstruct > 0
                [Pnumber_potentialobstruct, Snumber_potentialobstruct] = ind2sub([nplanes,nsources], ivpotential);
                possibleSPR_obstruct = [possibleSPR_obstruct;[ Snumber_potentialobstruct  Pnumber_potentialobstruct ii*ones(npotentialobstruct,1)]];
            end
        end
    end

    if doallSRcombinations == 0
        iv = find(possibleSPR_obstruct(:,1) ~= possibleSPR_obstruct(:,3));
        possibleSPR_obstruct(iv,:) = [];
        npotentialIS = size(possibleSPR,1);
        directsoundOK = eye(nreceivers);
        npotentialobstruct = size(possibleSPR_obstruct,1);
    else
        directsoundOK = ones(nreceivers,nsources);        
    end


    if npotentialobstruct > 0

        [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ...
            EDchkISvisible(sources(possibleSPR_obstruct(:,1),:),receivers(possibleSPR_obstruct(:,3),:),...
            planedata.planeeqs(possibleSPR_obstruct(:,2),4),planedata.planeeqs(possibleSPR_obstruct(:,2),1:3),...
            planedata.minvals(possibleSPR_obstruct(:,2),:),planedata.maxvals(possibleSPR_obstruct(:,2),:),...
            planedata.planecorners(possibleSPR_obstruct(:,2),:),planedata.corners,planedata.ncornersperplanevec(possibleSPR_obstruct(:,2)));

        obstructlist = possibleSPR_obstruct(hitplanes,:);    
        iv = sub2ind([nreceivers,nsources],obstructlist(:,3),obstructlist(:,1));    
        directsoundOK(iv)=0;

        if ~isempty(edgehits)
            edgehitlist = possibleSPR_obstruct(edgehits,:);
            iv = sub2ind([nreceivers,nsources],edgehitlist(:,3),edgehitlist(:,1));  
            % Before 23 Jan 2018, the first line was  used. This lead to that a
            % completely obscured receiver (which got marked by hitplanes)
            % might get a non-zero direct sound if there was an additional
            % edgehit, unrelated to the obscuring plane.
            %         directsoundOK(iv)=0.5;        
            directsoundOK(iv)=0.5*directsoundOK(iv);        
        end
        if ~isempty(cornerhits)
            cornerhitlist = possibleSPR_obstruct(cornerhits,:);
            iv = sub2ind([nreceivers,nsources],cornerhitlist(:,3),cornerhitlist(:,1));  
            % Before 23 Jan 2018, the first line was  used. This lead to that a
            % completely obscured receiver (which got marked by hitplanes)
            % might get a non-zero direct sound if there was an additional
            % cornerhit, unrelated to the obscuring plane.
    %         directsoundOK(iv)=0.75;
            directsoundOK(iv)=0.75*directsoundOK(iv);
        end
    end

    ivdirectsoundOK = find(directsoundOK);

    [Rnumber_directsoundOK,Snumber_directsoundOK] = ind2sub([nreceivers,nsources], ivdirectsoundOK);

    dirsoundamp = directsoundOK(ivdirectsoundOK);

    numberofcomponents(1) = length(ivdirectsoundOK);
else
    Snumber_directsoundOK = [];
    Rnumber_directsoundOK = [];
    dirsoundamp = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pack the output data in a struct

firstorderpathdata = struct;
firstorderpathdata.specreflIScoords = validIScoords;
firstorderpathdata.specrefllist     = [validsounumber validrecnumber specreflamp];
firstorderpathdata.diffpaths        = diffpaths;
firstorderpathdata.edgeisactive     = edgeisactive;
firstorderpathdata.directsoundlist  = [Snumber_directsoundOK Rnumber_directsoundOK dirsoundamp];
firstorderpathdata.ncomponents      = numberofcomponents;
