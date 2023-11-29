function [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
    EDfindconvexGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
    EDversionnumber,filehandlingparameters)
% EDfindconvexGApaths - Finds all the first-order specular and (possibly)
% first-order diffraction paths for a convex object.
% 
% From v0.4 of the EDtoolbox, the input parameter list changed.
%
% Input parameters:
%   planedata, edgedata, Sdata, Rdata       structs
% 	controlparameters                       struct (only fields .difforder
%                                           and .directsound are used)
%   EDversionnumber
%   filehandlingparameters  filehandlingparameters is a struct which
%                           contains the field showtext.
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
%   elapsedtimefindpaths
%                           This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimefindpaths has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        If an existing file was found that was reused,
%                           then the reused file name is given here. 
%                           If no existing file could be reused then this 
%                           variable is empty. 
%
% Uses functions  EDfindis EDchkISvisible EDrecycleresultfiles from EDtoolbox
% Uses function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 29 Nov. 2023
%
% [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
%    EDfindconvexGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
%    EDversionnumber,filehandlingparameters);

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
% 1 Feb 2018 Fixed a bug; the directsoundlist sometimes got a horizontal
% format.
% 2 Feb 2018 Fixed a small bug: if doaddsources was set to 0, and no
% direct sound obstructions were possible, or no specular reflection was possible,
% then an error occurred.
% 8 Feb 2018 Introduced the EDinputdatahash
% 6 Feb 2019 Fixed a bug which occurred if nsources > 1 and nreceivers > 1.
% The number of potential specular reflections was erroneously found to be the last
% number of possible reflections in the for-loop around line 160.
% 16 Mar 2021 Empirically established corner-hit scaling of the direct
% sound and specular reflections were introduced. Also, an extra test was
% introduced for edge hits and corner hits, so that the direct sound was not
% slipping through in special cases. 
% 28 Oct. 2022 Fixed a bug: the directsound visibility was not done
% correctly when "doaddsources = 0".
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 9 Oct. 2023 Corrected the description of input and output parameters
% 12 Oct. 2023 Made small modifications so a free-field case could be
% handled.
% 27 Oct. 2023 Substantial change to the input parameters.
% 30 Oct. 2023 Fine-tuned the EDinputdatahash
% 9 Nov. 2023 Added the piston coordinates and piston gauss order to the
% EDinputdatahash.
% 29 Nov. 2023 Adapted to the name change of the field pistongausspoints to
% pistongaussorder

t00 = clock;

showtext = filehandlingparameters.showtext;
difforder = controlparameters.difforder;
directsound = controlparameters.directsound;

EDinputdatastruct = struct('corners',planedata.corners,'planecorners',...
    planedata.planecorners,'offedges',edgedata.offedges,...
    'sources',Sdata.coordinates,'Snedgesubs',Sdata.nedgesubs,...
    'pistoncoordinates',Sdata.pistoncornercoordinates,...
    'pistongaussorder',Sdata.pistongaussorder,...
    'receivers',Rdata.coordinates,'Rnedgesubs',Rdata.nedgesubs,...
    'difforder',difforder,'directsound',directsound,'doallSRcombinations',...
    Sdata.doallSRcombinations,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. Should the data be saved in a file? 

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_paths',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_paths.mat'];

if foundmatch == 1
%        eval(['load ''',existingfilename,''''])
    eval(['load(''',existingfilename,''')'])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    elapsedtimefindpaths_new = etime(clock,t00);
    elapsedtimefindpaths = [elapsedtimefindpaths_new elapsedtimefindpaths];
    return
end

%-----------------------------------------------------------------
nplanes = size(planedata.planecorners,1);
nedges = size(edgedata.edgecorners,1);
nsources = size(Sdata.coordinates,1);
nreceivers = size(Rdata.coordinates,1);

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
        if Sdata.doallSRcombinations == 1
            for jj = 1:nsources
                tempvec = Rdata.vispartedgesfromr(:,ii).*Sdata.vispartedgesfroms(:,jj);
                edgeisactive = edgeisactive + tempvec;        
                diffpaths(ii,jj,:)=( tempvec>0 );
            end
        else
           jj = ii; 
           tempvec = Rdata.vispartedgesfromr(:,ii).*Sdata.vispartedgesfroms(:,jj);
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
% visplanesfroms has size [nplanes,nsources]
% visplanesfromr has size [nplanes,nreceivers]

if showtext >= 2
   disp(['      Finding first-order specular reflection paths']) 
end

min_number_elements = min([nsources nreceivers nplanes]);

if nplanes > 0
    if nsources == min_number_elements
        possibleSPR = [];
        for ii = 1:nsources
            visplanesfromoneS = Sdata.visplanesfroms(:,ii);
            tempmatrix = Rdata.visplanesfromr.*visplanesfromoneS(:,ones(1,nreceivers));
            ivpotential = find(tempmatrix);
            npotentialIS = length(ivpotential);
            if npotentialIS > 0
                [Pnumber_potentialIS,Rnumber_potentialIS] = ind2sub([nplanes,nreceivers], ivpotential);
                possibleSPR = [possibleSPR;[ii*ones(npotentialIS,1) Pnumber_potentialIS Rnumber_potentialIS]];
            end
        end
    else 
        % 27 Nov. 2023 preallocated space for cases when nsources =
        % nreceivers because quite possibly, there will be many cases
        if nsources == nreceivers
            possibleSPR = zeros(nreceivers^2,3);
            startrow = 1;            
            for ii = 1:nreceivers
                visplanesfromoneR = Rdata.visplanesfromr(:,ii);
                % 30 Oct. 2023 Fixed a bug on the line below. visplanesfromoneR
                % is not a field in Rdata. This section has not been run often.
                %            tempmatrix = Sdata.visplanesfroms.*Rdata.visplanesfromoneR(:,ones(1,nsources));
                tempmatrix = Sdata.visplanesfroms.*visplanesfromoneR(:,ones(1,nsources));
                ivpotential = find(tempmatrix);
                npotentialIS = length(ivpotential);
                if npotentialIS > 0
                    [Pnumber_potentialIS, Snumber_potentialIS ] = ind2sub([nplanes,nsources], ivpotential);
                    matrixaddition = [ Snumber_potentialIS  Pnumber_potentialIS ii*ones(npotentialIS,1)];
                    nnewrows = size(matrixaddition,1);
                    endrow = startrow + nnewrows -1;
                    if endrow <= size(possibleSPR,1)
                        possibleSPR(startrow:endrow,:) =  matrixaddition;
                    else
                        error('ERROR: need to expand possibleSPR')
                    end
                    startrow = endrow + 1;
                end
            end
        else
            possibleSPR = [];
            for ii = 1:nreceivers
                visplanesfromoneR = Rdata.visplanesfromr(:,ii);
                % 30 Oct. 2023 Fixed a bug on the line below. visplanesfromoneR
                % is not a field in Rdata. This section has not been run often.
                %            tempmatrix = Sdata.visplanesfroms.*Rdata.visplanesfromoneR(:,ones(1,nsources));
                tempmatrix = Sdata.visplanesfroms.*visplanesfromoneR(:,ones(1,nsources));
                ivpotential = find(tempmatrix);
                npotentialIS = length(ivpotential);
                if npotentialIS > 0
                    [Pnumber_potentialIS, Snumber_potentialIS ] = ind2sub([nplanes,nsources], ivpotential);
                    possibleSPR = [possibleSPR;[ Snumber_potentialIS  Pnumber_potentialIS ii*ones(npotentialIS,1)]];            
                end
            end
        end
    end
else
    possibleSPR = [];
end

if Sdata.doallSRcombinations == 0
    if ~isempty(possibleSPR)
        iv = find(possibleSPR(:,1) ~= possibleSPR(:,3));
        possibleSPR(iv,:) = [];
        npotentialIS = size(possibleSPR,1);
    else
        npotentialIS = 0;
    end
else
   npotentialIS = size(possibleSPR,1);    
end

if npotentialIS > 0
    coords_potentialIS = EDfindis(Sdata.coordinates(possibleSPR(:,1),:),possibleSPR(:,2),planedata.planeeqs);
    [hitplanes,hitpoints,edgehits,edgehitpoints,edgehitnumbers,cornerhits,cornerhitpoints,cornerhitnumbers] = ...
        EDchkISvisible(coords_potentialIS,Rdata.coordinates(possibleSPR(:,3),:),...
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
        % 16 March 2021: for a 90 degree corner (thin sheet or cube),
        % a smooth transition results across a corner if the specular
        % reflection is scaled to 1/3 (rather than 1/2) 
        specreflamp = [specreflamp;1/3*ones(size(cornerhits))];        
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
% visplanesfroms has size [nplanes,nsources]
% visplanesfromr has size [nplanes,nreceivers]

if directsound ~= 0
    if showtext >= 2
       disp(['      Finding direct sound paths']) 
    end

    if nplanes == 0
        possibleSPR_obstruct = [];
    else
        if nsources == min_number_elements
            possibleSPR_obstruct = [];
            for ii = 1:nsources
                visplanesfromoneS = Sdata.visplanesfroms(:,ii);
                ivpotential = find(Rdata.visplanesfromr ~= visplanesfromoneS(:,ones(1,nreceivers)));
                npotentialobstruct = length(ivpotential);
                if npotentialobstruct > 0
                    [Pnumber_potentialobstruct, Rnumber_potentialobstruct] = ind2sub([nplanes,nreceivers], ivpotential);
                    possibleSPR_obstruct = [possibleSPR_obstruct;[ii*ones(npotentialobstruct,1) Pnumber_potentialobstruct Rnumber_potentialobstruct]];
                end
            end
        else 
            possibleSPR_obstruct = [];
            for ii = 1:nreceivers
                visplanesfromoneR = Rdata.visplanesfromr(:,ii);
                % 30 Oct. 2023 Fixed a bug on the line below. visplanesfromoneR
                % is not a field in Rdata. This section has not been run often.
%                ivpotential = find(Sdata.visplanesfroms ~= Rdata.visplanesfromoneR(:,ones(1,nsources)) );
                ivpotential = find(Sdata.visplanesfroms ~= visplanesfromoneR(:,ones(1,nsources)) );
                npotentialobstruct = length(ivpotential);
                if npotentialobstruct > 0
                    [Pnumber_potentialobstruct, Snumber_potentialobstruct] = ind2sub([nplanes,nsources], ivpotential);
                    possibleSPR_obstruct = [possibleSPR_obstruct;[ Snumber_potentialobstruct  Pnumber_potentialobstruct ii*ones(npotentialobstruct,1)]];
                end
            end
        end
    end

    if Sdata.doallSRcombinations == 0
        directsoundOK = eye(nreceivers);
        if ~isempty(possibleSPR_obstruct)
            iv = find(possibleSPR_obstruct(:,1) ~= possibleSPR_obstruct(:,3));
            possibleSPR_obstruct(iv,:) = [];
            npotentialIS = size(possibleSPR,1);
            npotentialobstruct = size(possibleSPR_obstruct,1);
        else
            npotentialobstruct = 0;
        end
    else
        directsoundOK = ones(nreceivers,nsources);
        npotentialobstruct = size(possibleSPR_obstruct,1);
    end

    if npotentialobstruct > 0

        [hitplanes,hitpoints,edgehits,edgehitpoints,edgehitnumbers,cornerhits,cornerhitpoints,cornerhitnumbers] = ...
            EDchkISvisible(Sdata.coordinates(possibleSPR_obstruct(:,1),:),Rdata.coordinates(possibleSPR_obstruct(:,3),:),...
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
            % 16 Mar 2021: for the direct sound to be (half) visible, the
            % edge must be seen by the source and the receiver. It is a bit
            % inefficient to run this after the direct sound amplitude has
            % already been modified, but it works.
            ivedgehit = find(edgehitnumbers);
            for kk = 1:length(ivedgehit)
               corners_of_one_plane = planedata.planecorners...
                   (possibleSPR_obstruct(ivedgehit(kk),2),:);
               corners_of_one_plane = [corners_of_one_plane corners_of_one_plane(1)];
               corners_of_one_edge = corners_of_one_plane(edgehitnumbers(ivedgehit(kk)):edgehitnumbers(ivedgehit(kk))+1);              
               [lia,locb] = ismember(corners_of_one_edge,edgedata.edgecorners,'rows');
               if lia == 0
                [lia,locb] = ismember(fliplr(corners_of_one_edge),edgedata.edgecorners,'rows');
               end
               edgenumber_that_was_hit = locb(1);
               if Sdata.vispartedgesfroms(edgenumber_that_was_hit)*Rdata.vispartedgesfromr(edgenumber_that_was_hit) == 0
                    directsoundOK(iv(kk))=0; 
               end
            end
        end
        if ~isempty(cornerhits)
            cornerhitlist = possibleSPR_obstruct(cornerhits,:);
            iv = sub2ind([nreceivers,nsources],cornerhitlist(:,3),cornerhitlist(:,1));  
            % 16 March 2021: for a 90 degree corner (thin sheet or cube),
            % a smooth transition results across a corner if the direct
            % sound is scaled to 2/3 (rather than 1/2)
            directsoundOK(iv)=2/3*directsoundOK(iv);
            % 16 Mar 2021: for the direct sound to be (partly) visible, the
            % corner must be seen by the source and the receiver. It is a bit
            % inefficient to run this after the direct sound amplitude has
            % already been modified, but it works.
            ivcornerhit = find(cornerhitnumbers);
            for kk = 1:length(ivcornerhit)
                corners_of_one_plane = planedata.planecorners...
                   (possibleSPR_obstruct(ivcornerhit(kk),2),:);
                the_hit_corner = corners_of_one_plane(cornerhitnumbers(ivcornerhit(kk)));
                connectededges_to_hit_corner = find( sum( (edgedata.edgecorners == the_hit_corner).' ).');
                if sum(Sdata.vispartedgesfroms(connectededges_to_hit_corner).*Rdata.vispartedgesfromr(connectededges_to_hit_corner)) == 0
                    directsoundOK(iv(kk))=0;
                end
            end
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
firstorderpathdata.specrefllist     = [validsounumber(:) validrecnumber(:) specreflamp(:)];
firstorderpathdata.diffpaths        = diffpaths;
firstorderpathdata.edgeisactive     = edgeisactive;
firstorderpathdata.directsoundlist  = [Snumber_directsoundOK(:) Rnumber_directsoundOK(:) dirsoundamp(:)];
firstorderpathdata.ncomponents      = numberofcomponents;

elapsedtimefindpaths = etime(clock,t00);
outpar = elapsedtimefindpaths;

if filehandlingparameters.savepathsfile == 1
	eval(['save(''',desiredname,''',''firstorderpathdata'',''EDinputdatahash'',''elapsedtimefindpaths'');'])
end
