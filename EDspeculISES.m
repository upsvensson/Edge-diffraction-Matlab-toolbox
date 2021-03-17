function 	[validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDspeculISES(corners,planecorners,planeeqs,planenvecs,...
    S,R,lengthNspecmatrix,specorder,visplanesfromoner,planeisthin,canplaneobstruct,minvals,maxvals,ncornersperplanevec,planeseesplane,rearsideplane,modeltype,showtext)
% ESIE2speculISES - Finds the valid specular reflections by checking visibility and obstruction.
% Finds the valid specular reflections, given a list of
% potential IS combinations, by checking visibility and obstruction.
%
% Input parameters:
%   corners, planecorners,planeeqs,planenvecs          
%   S, R
%   lengthNspecmatrix,
%   specorder
%   visplanesfromoner
%   planeisthin,
%   canplaneobstruct,minvals,maxvals,ncornersperplanevec,planeseesplane,rearsideplane,modeltype,showtext
%
% GLOBAL parameters:
%   POTENTIALISES, ISCOORDS, ORIGINSFROM, IVNSPECMATRIX 
%
% Output parameters
%   validISlist         Matrix [nreflections,specorder] of all reflections, with one row for each reflection.
%				        The reflecting plane numbers are given, one for each column.
%   validIScoords       Matrix [nreflections,3] of all image source coordinates.
%   allreflpoints       Matrix [nreflections,(specorder-1)*3] of all hit points in planes.
%   listguide           Matrix [nactiveorders,3] which for each row gives
%                       1. the number of valid reflections for one reflection order
%                       2. the first, and 3. the last
%                       row in validISlist etc which contain reflections of that order.
%                       Which reflection order each row refers to is given
%                       by listofreflorder.
%   listofreflorder     List [nactiveorder,3] which for each row gives the specular
%                       reflection order that the corresponding row in listguide refers to.
%
% Uses functions EDchkISvisible EDcheckobstrpaths
%
% Peter Svensson (peter.svensson@ntnu.no) 16 March 2021
%
% [validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDspeculISES(corners,planecorners,planeeqs,planenvecs,...
%     S,R,lengthNspecmatrix,specorder,visplanesfromoner,planeisthin,canplaneobstruct,minvals,maxvals,ncornersperplanevec,planeseesplane,rearsideplane,modeltype,showtext);

% 9 Nov. 2004 Functioning version
% 28 Nov. 2017 Copied to EDtoolbox. Introduced the input parameter
% modeltype: for a convex scatterer, obstruction test is not needed for
% the first-order specular reflection. Also introduced the non-global input
% parameter showtext.
% 29 Nov. 2017 Changed so that a thin-plate model doesn't need obstruction
%              check for the specular reflections.
% 16 Mar 2021 Adapted to a change in EDchkISvisible.

global POTENTIALISES ISCOORDS ORIGINSFROM IVNSPECMATRIX

% [ncorners,slask] = size(corners);
% [nplanes,ncornersperplane] = size(planecorners);

% [n1,n2] = size(POTENTIALISES);
n2 = size(POTENTIALISES,2);

if n2 < specorder
    if showtext >= 2
        disp('   WARNING: The ISEStree-file does not contain high enough order to support')
        disp(['            spec. refl. order ',int2str(specorder)])
    end
    specorder_maxpossible = n2;
else
    specorder_maxpossible = specorder;
end

%----------------------------------------------
% Find all IS that can be propagated to higher order.
% These must come via planes that S can see, and will be stored in xis1.

% thinlist = find(planeisthin==1);

%   ###########################################
%   #                                         #
%   #         First order or higher           #
%   #                                         #
%   ###########################################

validISlist = [];
validIScoords = [];
allreflpoints = [];
listguide        = zeros(specorder,3);
listofreflorder = zeros(specorder,1);
listguide(1,2)   = 1;

obstructtestneeded = (sum(canplaneobstruct)~=0) && strcmp(modeltype,'convex_ext')==0 && isempty(strfind(modeltype,'plate'));

for norder = 1:specorder_maxpossible

    if showtext >= 3
		disp(['      Reflection order ',int2str(norder)])
    end
    
    % Start with all the theoretically possible IS
    masterivlist = IVNSPECMATRIX(1:lengthNspecmatrix(norder),norder);
        
    possiblecombs = POTENTIALISES(masterivlist,1:norder);
    
    % Pick out only those IS that come via a last refl plane that the
    % receiver is in front of
%    okcombs = find(visplanesfromoner(possiblecombs(:,norder))==1);
    okcombs = visplanesfromoner(possiblecombs(:,norder))==2;
    masterivlist = masterivlist(okcombs);
    possiblecombs = POTENTIALISES(masterivlist,1:norder); 
    possibleIS = ISCOORDS(masterivlist,:);
    nis = length(masterivlist);
 
    if showtext >= 3
		disp(['         ',int2str(nis),' IS found initially for order ',int2str(norder)])
        if showtext >= 5
           possiblecombs    
        end
    end

    % Go through iteratively the reflection points, and check if they are
    % inside the reflection plane. Start with the last reflection plane and
    % work backwards.
    for jj = norder:-1:1
        if nis > 0

            if jj == norder
                fromcoords = possibleIS;
                tocoords = R; 
               else    
                
                eval(['tocoords = reflpoints',int2str(jj+1),';'])    
                ivreflist = ORIGINSFROM(masterivlist);
                for kk = jj:norder-2
                    ivreflist = ORIGINSFROM(ivreflist);
                end
                fromcoords = ISCOORDS(ivreflist,:);
            end
            
%             [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ESIE2chkISvisible(fromcoords,tocoords,planeeqs(possiblecombs(:,jj),:),planenvecs(possiblecombs(:,jj),:),minvals(possiblecombs(:,jj),:),...
% 				maxvals(possiblecombs(:,jj),:),planecorners(possiblecombs(:,jj),:),corners,ncornersperplanevec(possiblecombs(:,jj)));
%             [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ESIE2chkISvisible(fromcoords,tocoords,planeeqs(possiblecombs(:,jj),4),planenvecs(possiblecombs(:,jj),:),minvals(possiblecombs(:,jj),:),...
% 				maxvals(possiblecombs(:,jj),:),planecorners(possiblecombs(:,jj),:),corners,ncornersperplanevec(possiblecombs(:,jj)));
            [hitplanes,reflpoints,edgehits,edgehitpoints,edgehitnumbers,cornerhits,cornerhitpoints,cornerhitnumbers] = ...
                EDchkISvisible(fromcoords,tocoords,planeeqs(possiblecombs(:,jj),4),planenvecs(possiblecombs(:,jj),:),minvals(possiblecombs(:,jj),:),...
				maxvals(possiblecombs(:,jj),:),planecorners(possiblecombs(:,jj),:),corners,ncornersperplanevec(possiblecombs(:,jj)));
            if ~isempty(edgehits) || ~isempty(cornerhits)
                disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                disp('         handled correctly yet.')
            end
            eval(['reflpoints',int2str(jj),' = reflpoints;'])
            
            masterivlist = masterivlist(hitplanes);
            possiblecombs = POTENTIALISES(masterivlist,1:norder);  
            possibleIS = ISCOORDS(masterivlist,:);
            if jj < norder
                for kk = jj+1:norder
                   eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(hitplanes,:);'])
                end
            end
            nis = length(masterivlist);
	
            if showtext >= 3
				disp(['         ',int2str(nis),' IS survived the visibility test in refl plane number ',int2str(jj)])
                if showtext >= 5
                    possiblecombs    
                end
            end
        end    

    end

    if obstructtestneeded && nis > 0
        % Check obstructions for all the paths: S -> plane1 -> plane2 -> ...
        % -> planeN -> R

        for jj = 1:norder+1
 
            if nis > 0
                if jj==1
                    fromcoords = S;    
                    startplanes = [];    
                else
                    startplanes = possiblecombs(:,jj-1);
                    eval(['fromcoords = reflpoints',int2str(jj),';'])
                end
                if jj == norder+1
                    tocoords = R;
                    endplanes = [];    
                else
                    eval(['tocoords = reflpoints',int2str(jj),';'])    
                    endplanes = possiblecombs(:,jj);    
                end
                [nonobstructedpaths,nobstructions] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);

                if nobstructions > 0
                    masterivlist = masterivlist(nonobstructedpaths);
                    possiblecombs = POTENTIALISES(masterivlist,1:norder);  
                    possibleIS = ISCOORDS(masterivlist,:);
                    for kk = 1:norder
                        eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(nonobstructedpaths,:);'])    
                    end
                    nis = length(masterivlist);
                end
                
                if showtext >= 3
			        disp(['         ',int2str(nis),' IS survived the obstruction test for path ',int2str(jj)])
                    if showtext >= 5
                       possiblecombs    
                    end
                end
            end
        
        end
    
    end
    if nis > 0
        if ~isempty(validISlist)
            [n1,n2] = size(validISlist);
            n4 = size(possiblecombs,2);
            validISlist = [[validISlist zeros(n1,n4-n2)];(possiblecombs)];
        else
            validISlist = possiblecombs;
        end
        validIScoords = [validIScoords;possibleIS];

        newestreflpoints = [];
        for kk = 1:norder
            eval(['newestreflpoints = [newestreflpoints reflpoints',int2str(kk),'];'])    
        end
        
        if ~isempty(allreflpoints)
            [n1,n2] = size(allreflpoints);
            n4 = size(newestreflpoints,2);
            allreflpoints = [[allreflpoints zeros(n1,n4-n2)];newestreflpoints];    
        else
            allreflpoints = newestreflpoints;    
        end
    end

    listguide(norder,1) = nis;
    listguide(norder,3) = listguide(norder,2)+nis-1;
    if norder < specorder
        listguide(norder+1,2) = listguide(norder,3)+1;    
    end
    listofreflorder(norder) = norder;
    
end

iv = find(listguide(:,3)==0);
listguide(iv,2) = zeros(size(iv));

[n1,n2] = size(validISlist);
if n2 < specorder
    validISlist = [validISlist zeros(n1,specorder-n2)];    
end

[n1,n2] = size(allreflpoints);
if n2 < specorder*3
    allreflpoints = [allreflpoints zeros(n1,specorder*3-n2)];    
end

iv = find(listguide(:,1)==0);
listguide(iv,:) = [];
listofreflorder(iv) = [];
