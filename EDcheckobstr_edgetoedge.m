function [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_edgetoedge(canplaneobstruct,...
    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane)
% EDcheckobstr_edgetoedge - Checks obstruction for a list of edge-to-edge paths.
% Checks if the paths from N starting points (fromcoords) to
% N ending points (tocoords) pass through any of M planes (those marked by 1 in
% canplaneobstruct), i.e., if they are obstructed.
%
% Input parameters:
%   canplaneobstruct,planeseesplane,planeeqs,planenvecs,minvals,maxvals,...
%   planecorners,corners,ncornersperplanevec,rearsideplane
%                       Data that should have been passed on from the
%                       eddatafile.
%   global:
%   REFTOFROMCOSHO  REFTOTOCOSHO STARTPLANES ENDPLANES
%   BIGTOCOORDS BIGFROMCOORDS BIGSTARTPLANES BIGENDPLANES
%
% Output parameters:
%   nonobstructedpaths  A list, [N_nobstructions,1] of the indices of paths that are not
%                       obstructed. The indices refer to the input
%                       lists fromcoords, tocoords, startplanes
%                       endplanes
%   nobstructions       The number of paths (of the N input paths) that
%                       were not obstructed.
%   edgehits            A list, [N_edgehits,1] of the indicies of paths that hit a
%                       plane right at an edge.
%   cornerhits          A list, [N_cornerhits,1] of the indicies of paths that hit a
%                       plane right at a corner.   
%
% Uses EDpoinpla
%
% Peter Svensson (peter.svensson@ntnu.no) 14 March 2021
%
% [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_edgetoedge(canplaneobstruct,planeseesplane,...
%    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane)

% 22 Sep. 2005 Functioning version
% 27 Nov. 2017 Copied from ESIE2toolbox
% 28 Nov. 2017 Removed input parameter planeseesplane after Matlab
%               recommendation
% 14 Mar 2021 Adapted to change of EDpoinpla

global  REFTOFROMCOSHO  REFTOTOCOSHO STARTPLANES ENDPLANES
global BIGTOCOORDS BIGFROMCOORDS BIGSTARTPLANES BIGENDPLANES

nstartpoints = size(REFTOFROMCOSHO,1);
nendpoints = size(REFTOTOCOSHO,1);
npaths = max(max([nstartpoints nendpoints]));
nplanes = length(canplaneobstruct);

% The startplanes and endplanes can either be one-column or two-column
% lists depending on if paths are between specular reflections or edges.

nplanecols1 = size(STARTPLANES,2);
nplanecols2 = size(ENDPLANES,2);

% In addition, we only need to check planes for which canplaneobstruct = 1.

if nplanes < 256
    maxlistofplanestocheck = uint8(find( canplaneobstruct>0).');
elseif nplanes < 65536
    maxlistofplanestocheck = uint16(find( canplaneobstruct>0).');
else
    maxlistofplanestocheck = uint32(find( canplaneobstruct>0).');    
end
nplanestocheck = length(maxlistofplanestocheck);

if nplanestocheck > 0

	ntot = nplanestocheck*npaths;
	
	% Make one long list of planes to check obstruction for. It will be aligned
	% with expanded lists of path startpoints and endpoints in addition to startplanes and
	% endplanes. We need to construct a ref. list because some of the
	% combinations will be thrown away and we must be able to refer back to the
	% rectangular, complete matrix of size [nplanestocheck,npaths].
    
	bigplanelist = maxlistofplanestocheck(:,ones(1,npaths));
	bigplanelist = reshape(bigplanelist,ntot,1);
	
	reftoorigmatrix = uint32([1:npaths*nplanestocheck].');
		
	% Don't check the one or two planes that are involved in each path for obstruction
	% or, their respective rearside planes, if they happen to be thin.
	
	if ~isempty(BIGSTARTPLANES)
        if nplanecols1 == 1
            iv = find(bigplanelist==BIGSTARTPLANES | bigplanelist==rearsideplane(BIGSTARTPLANES));
            bigplanelist(iv) = [];
            BIGSTARTPLANES(iv) = [];
            reftoorigmatrix(iv) = [];
            if ~isempty(BIGENDPLANES)
                BIGENDPLANES(iv,:) = [];
            end
            BIGTOCOORDS(iv,:) = [];
            BIGFROMCOORDS(iv,:) = [];
        else            
            
            iv = find(bigplanelist==BIGSTARTPLANES(:,1) | bigplanelist==BIGSTARTPLANES(:,2) | ...
                      bigplanelist==rearsideplane(BIGSTARTPLANES(:,1)) | bigplanelist==rearsideplane(BIGSTARTPLANES(:,2)));
            bigplanelist(iv) = [];
            BIGSTARTPLANES(iv,:) = [];
            reftoorigmatrix(iv) = [];
            if ~isempty(BIGENDPLANES)
                BIGENDPLANES(iv,:) = [];
            end
            BIGTOCOORDS(iv,:) = [];
            BIGFROMCOORDS(iv,:) = [];
        end
	end
	
	if ~isempty(BIGENDPLANES)
        if nplanecols2 == 1
            iv = find(bigplanelist==BIGENDPLANES | bigplanelist==rearsideplane(BIGENDPLANES));
            bigplanelist(iv) = [];
            if ~isempty(BIGSTARTPLANES)
                BIGSTARTPLANES(iv,:) = [];
            end
            BIGENDPLANES(iv,:) = [];
            reftoorigmatrix(iv) = [];
            BIGTOCOORDS(iv,:) = [];
            BIGFROMCOORDS(iv,:) = [];
        else
            iv = find(bigplanelist==BIGENDPLANES(:,1)           | bigplanelist==BIGENDPLANES(:,2) | ...
                      bigplanelist==rearsideplane(BIGENDPLANES(:,1)) | bigplanelist==rearsideplane(BIGENDPLANES(:,2)));
            bigplanelist(iv) = [];
            if ~isempty(BIGSTARTPLANES)
                BIGSTARTPLANES(iv,:) = [];
            end
            BIGENDPLANES(iv,:) = [];
            reftoorigmatrix(iv) = [];
            BIGTOCOORDS(iv,:) = [];
            BIGFROMCOORDS(iv,:) = [];
        end
	end
    
% [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ESIE2chkISvisible(BIGFROMCOORDS,BIGTOCOORDS,planeeqs(bigplanelist,4),planenvecs(bigplanelist,:),minvals(bigplanelist,:),...
%        maxvals(bigplanelist,:),planecorners(bigplanelist,:),corners,ncornersperplanevec(bigplanelist));
%
% 20050922: The entire ESIE2chkISvisible is pasted in here.

%##########################################################################
%##########################################################################
%##########################################################################
%##########################################################################


% nsou = size(BIGFROMCOORDS,1);
nplanes = length(bigplanelist);

nrec = size(BIGTOCOORDS,1);
if nrec == 1
    onesvec = uint8(ones(nplanes,1));
    BIGTOCOORDS = BIGTOCOORDS(onesvec,:);
    clear onesvec    
end

planenvecsexpanded = planenvecs(bigplanelist,:).';

%----------------------------------------------------------------
% First we check if the IS and R are on the same side of the
% respective planes. Then the path can not pass through the plane.

tempmatrix = corners(planecorners(bigplanelist,1),:);
sseesplanes = sum(  ((BIGFROMCOORDS - tempmatrix).').*planenvecsexpanded   ).';
rseesplanes = sum(  ((BIGTOCOORDS - tempmatrix).').*planenvecsexpanded   ).';
iv1 = uint32(find( (sseesplanes>= (-eps*30) )~=(rseesplanes>= (-eps*30)  ) ));
clear sseesplanes rseesplanes


%--------------------------------------------------------------------------
% Next tests (if there were any planes that had IS and R at different sides)
% 1. Are the vectors IS -> R parallel to the planes?
% 2. Are the reflpoints inside the planes?

if ~isempty(iv1)
% 	npossible = length(iv1);
    % Below tempmatrix is the direction vector
	tempmatrix = BIGTOCOORDS - BIGFROMCOORDS;
    BIGTOCOORDS = [];
    dirveclengths = sqrt(sum(tempmatrix.'.^2)).' + eps*100;
    tempmatrix = tempmatrix./dirveclengths(:,ones(1,3));

	% If test == 0, then the vector S -> R runs along the
	% plane.
	
	iv1 = iv1(sum(   planenvecsexpanded(:,iv1).*(tempmatrix(iv1,:).')   ).'~=0);
    
	if ~isempty(iv1)

		% The last test is if the hitpoint is inside the plane
	
        % Step 1: calculate the hit point using u = the line parameter (with
        % unit meters) from the source towards the receiver.

        udir = ( planeeqs(bigplanelist(iv1),4) - sum( planenvecsexpanded(:,iv1).*(BIGFROMCOORDS(iv1,:).') ).' );        
        
        udir = udir./( sum(   planenvecsexpanded(:,iv1).*(tempmatrix(iv1,:).')   ).');
        
        % Step 2: check that the hit point is at a shorter distance than
        % the receiver point, and that the hit point is not in the opposite
        % direction than the receiver.

        iv2 = find( udir<(dirveclengths(iv1)*1.001) & udir>=0 );
        clear dirveclengths
        if length(iv2) < length(iv1)
            iv1 = iv1(iv2);
            udir = udir(iv2);
            clear iv2
        end
        
        if ~isempty(iv1)
            % Step 3: calculate the actual xyz coordinates of the hit points
            % Now tempmatrix gets the values of the hit points
            tempmatrix = BIGFROMCOORDS(iv1,:) + udir(:,ones(1,3)).*tempmatrix(iv1,:);        
        
            clear udir
            BIGFROMCOORDS = [];
            
%             [hitvec,edgehitvec,cornerhitvec] = EDpoinpla(tempmatrix,iv1,minvals(bigplanelist,:),maxvals(bigplanelist,:),planecorners(bigplanelist,:),corners,ncornersperplanevec(bigplanelist),planenvecsexpanded.');
            [hitvec,edgehitvec,edgehitnumbers,cornerhitvec,cornerhitnumbers] = ...
                EDpoinpla(tempmatrix,iv1,minvals(bigplanelist,:),maxvals(bigplanelist,:),planecorners(bigplanelist,:),corners,ncornersperplanevec(bigplanelist),planenvecsexpanded.');
        
    	   	hitplanes = [];
%             reflpoints = [];
            edgehits = [];
%          edgehitpoints = [];
            cornerhits = [];
%             cornerhitpoints = [];
            if any(any(hitvec)) ~=0
                ivhit = hitvec==1;
			    hitplanes = iv1( ivhit );  
%              reflpoints = tempmatrix(ivhit,:);
            end
            if ~isempty(edgehitvec)
                ivhit = edgehitvec==1;
                edgehits      = iv1( ivhit );
%                 edgehitpoints = tempmatrix(ivhit,:); 
            end
            if ~isempty(cornerhitvec)
                ivhit = cornerhitvec==1;
                cornerhits  = iv1( ivhit );
%                 cornerhitpoints = tempmatrix(ivhit,:);    
            end
        else
    		hitplanes = [];
%             reflpoints = [];
            edgehits = [];
%             edgehitpoints = [];
            cornerhits = [];
%             cornerhitpoints = [];        
        end
	else
		hitplanes = [];
%         reflpoints = [];
        edgehits = [];
%         edgehitpoints = [];
        cornerhits = [];
%         cornerhitpoints = [];
	end
else
	hitplanes = [];   
%     reflpoints = [];
    edgehits = [];
%     edgehitpoints = [];
    cornerhits = [];
%     cornerhitpoints = [];
end

%##########################################################################
%##########################################################################
%##########################################################################
%##########################################################################

    zerosvec1 = zeros(nplanestocheck,npaths);

    obstructionoccurrence = zerosvec1;
	obstructionoccurrence(reftoorigmatrix(hitplanes)) = ones(size(hitplanes));	
	if nplanestocheck > 1
        obstructionoccurrence = sum(obstructionoccurrence);
	end
	nonobstructedpaths = find(obstructionoccurrence==0);
    
    if ~isempty(edgehits)
        edgehitoccurrence = zerosvec1;
        edgehitoccurrence(reftoorigmatrix(edgehits)) = ones(size(edgehits));
        if nplanestocheck > 1
            edgehitoccurrence = sum(edgehitoccurrence);
        end
        edgehits = find(edgehitoccurrence~=0);
    end

    if ~isempty(cornerhits)
        cornerhitoccurrence = zerosvec1;
        cornerhitoccurrence(reftoorigmatrix(cornerhits)) = ones(size(cornerhits));
        if nplanestocheck > 1
            cornerhitoccurrence = sum(cornerhitoccurrence);
        end
        cornerhits = find(cornerhitoccurrence~=0);
    end

    nobstructions = npaths-length(nonobstructedpaths);
else
%     npaths
    nonobstructedpaths = [1:npaths].';
    pause
    nobstructions = 0;
    edgehits = [];
    cornerhits = [];
end
