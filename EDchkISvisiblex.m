function [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDchkISvisiblex(bigplanelist,planeeqs_lastvalue,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec)
% EDchkISvisiblex - Checks if paths from a set of IS to a set of R pass through resp. refl. planes. Extended version
% EDchkISvisiblex checks if the paths between a number of IS and a single R, or a number of R, 
% pass through their respective reflecting planes.
%
% Input parameters:
%   bigplanelist    List, [nIS,1], of the plane numbers in a long list of
%               combinations to check.
%   planeeqs_lastvalue  List [nplanes,1] of the fourth value of the plane
%               equations.
%   planenvecs,minvals, maxvals, planecorners, corners, ncornersperplanevec
%               Data that should have been taken from the corresponding
%               variables in the eddatafile.
%               NB!! The matrices planeeqs, minvals, maxvals
%               have been rearranged so that they have nIS rows, and each
%               row contain the data for one specific plane: the reflecting
%               plane for the IS.
%   BIGFROMCOORDSSHORTLIST,REFTOFROMSHORTLIST,BIGTOCOORDSSHORTLIST,REFTOTOSHORTLIST
%               (GLOBAL)
%
% Output parameters:
%	 hitplanes		List, [nhits,1], of the planes that have been hit. The values
%	                in this list refer to the index numbers of the input list ISlist
%    hitpoints      List, [nhits,3] of the hitpoint coordinates, for the planes that have been hit
%    edgehits       List, [nedgehits,1] of the planes that were hit right at an edge.
%                   These combinations were marked as hit in the list hitplanes, 
%                   but the extra information in edgehits can possibly be used.
%    edgehitpoints  List, [nedgehits,3] of the hitpoint coordinates for the
%                   planes that were hit at an edge.
%    cornerhits     List, [ncornerhits,1] of the planes that were hit right at a corner.
%                   These combinations were marked as hit in the list hitplanes, 
%                   but the extra information in edgehits can possibly be used.
%    cornerhitpoints  List, [ncornerhits,3] of the hitpoint coordinates for the
%                   planes that were hit at a corner.
%
% Uses the function EDpoinplax.
%
% Peter Svensson (peter.svensson@ntnu.no) 29 Nov. 2016
% 
% [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDchkISvisiblex(bigplanelist,planeeqs_lastvalue,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec)

% 16 Jan. 2005 Functioning version
% 27 Nov. 2017 Copied from ESIE2toolbox
% 29 Nov. 2017 Called ESIE2poinplax - changed to ED

global BIGFROMCOORDSSHORTLIST REFTOFROMSHORTLIST BIGTOCOORDSSHORTLIST REFTOTOSHORTLIST

if isempty(bigplanelist)
    hitplanes = [];
    hitpoints = [];
    edgehits = [];
    edgehitpoints = [];
    cornerhits = [];
    cornerhitpoints = [];
    return
end

% nsou = size(REFTOFROMSHORTLIST,1);
npaths = size(bigplanelist,1);

if size(planeeqs_lastvalue,2) > 1
    error('Fix the call to EDchkISvisiblex!!!')    
end

nrec = size(REFTOTOSHORTLIST,1);

planenvecs = planenvecs.';

%----------------------------------------------------------------
% First we check if the IS and R are on the same side of the
% respective planes. Then the path can not pass through the plane.

if nrec == 1
    iv1 = uint32(find( (sum(  ((BIGFROMCOORDSSHORTLIST(REFTOFROMSHORTLIST,:) - corners(planecorners(bigplanelist,1),:)).').*planenvecs(:,bigplanelist)   ).' >= (-eps*30) )~=( sum(  ((BIGTOCOORDSSHORTLIST(ones(npaths,1),:) - corners(planecorners(bigplanelist,1),:)).').*planenvecs(:,bigplanelist)   ).'>= (-eps*30)  ) ));
else
    iv1 = uint32(find( (sum(  ((BIGFROMCOORDSSHORTLIST(REFTOFROMSHORTLIST,:) - corners(planecorners(bigplanelist,1),:)).').*planenvecs(:,bigplanelist)   ).' >= (-eps*30) )~=( sum(  ((BIGTOCOORDSSHORTLIST(REFTOTOSHORTLIST,:) - corners(planecorners(bigplanelist,1),:)).').*planenvecs(:,bigplanelist)   ).'>= (-eps*30)  ) ));    
end

%--------------------------------------------------------------------------
% Next tests (if there were any planes that had IS and R at different sides)
% 1. Are the vectors IS -> R parallel to the planes?
% 2. Are the hitpoints inside the planes?

if ~isempty(iv1)
% 	npossible = length(iv1);
    % Below tempmatrix is the direction vector

    if nrec == 1
        tempmatrix = BIGTOCOORDSSHORTLIST(ones(npaths,1),:) - BIGFROMCOORDSSHORTLIST(REFTOFROMSHORTLIST,:);
    else
        tempmatrix = BIGTOCOORDSSHORTLIST(REFTOTOSHORTLIST,:) - BIGFROMCOORDSSHORTLIST(REFTOFROMSHORTLIST,:);        
    end
    dirveclengths = sqrt(sum(tempmatrix.'.^2)).' + eps*100;
    tempmatrix = tempmatrix./dirveclengths(:,ones(1,3));

	% If test == 0, then the vector S -> R runs along the
	% plane.
	
	iv1 = iv1(sum(   planenvecs(:,bigplanelist(iv1)).*(tempmatrix(iv1,:).')   ).'~=0);
    
	if ~isempty(iv1)

		% The last test is if the hitpoint is inside the plane

        % Step 1: calculate the hit point using u = the line parameter (with
        % unit meters) from the source towards the receiver.
        ISlist = BIGFROMCOORDSSHORTLIST(REFTOFROMSHORTLIST(iv1),:);
        udir = ( planeeqs_lastvalue(bigplanelist(iv1)) - sum( planenvecs(:,bigplanelist(iv1)).*(ISlist.') ).' );        
        clear planeeqs_lastvalue
        udir = udir./( sum(   planenvecs(:,bigplanelist(iv1)).*(tempmatrix(iv1,:).')   ).');
        
        % Step 2: check that the hit point is at a shorter distance than
        % the receiver point, and that the hit point is not in the opposite
        % direction than the receiver.
        iv2 = find(udir < -eps*1e4);
        if ~isempty(iv2)
            disp('Lowest value (should always be > 0)')
            min(udir(iv2))
             error('WARNING!!! Check EDchkISvisiblex; some hit point is in the wrong direction from the receiver. Check the plane warping tolerances etc.')
        end
        
        iv2 = find(   (udir - dirveclengths(iv1)) > eps*1e4);
        if ~isempty(iv2)
            disp('Lowest value (should always be > eps*1e4):')
            min(udir(iv2) - dirveclengths(iv1(iv2)))
            disp(['Problematic hit point is related to plane ',int2str(bigplanelist(iv1(iv2(1)))),' and IS w coords: '])
            ISlist(iv2(1),:)
             error('WARNING!!! Check EDchkISvisiblex; some hit point was at a shorter distance than the IS-receiver distance. Check the plane warping tolerances etc.')
        end
        
        
        % Step 3: calculate the actual xyz coordinates of the hit points
        % Now tempmatrix gets the values of the hit points
        tempmatrix = ISlist + udir(:,ones(1,3)).*tempmatrix(iv1,:);
        
        clear ISlist udir
        [hitvec,edgehitvec,cornerhitvec] = EDpoinplax(bigplanelist,tempmatrix,iv1,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs.');
        
	   	hitplanes = [];
        hitpoints = [];
        edgehits = [];
        edgehitpoints = [];
        cornerhits = [];
        cornerhitpoints = [];
        if any(any(hitvec)) ~=0
            ivhit = find(hitvec==1);
			hitplanes = iv1( ivhit );  
            hitpoints = tempmatrix(ivhit,:);
        end
        if ~isempty(edgehitvec)
            ivhit = find(edgehitvec==1);
            edgehits      = iv1( ivhit );
            edgehitpoints = tempmatrix(ivhit,:); 
        end
        if ~isempty(cornerhitvec)
            ivhit = find(cornerhitvec==1);
            cornerhits  = iv1( ivhit );
            cornerhitpoints = tempmatrix(ivhit,:);
        end
        
	else
		hitplanes = [];
        hitpoints = [];
        edgehits = [];
        edgehitpoints = [];
        cornerhits = [];
        cornerhitpoints = [];
	end
else
	hitplanes = [];   
    hitpoints = [];
    edgehits = [];
    edgehitpoints = [];
    cornerhits = [];
    cornerhitpoints = [];
end
