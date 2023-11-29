function [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_pointtoedge(fromcoordsshortlist,...
    reftoshortlist,tocoords,canplaneobstruct,planedata)
% EDcheckobstr_pointtoedge - Checks obstruction for a list of paths from points to edges.
% Checks if the paths from N starting points (reftoshortlist -> fromcoordsshortlist) to
% N ending points on edges (tocoords) pass through any of M planes (those marked by 1 in
% canplaneobstruct), i.e., if they are obstructed.
%
% Input parameters:
%   fromcoordsshortlist Matrix with the unique coordinates of the N
%                       starting points.
%   reftoshortlist      List, [N,1], with the N indices to fromcoordsshortlist.
%   tocoords            Matrix, [N,3], with the coordinates of the
%                       N ending points (on edges) of the N paths to check.
%   canplaneobstruct    Special version (not same as in planedata struct)
%   planedata           struct
%
% Output parameters:
%   nonobstructedpaths  A list, [N_nobstructions,1] of the indices of paths that are not
%                       obstructed. The indices refer to the input
%                       lists fromcoords, tocoords, startplanes,
%                       endplanes
%   nobstructions       The number of paths (of the N input paths) that
%                       were not obstructed.
%   edgehits            A list, [N_edgehits,1] of the indicies of paths that hit a
%                       plane right at an edge.
%   cornerhits          A list, [N_cornerhits,1] of the indicies of paths that hit a
%                       plane right at a corner.   
%
% Uses function EDpoinpla
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2017
%
% [nonobstructedpaths,nobstructions,edgehits,cornerhits] =
% EDcheckobstr_pointtoedge(fromcoordsshortlist,reftoshortlist,tocoords,...
% canplaneobstruct,planedata);

% 22 Sep. 2005 Functioning version
% 27 Nov. 2017 Copied from ESIE2toolbox

npaths = size(tocoords,1);
nplanes = size(canplaneobstruct,2);

% We only need to check planes for which canplaneobstruct = 1.

if nplanes < 256
    maxlistofplanestocheck = uint8( find( sum(canplaneobstruct) ).' );
elseif nplanes < 65536
    maxlistofplanestocheck = uint16( find( sum(canplaneobstruct) ).' );
else
    maxlistofplanestocheck = uint32( find( sum(canplaneobstruct) ).' );
end
    
nplanestocheck = length(maxlistofplanestocheck);

if nplanestocheck > 0
	ntot = nplanestocheck*npaths;
    
%     nmaxcorners = double(max(planedata.ncornersperplanevec(maxlistofplanestocheck)));
%     planecorners = planedata.planecorners(:,1:nmaxcorners);
    
    % For the obstruction test, we create expanded matrices/lists where the
    % plane-number is the fastest varying variable. We calculate
    % the direction vectors, from point to point, before the expanding.

    dirvec = tocoords - fromcoordsshortlist(reftoshortlist,:);
    clear tocoords
    dirveclengths = sqrt(sum(dirvec.'.^2)).' + eps*100;
    dirvec = dirvec./dirveclengths(:,ones(1,3));

    % Create the expanded matrices, by duplicating each path with the
    % number of planes that need to be checked obstruction against.
    % Immediately select only those individual combinations where
    % canplaneobstruct = 1. This selection is done both for bigplanelist
    % and for the reftoshortlist, and for the reftodirvec (which points to
    % the shorter dirvec matrix).

    iv1 = (find(canplaneobstruct(:,maxlistofplanestocheck)));
    
    bigplanelist = maxlistofplanestocheck(ceil((iv1)/npaths));
    if ntot < 65536
        iv1 = uint16(iv1);
    elseif ntot < 4.29e9
        iv1 = uint32(iv1);    
    end	

    reftoshortlist = reftoshortlist(:,ones(1,nplanestocheck));
    reftoshortlist = reftoshortlist(iv1);

    if npaths*nplanestocheck < 65536
        reftodirvec = uint16([1:npaths].');
    else
        reftodirvec = uint32([1:npaths].');
    end        
    reftodirvec = reftodirvec(:,ones(1,nplanestocheck));
    reftodirvec = reftodirvec(iv1);

    %--------------------------------------------------------------------------
    % Tests
    % (Earlier, it was first tested if the IS and R were on the same side of a plane. This is checked at another place).
    % (Earlier, test 2 was: Are the vectors IS -> R parallel to the planes?
    % This shold not be needed now).
    %
    % Are the hitpoints inside the planes?
    % 1. Calculate the hit points' coordinates. 
    % 2. Check that the hit point is not further away than the receiver
    %    point, and that the hit point is not in the wrong direction!
    % 3. Check if the hit points are inside the plane that is checked.

    % Step 1
    udir_bigplanenvecs = ( planedata.planeeqs(bigplanelist,4) - sum( planedata.planeeqs(bigplanelist,1:3).'.*(fromcoordsshortlist(reftoshortlist,:).') ).' )./( sum(   planedata.planeeqs(bigplanelist,1:3).'.*(dirvec(reftodirvec,:).')   ).');

    % Step 2
    iv2 = find( udir_bigplanenvecs<dirveclengths(reftodirvec) & udir_bigplanenvecs>0 );
    clear dirveclengths
    iv1 = iv1(iv2);
    bigplanelist = bigplanelist(iv2);
    reftoshortlist = reftoshortlist(iv2);
    reftodirvec = reftodirvec(iv2);
    udir_bigplanenvecs = udir_bigplanenvecs(iv2,:);
    clear iv2

    if ~isempty(iv1)
        % Step 3
        xhit_list = fromcoordsshortlist(reftoshortlist,:) + udir_bigplanenvecs(:,ones(1,3)).*dirvec(reftodirvec,:);
    
        clear udir_bigplanenvecs dirvec reftoshortlist


        [hitvec,edgehit,cornerhit] = EDpoinpla(xhit_list,bigplanelist,planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.planeeqs(:,1:3));
        clear xhit_list
    
        % We recycle canplaneobstruct since we need to compile how many hits
        % that occurred.
        canplaneobstruct = zeros(size(canplaneobstruct(:,maxlistofplanestocheck)));
        ivhits = find(hitvec);
        if ~isempty(ivhits)
%             indexvec = iv1(ivhits);
            canplaneobstruct(iv1(ivhits)) = canplaneobstruct(iv1(ivhits))+1;
        end    
        edgehits = [];
        ivhits = find(edgehit);
        if ~isempty(ivhits)
            disp('WARNING: Edge hit. This is treated as an obstruction.')
%             indexvec = iv1(ivhits);
            canplaneobstruct(iv1(ivhits)) = canplaneobstruct(iv1(ivhits))+1;
        end

        cornerhits = [];
        ivhits = find(cornerhit);
        if ~isempty(ivhits)
            disp('WARNING: Corner hit. This is treated as an obstruction.')
%             indexvec = iv1(ivhits);
            canplaneobstruct(iv1(ivhits)) = canplaneobstruct(iv1(ivhits))+1;
        end

        % Here was a bug: erroneous summing when there was a single plane that
        % could obstruct (PS 041127)
        if size(canplaneobstruct,2) > 1
            canplaneobstruct = 1-sign(sum(canplaneobstruct.'));
        else
            canplaneobstruct = 1-canplaneobstruct;            
        end
        nonobstructedpaths = find(canplaneobstruct);

        nobstructions = npaths-length(nonobstructedpaths);
    else
        nonobstructedpaths = [1:npaths].';   
        nobstructions = 0;
        edgehits = [];
        cornerhits = [];
    end
else
    nonobstructedpaths = [1:npaths].';   
    nobstructions = 0;
    edgehits = [];
    cornerhits = [];
end
