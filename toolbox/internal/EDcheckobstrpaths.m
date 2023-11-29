function [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane)
% EDcheckobstrpaths - Checks obstruction for a list of paths, or check if specular reflections are valid.
% Checks if the paths from N starting points (fromcoords) to
% N ending points (tocoords) pass through any of M planes (those marked by 1 in
% canplaneobstruct), i.e., if they are obstructed.
% 
% Instead of testing for N to N points, it is possible to specify
% 1 starting point and N ending points or the other way around.
% Also, it is possible to specify, for each individual starting or ending point,
% one or two planes that the respective point is directly at,
% because then, those planes will be excluded from the list of planes to check for.
%
% Input parameters:
%   fromcoords          Matrix, [N,3] or [1,3], with the coordinates of the
%                       N (or 1) starting points of the N paths to check.
%   tocoords            Matrix, [N,3] or [1,3], with the coordinates of the
%                       N (or 1) ending points of the N paths to check.
%   startplanes         Matrix, [N,1] or [N,2] or [0,0], with the plane
%                       numbers that the starting points are positioned directly at.
%                       The size is determined by:
%                       [N,1] fromcoords are specular reflection points
%                       [N,2] fromcoords are edge points
%                       [0,0] fromcoords are image sources or receiver
%                       positions.
%   endplanes           Matrix with the plane numbers that the ending points are
%                       positioned directly at.
%   canplaneobstruct,planeseesplane,planeeqs,planenvecs,minvals,maxvals,..
%   planecorners,corners,ncornersperplanevec,rearsideplane
%                       Data that should have been passed on from the
%                       eddata struct. 
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
% Uses function EDchkISvisiblex
%
% Peter Svensson (peter.svensson@ntnu.no) 29 Nov. 2017
%
% [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
%    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);

% 21 June 2006 Functioning version
% 27 Nov. 2017 Copied to EDtoolbox
% 28 Nov. 2017 Cleaned code a bit.
% 29 Nov. 2017 Fixed a problem with mismatch between planesareseen (logical) 
%              and canplaneobstruct (int8)

global BIGFROMCOORDSSHORTLIST REFTOFROMSHORTLIST BIGTOCOORDSSHORTLIST REFTOTOSHORTLIST

nstartpoints = size(fromcoords,1);
nendpoints = size(tocoords,1);
npaths = max(max([nstartpoints nendpoints]));
nplanes = length(canplaneobstruct);

% The startplanes and endplanes can either be one-column or two-column
% lists depending on if paths are between specular reflections or edges.

[nplanes1,nplanecols1] = size(startplanes);
[nplanes2,nplanecols2] = size(endplanes);

% If startplanes and/or endplanes have been specified, then we should only
% need to check planes that are seen by the start- or end-plane. 

if ~isempty(startplanes) && ~isempty(endplanes)
    if nplanecols1 == 1 && nplanecols2 == 1
        planesareseen = ( sum( planeseesplane(endplanes,:) > 0 ) + sum( planeseesplane(startplanes,:) > 0 ) ) > 0;
    elseif nplanecols1 == 2 && nplanecols2 == 2
        planesareseen = ( sum( planeseesplane(endplanes(:,1),:) > 0 ) + sum( planeseesplane(endplanes(:,2),:) > 0 ) + ...
            sum( planeseesplane(startplanes(:,1),:) > 0 ) + sum( planeseesplane(startplanes(:,2),:) > 0 ) ) > 0;
    elseif nplanecols1 == 2 && nplanecols2 == 1
        planesareseen = ( sum( planeseesplane(endplanes,:) > 0 ) + ...
            sum( planeseesplane(startplanes(:,1),:) > 0 ) + sum( planeseesplane(startplanes(:,2),:) > 0 ) ) > 0;
    elseif nplanecols1 == 1 && nplanecols2 == 2
        planesareseen = ( sum( planeseesplane(endplanes(:,1),:) > 0 ) + sum( planeseesplane(endplanes(:,2),:) > 0 ) + ...
           sum( planeseesplane(startplanes,:) > 0 ) ) > 0;
    end
else
    if ~isempty(startplanes)
        if nplanes1 > 1
            planesareseen = sum( planeseesplane(startplanes,:) > 0 ) > 0;       
        else
            planesareseen = planeseesplane(startplanes,:) > 0;
        end
    elseif ~isempty(endplanes)
        if nplanes2 > 1
            planesareseen = sum( planeseesplane(endplanes,:) > 0 ) > 0;      
        else
            planesareseen = planeseesplane(endplanes,:) > 0;                  
        end
    else
        planesareseen = ones(1,nplanes);
    end
end

% In addition, we only need to check planes for which canplaneobstruct = 1.
% planesareseen has logical values
% canplaneobstruct has int8 values

if nplanes <= 255
    maxlistofplanestocheck = uint8(find((planesareseen & canplaneobstruct==1)>0).');
elseif nplanes <= 65535
    maxlistofplanestocheck = uint16(find((planesareseen & canplaneobstruct==1)>0).');    
else
    maxlistofplanestocheck = uint32(find((planesareseen & canplaneobstruct==1)>0).');
end

nplanestocheck = length(maxlistofplanestocheck);

if nplanestocheck > 0
	onesvec1 = ones(1,nplanestocheck);
	onesvec2 = ones(npaths*nplanestocheck,1);
	ntot = nplanestocheck*npaths;
    
    % Make one long list of planes to check obstruction for. It will be aligned
	% with expanded lists of path startpoints and endpoints in addition to startplanes and
	% endplanes. We need to construct a ref. list because some of the
	% combinations will be thrown away and we must be able to refer back to the
	% rectangular, complete matrix of size [nplanestocheck,npaths].

    bigplanelist = maxlistofplanestocheck(:,ones(1,npaths));
	bigplanelist = reshape(bigplanelist,ntot,1);
	
	reftoorigmatrix = uint32([1:npaths*nplanestocheck].');
    
    if nstartpoints == 1 && nendpoints ~= 1       
        speedyalgorithm = 1;
    else
        speedyalgorithm = 0;        
    end

    % Make one long matrix of path endpoints and path endplanes
	    
	if nendpoints > 1
        if speedyalgorithm == 0
%             bigtocoords = reshape(repmat(tocoords.',[nplanestocheck,1]),3,ntot).';
            BIGTOCOORDSSHORTLIST = tocoords;
%             REFTOTOSHORTLIST = uint32([1:npaths*nplanestocheck].');
            if npaths < 65536
                REFTOTOSHORTLIST = uint16(1:npaths);
            else
                REFTOTOSHORTLIST = uint32(1:npaths);                
            end
          REFTOTOSHORTLIST = reshape(REFTOTOSHORTLIST(ones(nplanestocheck,1),:),ntot,1);
        end
        
        if isempty(endplanes)
            bigendplanes = [];
        else
            if nplanecols2 == 1
                bigendplanes = endplanes(:,onesvec1);
                bigendplanes = reshape(bigendplanes.',ntot,1);
            else               
                bigendplanes = reshape(repmat(endplanes.',[nplanestocheck,1]),2,ntot).';
            end
        end
	else
        if speedyalgorithm == 0
%             bigtocoords = tocoords(onesvec2,:);
            BIGTOCOORDSSHORTLIST = tocoords;
            REFTOTOSHORTLIST = int8(onesvec2);
        end
        if isempty(endplanes)
            bigendplanes = [];
        else
            if length(endplanes) > 2
                error('ERROR: Some strange error here!')    
            end
            if nplanecols2 == 1
                bigendplanes = endplanes(onesvec2,:);
            else
                bigendplanes = [endplanes(1) endplanes(2)];
                bigendplanes = bigendplanes(onesvec2,:);
            end
        end
	end
	
	% Make one long matrix of path startplanes and startpoints
	
	if nstartpoints > 1
        if speedyalgorithm == 0
%             bigfromcoords = reshape(repmat(fromcoords.',[nplanestocheck,1]),3,ntot).';
            BIGFROMCOORDSSHORTLIST = fromcoords;
            if npaths < 65536
                REFTOFROMSHORTLIST = uint16(1:npaths);
            else
                REFTOFROMSHORTLIST = uint32(1:npaths);                
            end
          REFTOFROMSHORTLIST = reshape(REFTOFROMSHORTLIST(ones(nplanestocheck,1),:),ntot,1);
        end
        if isempty(startplanes)
            bigstartplanes = [];
        else
            if nplanecols1 == 1
                bigstartplanes = startplanes(:,onesvec1);
                bigstartplanes = reshape(bigstartplanes.',ntot,1);
            else
                bigstartplanes = reshape(repmat(startplanes.',[nplanestocheck,1]),2,ntot).';
            end  
        end    
	else
        if speedyalgorithm == 0
%             bigfromcoords = fromcoords(onesvec2,:);
            BIGFROMCOORDSSHORTLIST = fromcoords;
            REFTOFROMSHORTLIST = int8(onesvec2);
        end
        if isempty(startplanes)
            bigstartplanes = [];
        else
            if nplanecols1 == 1
                bigstartplanes = startplanes(onesvec2,:);
            else
                bigstartplanes = [startplanes(1) startplanes(2)];
                bigstartplanes = bigstartplanes(onesvec2,:);
            end
        end
	end

    if speedyalgorithm == 0
        clear onesvec2
    end
    
    % Now we can clear all combinations where neither the startplane or the
    % endplane sees the obstructing plane!
    %
    % After a number of combinations have been cleared, the fromcoords and
    % tocoords can be expanded and then pruned.
    %
    % This is the speedy algorithm...

    if nstartpoints == 1 && nendpoints ~= 1      
        if size(bigendplanes,2) == 1
%             refto_planeseesplane = uint32(double(bigendplanes) + (double(bigplanelist)-1)*nplanes);
%             iv = uint32(find(planeseesplane(refto_planeseesplane)==0));
%             clear refto_planeseesplane
            iv = uint32(find(planeseesplane(uint32(double(bigendplanes) + (double(bigplanelist)-1)*nplanes))==0));
            bigplanelist(iv) = [];
            bigendplanes(iv,:) = [];
            reftoorigmatrix(iv) = [];
            
%             bigfromcoords = fromcoords(onesvec2(1:length(reftoorigmatrix)),:);
            BIGFROMCOORDSSHORTLIST = fromcoords;
            REFTOFROMSHORTLIST = int8(onesvec2(1:length(reftoorigmatrix)));
            clear onesvec2
            
%             bigtocoords = reshape(repmat(tocoords.',[nplanestocheck,1]),3,ntot).';     
                       
%             bigtocoords(iv,:) = [];
            BIGTOCOORDSSHORTLIST = tocoords;
%             npaths
%             nplanestocheck
%             REFTOTOSHORTLIST = uint32([1:npaths*nplanestocheck].');
            if npaths < 65536
                REFTOTOSHORTLIST = uint16(1:npaths);
            else
                REFTOTOSHORTLIST = uint32(1:npaths);                
            end
            REFTOTOSHORTLIST = reshape(REFTOTOSHORTLIST(ones(nplanestocheck,1),:),ntot,1);
            REFTOTOSHORTLIST(iv) = [];
            clear iv
        else
            ref1to_planeseesplane = uint32(double(bigendplanes,1) + (double(bigplanelist)-1)*nplanes);
            ref2to_planeseesplane = uint32(double(bigendplanes,2) + (double(bigplanelist)-1)*nplanes);
            iv = uint32(find(planeseesplane(ref1to_planeseesplane)==0 && planeseesplane(ref2to_planeseesplane)==0));
            clear refto_planeseesplane
            bigplanelist(iv) = [];
            bigendplanes(iv,:) = [];
            reftoorigmatrix(iv) = [];
            
%             bigfromcoords = fromcoords(onesvec2(1:length(reftoorigmatrix)),:);
            BIGFROMCOORDSSHORTLIST = fromcoords;
            REFTOFROMSHORTLIST = int8(onesvec2(1:length(reftoorigmatrix)));
            clear onesvec2
            
%             bigtocoords = reshape(repmat(tocoords.',[nplanestocheck,1]),3,ntot).';            
%             bigtocoords(iv,:) = [];
            BIGTOCOORDSSHORTLIST = tocoords;
%             npaths
%             nplanestocheck
%             REFTOTOSHORTLIST = uint32([1:npaths*nplanestocheck].');
            if npaths < 65536
                REFTOTOSHORTLIST = uint16(1:npaths);
            else
                REFTOTOSHORTLIST = uint32(1:npaths);                
            end
            REFTOTOSHORTLIST = reshape(REFTOTOSHORTLIST(ones(nplanestocheck,1),:),ntot,1);
            REFTOTOSHORTLIST(iv) = [];
            clear iv
        end
    end
    if nstartpoints ~= 1 && nendpoints ~= 1       
        
    end
    
	% Don't check the one or two planes that are involved in each path for obstruction
	% or, their respective rearside planes, if they happen to be thin.
	
	if ~isempty(bigstartplanes)
        if nplanecols1 == 1
            iv = find(bigplanelist==bigstartplanes | bigplanelist==rearsideplane(bigstartplanes));
            clear bigstartplanes
            bigplanelist(iv) = [];
%             bigstartplanes(iv) = [];
            reftoorigmatrix(iv) = [];
            %    if nendpoints > 1 && ~isempty(bigendplanes),
            if ~isempty(bigendplanes)
                bigendplanes(iv,:) = [];
            end
%             bigtocoords(iv,:) = [];
            REFTOTOSHORTLIST(iv) = [];
            REFTOFROMSHORTLIST(iv) = [];
%             bigfromcoords(iv,:) = [];
        else
            iv = find(bigplanelist==bigstartplanes(:,1) | bigplanelist==bigstartplanes(:,2) | ...
                      bigplanelist==rearsideplane(bigstartplanes(:,1)) | bigplanelist==rearsideplane(bigstartplanes(:,2)));
            clear bigstartplanes
            bigplanelist(iv) = [];
%             bigstartplanes(iv,:) = [];
            reftoorigmatrix(iv) = [];
            %    if nendpoints > 1 && ~isempty(bigendplanes),
            if ~isempty(bigendplanes)
                bigendplanes(iv,:) = [];
            end
%             bigtocoords(iv,:) = [];
            REFTOTOSHORTLIST(iv) = [];
            REFTOFROMSHORTLIST(iv) = [];
%             bigfromcoords(iv,:) = [];
        end
	end
	
	if ~isempty(bigendplanes)
        if nplanecols2 == 1
            
            iv = find(bigplanelist==bigendplanes | bigplanelist==rearsideplane(bigendplanes));
            
            clear bigendplanes
            bigplanelist(iv) = [];
		%     if nstartpoints > 1 && ~isempty(bigstartplanes),
%             if ~isempty(bigstartplanes),
%                 bigstartplanes(iv,:) = [];
%             end
%             bigendplanes(iv,:) = [];
            reftoorigmatrix(iv) = [];
%             bigtocoords(iv,:) = [];
            REFTOTOSHORTLIST(iv) = [];
            REFTOFROMSHORTLIST(iv) = [];
%             bigfromcoords(iv,:) = [];
        else
            iv = find(bigplanelist==bigendplanes(:,1)           | bigplanelist==bigendplanes(:,2) | ...
                      bigplanelist==rearsideplane(bigendplanes(:,1)) | bigplanelist==rearsideplane(bigendplanes(:,2)));
            clear bigendplanes
            bigplanelist(iv) = [];
%             if ~isempty(bigstartplanes),
%                 bigstartplanes(iv,:) = [];
%             end
%             bigendplanes(iv,:) = [];
            reftoorigmatrix(iv) = [];
%             bigtocoords(iv,:) = [];
            REFTOTOSHORTLIST(iv) = [];
            REFTOFROMSHORTLIST(iv) = [];
%             bigfromcoords(iv,:) = [];
        end
	end
        

% [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = ESIE2chkISvisiblex(bigplanelist,planeeqs(:,4),planenvecs,minvals,...
%         maxvals,planecorners,corners,ncornersperplanevec);
[hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDchkISvisiblex(bigplanelist,planeeqs(:,4),planenvecs,minvals,...
        maxvals,planecorners,corners,ncornersperplanevec);

   
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
    nonobstructedpaths = [1:npaths].';   
    nobstructions = 0;
    edgehits = [];
    cornerhits = [];
end
