function [hitvec,edgehit,cornerhit] = EDpoinplax(bigplanelist,xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs)
% EDpoinplax - Detects if one or more points are inside a number of finite planes. Extended version.
% If one point is given as input, it will be checked against all
% planes in a list of planes. If N points are given as inputs, a list of
% planes should have the N planes, and each point will be checked against
% its corresponding plane.
%
% Input parameters:
% NB!!! This version is a special version, so bigplanelist and xpoints
% actually have different sizes! planelist and xpoints have the same sizes
% however.
%   bigplanelist    List, [N,1], of the plane numbers in all the
%                   N points to check.
%   xpoints         Matrix, [N,3], of coordinates for the N points to check
%   planelist       List, [nplanes,1], of planes to check the N points
%                   against. If N~= 1, then nplanes must be equal to N.
%                   NB! 
%   minvals, maxvals, planecorners, corners, ncornersperplanevec, planenvecs
%                   Data that should have been taken from the corresponding
%                   variables in the eddatafile.
%                   NB!! All of these matrices except corners
%                   have been rearranged so that they have N rows, and each
%                   row contain the data for one specific plane, the one
%                   that the inside-check should be done for.
%
% Output parameters:
%   hitvec			List, [N,1], with 1 or 0 indicating whether a point is
%                   inside or outside the plane given by planelist.
%   edgehit         List, [N,1], with 1 or 0 indicating if a hit was right at the
%                   edge of a plane. These hits were not marked in hitvec.
%   cornerhit       List, [N,1], with 1 or 0 indicating if a hit was right at the
%                   corner of a plane. These hits were not marked in hitvec.
%
% Uses no special functions
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2017
% 
% [hitvec,edgehit] = EDpoinplax(bigplanelist,xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs);

% 4 Febr. 2010 First version in the ESIE toolbox
% 1 May 2017   Fixed one bug: the ray casting algorithm didn't work for
%              rays that hit a corner.
% 27 Nov. 2017 Copied from ESIE2toolbox

% global SHOWTEXT
SHOWTEXT = 1;

geomacc = 1e-12;

npoints = size(xpoints,1);
nplanestotest = length(planelist);
if npoints == 1
	xpoints = xpoints(ones(nplanestotest,1),:);
end

planecorners = [planecorners planecorners(:,1)];

%------------------------------------------------------------
% First test: are the points inside the cubic boxes?
%
% NB! The values in possibleones tell which entries in the lists xpoint
%     and planelist that passed the first test.

if nplanestotest <= 65535
    possibleones = uint16(1:nplanestotest);
else
    possibleones = uint32(1:nplanestotest);    
end
possibleones = possibleones(:);

iv = uint32(find(xpoints(:,1) > maxvals(bigplanelist(planelist),1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,1) < minvals(bigplanelist(planelist(possibleones)),1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) > maxvals(bigplanelist(planelist(possibleones)),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) < minvals(bigplanelist(planelist(possibleones)),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) > maxvals(bigplanelist(planelist(possibleones)),3)));
clear maxvals
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) < minvals(bigplanelist(planelist(possibleones)),3)));
possibleones(iv) = [];
clear minvals iv

hitvec = zeros(nplanestotest,1,'uint8');
hitvec(possibleones) = 1;

nposs = length(possibleones);

if SHOWTEXT >= 4
    disp(['         Of the ',int2str(npoints),' points,'])
    disp(['         ',int2str(length(possibleones)),' survived the cube test:'])      
end

edgehit   = zeros(nposs,1,'uint32');
cornerhit = zeros(nposs,1,'uint32');


%------------------------------------------------------------
% Second test: project onto two dimensions.
% Start by finding which dimension of the planes that has the strongest
% normal vector component. That dimension should be tossed.
%
% Easiest way to handle: make three subsets:
%   Combos that should be projected onto xy
%   Combos that should be projected onto xz
%   Combos that should be projected onto yz

if nposs>0   

    A = abs((planenvecs(bigplanelist(planelist(possibleones)),:)));
    maxA = max(A.').';
    markwhereismax = (A == maxA(:,[1 1 1]));
    iv = find(markwhereismax(:,1).*markwhereismax(:,2)~=0);
    if ~isempty(iv)
        markwhereismax(iv,2) = 0;
    end
    iv = find(markwhereismax(:,1).*markwhereismax(:,3)~=0);
    if ~isempty(iv)
        markwhereismax(iv,1) = 0;
    end
    iv = find(markwhereismax(:,2).*markwhereismax(:,3)~=0);
    if ~isempty(iv)
        markwhereismax(iv,2) = 0;
    end
    colno = [1 2 3];
    finalcolno = markwhereismax.*colno(ones(nposs,1),:);
    finalcolno = sum(finalcolno.').';

    yzsubset = find(finalcolno==1);
    xzsubset = find(finalcolno==2);
    xysubset = find(finalcolno==3);
    
    %---------------------------------------------
    % First the xysubset
    %
    % Create a ray that starts in xpoint and extends parallel to the x-axis
    % in the positive x-direction, that is:
    %       y = xpoints(:,2):
    %       xstart = xpoints(:,1);
    
    if ~isempty(xysubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(xysubset)),' xy-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(bigplanelist(planelist(possibleones(xysubset))));    
        edgenumbers = unique(numberofedgestocheck);

        yray = xpoints(possibleones(xysubset),2);
        xstart = xpoints(possibleones(xysubset),1);
        edgecrossings = zeros(size(xysubset));

        % Use a parametric representation for each edge:
        % x_edge = x_1 + t*(x_2 - x_1)
        % y_edge = y_1 + t*(y_2 - y_1)
        % Find t by setting y_ray = y_edge

        for ii = 1:max(edgenumbers)

            y1 = corners(planecorners(bigplanelist(planelist(possibleones(xysubset))),ii),2);
            y2 = corners(planecorners(bigplanelist(planelist(possibleones(xysubset))),ii+1),2);
            tedgecrossing = (yray - y1)./(y2-y1);

            x1 = corners(planecorners(bigplanelist(planelist(possibleones(xysubset))),ii),1);
            x2 = corners(planecorners(bigplanelist(planelist(possibleones(xysubset))),ii+1),1);

            xedge = x1 + tedgecrossing.*(x2-x1);

            edgecrossings = edgecrossings + (tedgecrossing > 0 & tedgecrossing < 1 & xedge > xstart).*(ii <= numberofedgestocheck);
            
%             closeones = find(abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc );
            closeones = find( (abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc) & xedge > xstart );
            if ~isempty(closeones)
%                disp(['               WARNING:  ',int2str(length(closeones)),' very close hits']) 
                edgecrossings(closeones) = edgecrossings(closeones) + ( ( (y1(closeones)>=yray) + (y2(closeones)>=yray) )== 1);                              
            end
        end
                 
        hitvec(possibleones(xysubset)) = (edgecrossings>0);

        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings>0))),' survived the xyplane projections test:'])  
        end
        
    end
        
    %---------------------------------------------
    % Then the xzsubset
    %
    % Create a ray that starts in xpoint and extends parallel to the x-axis
    % in the positive x-direction, that is:
    %       z = xpoints(:,3):
    %       xstart = xpoints(:,1);

    if ~isempty(xzsubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(xzsubset)),' xz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(bigplanelist(planelist(possibleones(xzsubset))));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(possibleones(xzsubset),3);
        xstart = xpoints(possibleones(xzsubset),1);
        edgecrossings = zeros(size(xzsubset));

        % Use a parametric representation for each edge:
        % x_edge = x_1 + t*(x_2 - x_1)
        % z_edge = z_1 + t*(z_2 - z_1)
        % Find t by setting z_ray = z_edge

        for ii = 1:max(edgenumbers)
            z1 = corners(planecorners(bigplanelist(planelist(possibleones(xzsubset))),ii),3);
            z2 = corners(planecorners(bigplanelist(planelist(possibleones(xzsubset))),ii+1),3);
            tedgecrossing = (zray - z1)./(z2-z1);

            x1 = corners(planecorners(bigplanelist(planelist(possibleones(xzsubset))),ii),1);
            x2 = corners(planecorners(bigplanelist(planelist(possibleones(xzsubset))),ii+1),1);

            xedge = x1 + tedgecrossing.*(x2-x1);

            edgecrossings = edgecrossings + (tedgecrossing > 0 & tedgecrossing < 1 & xedge > xstart).*(ii <= numberofedgestocheck);
%             closeones = find(abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc );
             closeones = find( (abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc) & xedge > xstart );
            if ~isempty(closeones)
%                disp(['               WARNING:  ',int2str(length(closeones)),' very close hits']) 
               edgecrossings(closeones) = edgecrossings(closeones) + ( ( (z1(closeones)>=zray) + (z2(closeones)>=zray) )== 1);               
             end

        end

        hitvec(possibleones(xzsubset)) = (edgecrossings>0);

        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings>0))),' survived the xzplane projections test:'])  
        end
    end
        
        
    %---------------------------------------------
    % Third the yzsubset
    %
    % Create a ray that starts in xpoint and extends parallel to the y-axis
    % in the positive y-direction, that is:
    %       z = xpoints(:,3):
    %       ystart = xpoints(:,2);

    if ~isempty(yzsubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(yzsubset)),' yz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(bigplanelist(planelist(possibleones(yzsubset))));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(possibleones(yzsubset),3);
        ystart = xpoints(possibleones(yzsubset),2);
        edgecrossings = zeros(size(yzsubset));

        % Use a parametric representation for each edge:
        % y_edge = y_1 + t*(y_2 - y_1)
        % z_edge = z_1 + t*(z_2 - z_1)
        % Find t by setting z_ray = z_edge

        for ii = 1:max(edgenumbers)
            z1 = corners(planecorners(bigplanelist(planelist(possibleones(yzsubset))),ii),3);
            z2 = corners(planecorners(bigplanelist(planelist(possibleones(yzsubset))),ii+1),3);
            tedgecrossing = (zray - z1)./(z2-z1);

            y1 = corners(planecorners(bigplanelist(planelist(possibleones(yzsubset))),ii),2);
            y2 = corners(planecorners(bigplanelist(planelist(possibleones(yzsubset))),ii+1),2);

            yedge = y1 + tedgecrossing.*(y2-y1);

            edgecrossings = edgecrossings + (tedgecrossing > 0 & tedgecrossing < 1 & yedge > ystart).*(ii <= numberofedgestocheck);
%             closeones = find(abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc );
            closeones = find( (abs(tedgecrossing)<geomacc | abs(tedgecrossing-1)<geomacc) & yedge > ystart );
            if ~isempty(closeones)
%                disp(['               WARNING:  ',int2str(length(closeones)),' very close hits']) 
                edgecrossings(closeones) = edgecrossings(closeones) + ( ( (z1(closeones)>=zray) + (z2(closeones)>=zray) )== 1);                               
            end

        end

        hitvec(possibleones(yzsubset)) = (edgecrossings>0);

        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings>0))),' survived the yzplane projections test:'])      
        end
    end
        
end
