function [hitvec,edgehit,edgehitnumbers,cornerhit,cornerhitnumbers] = EDpoinpla(xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs,geomacc,showtext)
% EDpoinpla - Detects if one or more points are inside a number of finite planes. 
% If one point is given as input, it will be checked against all
% planes in a list of planes. If N points are given as inputs, a list of
% planes should have the N planes, and each point will be checked against
% its corresponding plane.
%
% Uses 2D-projected ray casting, with the ray always cast in the positive
% x-direction (for planes that are projected onto the xy-plane etc).
%
% Input parameters:
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
%   geomacc         (optional) The value in meters, by which it is
%                   determined if a point belongs to a finite plane or not.
%   showtext        (optinal)  0 -> no text displayed on screen. Default: 0
%
% Output parameters:
%   hitvec			List, [N,1], with 1 or 0 indicating whether a point is
%                   inside or outside the plane given by planelist.
%   edgehit         List, [N,1], with 1 or 0 indicating if a hit was right at the
%                   edge of a plane. These hits were not marked in hitvec.
%   edgehitnumbers   List, [N,1], with edge numbers for the edgehits
%                    indicated by the "edgehit" list. The edge numbers are
%                    simply 1,2,3,... for each plane, gien by the order of 
%                    the corners for that plane.
%   cornerhit        List, [N,1], with 1 or 0 indicating if a hit was right at the
%                    corner of a plane. These hits were not marked in hitvec.
%   cornerhitnumbers List, [N,1], with corners numbers for the cornerhits
%                    indicated by the "cornerhit" list. The corner numbers are
%                    simply 1,2,3,... for each plane, gien by the order of 
%                    the corners for that plane.
%
% Uses no special subroutines
%
% Peter Svensson (peter.svensson@ntnu.no) 22 June 2021
% 
% [hitvec,edgehit,edgehitnumbers,cornerhit] = EDpoinpla(xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs,geomacc,showtext);

% 100204 Functioning version
% 111201 Fixed a bug: when the ray (from the hitpoint in the positive x-direction) passed exactly
%		 through a vertex, then there was an error. For those rare cases an extra ray is shot in a random direction.
% 27 Nov. 2017 Copied from ESIE2toolbox
% 15 Jan. 2018 Added the detection of edge hits and corner hits
% 15 Mar 2018 Fixed a bug which happened when different planes had
%             different numbers of corners.
% 18 June 2018 Introduced geomacc as an optional input parameter with
% default value 1e-9.
% 14 March 2021 Fixed serious errors that happened with edges parallel to
% xy-axes. Also, edge hits and corner hits were not registered correctly.
% 22 June 2021 Fixed a bug around line 83, which strangely has not lead to
% an error before.

if nargin < 10
    showtext = 0;
end
if nargin < 9
    geomacc = 1e-12;
end

npoints = size(xpoints,1);
nplanestotest = length(planelist);
if npoints == 1
	xpoints = xpoints(ones(nplanestotest,1),:);
end

% Fixed bug 15 Mar 2018. Before, the matrix planecorners could have two
% columns of zeros at the end, but the repetition of the first corner
% number happened after those zeros. So, now we remove those zero columns
% first.
%
% Fixed bug 22 June 2021. nplanestotest does not tell how many rows
% planecorners has, so the line below was changed.

% if nplanestotest == 1
if size(planecorners,1) == 1
   dataincolumns = sign(planecorners); 
else
    dataincolumns = any(planecorners);
end
ncolumnswithdata = sum(double(dataincolumns));
if ncolumnswithdata < size(planecorners,2)
   planecorners = planecorners(:,1:ncolumnswithdata);
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

iv = uint32(find(xpoints(:,1) > maxvals(planelist,1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,1) < minvals(planelist(possibleones),1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) > maxvals(planelist(possibleones),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) < minvals(planelist(possibleones),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) > maxvals(planelist(possibleones),3)));
clear maxvals
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) < minvals(planelist(possibleones),3)));
possibleones(iv) = [];
clear minvals iv

hitvec = uint8(zeros(nplanestotest,1));
hitvec(possibleones) = 1;
edgehit          = zeros(nplanestotest,1,'uint32');
edgehitnumbers   = zeros(nplanestotest,1,'uint32');
cornerhit        = zeros(nplanestotest,1,'uint32');
cornerhitnumbers = zeros(nplanestotest,1,'uint32');

nposs = length(possibleones);

if showtext >= 4
    disp(['         Of the ',int2str(npoints),' points,'])
    disp(['         ',int2str(length(possibleones)),' survived the cube test:'])   
    disp(num2str(possibleones.'))
end

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

    A = abs((planenvecs(planelist(possibleones),:)));
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
    
    yzsubset = possibleones(find(finalcolno==1));
    xzsubset = possibleones(find(finalcolno==2));
    xysubset = possibleones(find(finalcolno==3));
    
    %---------------------------------------------
    % First the xysubset
    %
    % Create a ray that starts in xpoint and extends parallel to the x-axis
    % in the positive x-direction, that is:
    %       y = xpoints(:,2):
    %       xstart = xpoints(:,1);
    
    if ~isempty(xysubset)
        if showtext >= 4
            disp(['            ',int2str(length(xysubset)),' xy-projected points:'])
        end
        
        numberofedgestocheck = ncornersperplanevec(planelist(xysubset));    
        edgenumbers = unique(numberofedgestocheck);

        yray = xpoints(xysubset,2);
        xstart = xpoints(xysubset,1);
        edgecrossings = zeros(size(xysubset));
        
        closetoedge = zeros(size(xysubset));
        addto_closetoedge = zeros(size(xysubset));
        edgenumbers_that_were_hit = zeros(size(xysubset));
        
        closetocorner = zeros(size(xysubset));
        addto_closetocorner = zeros(size(xysubset)); 
        edgeswithcorners_that_were_hit = zeros(length(xysubset),2);
        insidehorizontaledge = zeros(size(xysubset));
        smallvertdistance = zeros(size(xysubset));
        closetocornerofhorizontaledge = zeros(size(xysubset));
        
        % Use a parametric representation for each edge:
        % x_edge = x_1 + t*(x_2 - x_1)
        % y_edge = y_1 + t*(y_2 - y_1)
        % Find t by setting y_ray = y_edge

        for ii = 1:double(max(edgenumbers))

            y1 = corners(planecorners(planelist(xysubset),ii),2);
            y2 = corners(planecorners(planelist(xysubset),ii+1),2);
            
            horizontaledges = (y1==y2);
            nonhorizontaledges = 1 - horizontaledges;
            ivhor = find(horizontaledges);
            
            tedgecrossing = (yray - y1)./(y2-y1);

            talmostzero = abs(tedgecrossing)<geomacc;
            talmostone  = abs(tedgecrossing-1)<geomacc;
            tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
            tendpoint_countablehit = (talmostzero==1 & y2 < 0) + (talmostone==1 & y1 < 0);

            x1 = corners(planecorners(planelist(xysubset),ii),1);
            x2 = corners(planecorners(planelist(xysubset),ii+1),1);
            xedge = x1 + tedgecrossing.*(x2-x1);
            xedgeveryclosetostart = abs(xedge-xstart)<geomacc;

            edgecrossings = edgecrossings + ...
                (tinside & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                (tendpoint_countablehit & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;

            addto_closetocorner = (xedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck);
            edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));

            addto_closetoedge = (xedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck);
            closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
            if any(addto_closetoedge)
               edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii; 
            end
            if ~isempty(ivhor)
                insidehorizontaledge = insidehorizontaledge*0;
                smallvertdistance = smallvertdistance*0;
                smallvertdistance(ivhor) = abs(yray(ivhor)-y1(ivhor))<geomacc;               
                insidehorizontaledge(ivhor) = ...
                    ( ( xstart(ivhor)-x1(ivhor)>geomacc & x2(ivhor)-xstart(ivhor)>geomacc ) | ...
                      ( xstart(ivhor)-x2(ivhor)>geomacc & x1(ivhor)-xstart(ivhor)>geomacc ) ).* ...
                      smallvertdistance(ivhor);
                addto_closetoedge = (horizontaledges == 1 & insidehorizontaledge == 1 );
                closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                if any(addto_closetoedge)
                    edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                end
                closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                closetocornerofhorizontaledge(ivhor) = ...
                    (abs( xstart(ivhor)-x1(ivhor) ) < geomacc | ...
                     abs( xstart(ivhor)-x2(ivhor) ) < geomacc ).* ...
                      smallvertdistance(ivhor);
                addto_closetocorner = (horizontaledges == 1 & closetocornerofhorizontaledge == 1 );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));

            end 

        end
        
        hitvec(xysubset) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
        
        cornerhit(xysubset) = sign(closetocorner);
        iv_replace_highest_edgenumber_with_zero = ...
            find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                 edgeswithcorners_that_were_hit(:,1)== 1);
        edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
        cornerhitnumbers(xysubset) = max(edgeswithcorners_that_were_hit.').';

        edgehit(xysubset) = (closetoedge>0 );
        edgehitnumbers(xysubset) = edgenumbers_that_were_hit;
        
        if showtext >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the xyplane projections test:'])  
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
        if showtext >= 4
            disp(['            ',int2str(length(xzsubset)),' xz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(planelist(xzsubset));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(xzsubset,3);
        xstart = xpoints(xzsubset,1);
        edgecrossings = zeros(size(xzsubset));
        closetoedge = zeros(size(xzsubset));
        addto_closetoedge = zeros(size(xzsubset));
        edgenumbers_that_were_hit = zeros(size(xzsubset));
        
        closetocorner = zeros(size(xzsubset));
        addto_closetocorner = zeros(size(xzsubset)); 
        edgeswithcorners_that_were_hit = zeros(length(xzsubset),2);
        insidehorizontaledge = zeros(size(xzsubset));
        smallvertdistance = zeros(size(xzsubset));
        closetocornerofhorizontaledge = zeros(size(xzsubset));
        
        % Use a parametric representation for each edge:
        % x_edge = x_1 + t*(x_2 - x_1)
        % z_edge = z_1 + t*(z_2 - z_1)
        % Find t by setting z_ray = z_edge
      
       for ii = 1:double(max(edgenumbers))

            z1 = corners(planecorners(planelist(xzsubset),ii),3);
            z2 = corners(planecorners(planelist(xzsubset),ii+1),3);
            
            horizontaledges = (z1==z2);
            nonhorizontaledges = 1 - horizontaledges;
            ivhor = find(horizontaledges);
            
            tedgecrossing = (zray - z1)./(z2-z1);

            talmostzero = abs(tedgecrossing)<geomacc;
            talmostone  = abs(tedgecrossing-1)<geomacc;
            tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
            tendpoint_countablehit = (talmostzero==1 & z2 < 0) + (talmostone==1 & z1 < 0);
    
            x1 = corners(planecorners(planelist(xzsubset),ii),1);
            x2 = corners(planecorners(planelist(xzsubset),ii+1),1);
            xedge = x1 + tedgecrossing.*(x2-x1);
            xedgeveryclosetostart = abs(xedge-xstart)<geomacc;

            edgecrossings = edgecrossings + ...
                (tinside & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                (tendpoint_countablehit & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;

            addto_closetocorner = (xedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck);
            edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));

            addto_closetoedge = (xedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck);
            closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
            if any(addto_closetoedge)
               edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii; 
            end
            if ~isempty(ivhor)
                insidehorizontaledge = insidehorizontaledge*0;
                smallvertdistance = smallvertdistance*0;
                smallvertdistance(ivhor) = abs(zray(ivhor)-z1(ivhor))<geomacc;               
                insidehorizontaledge(ivhor) = ...
                    ( ( xstart(ivhor)-x1(ivhor)>geomacc & x2(ivhor)-xstart(ivhor)>geomacc ) | ...
                      ( xstart(ivhor)-x2(ivhor)>geomacc & x1(ivhor)-xstart(ivhor)>geomacc ) ).* ...
                      smallvertdistance(ivhor);
                addto_closetoedge = (horizontaledges == 1 & insidehorizontaledge == 1 );
                closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                if any(addto_closetoedge)
                    edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                end
                closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                closetocornerofhorizontaledge(ivhor) = ...
                    (abs( xstart(ivhor)-x1(ivhor) ) < geomacc | ...
                     abs( xstart(ivhor)-x2(ivhor) ) < geomacc ).* ...
                      smallvertdistance(ivhor);
                addto_closetocorner = (horizontaledges == 1 & closetocornerofhorizontaledge == 1 );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
            end                        
       end
       
       hitvec(xzsubset) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
        
        cornerhit(xzsubset) = sign(closetocorner);
        iv_replace_highest_edgenumber_with_zero = ...
            find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                 edgeswithcorners_that_were_hit(:,1)== 1);
        edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
        cornerhitnumbers(xzsubset) = max(edgeswithcorners_that_were_hit.').';

        edgehit(xzsubset) = (closetoedge>0 );
        edgehitnumbers(xzsubset) = edgenumbers_that_were_hit;

        if showtext >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the xzplane projections test:'])  
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
        if showtext >= 4
            disp(['            ',int2str(length(yzsubset)),' yz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(planelist(yzsubset));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(yzsubset,3);
        ystart = xpoints(yzsubset,2);

        edgecrossings = zeros(size(yzsubset));
        closetoedge = zeros(size(yzsubset));
        addto_closetoedge = zeros(size(yzsubset));
        edgenumbers_that_were_hit = zeros(size(yzsubset));
        
        closetocorner = zeros(size(yzsubset));
        addto_closetocorner = zeros(size(yzsubset)); 
        edgeswithcorners_that_were_hit = zeros(length(yzsubset),2);
        insidehorizontaledge = zeros(size(yzsubset));
        smallvertdistance = zeros(size(yzsubset));
        closetocornerofhorizontaledge = zeros(size(yzsubset));
        
        % Use a parametric representation for each edge:
        % y_edge = y_1 + t*(y_2 - y_1)
        % z_edge = z_1 + t*(z_2 - z_1)
        % Find t by setting z_ray = z_edge        
      
       for ii = 1:double(max(edgenumbers))
            z1 = corners(planecorners(planelist(yzsubset),ii),3);
            z2 = corners(planecorners(planelist(yzsubset),ii+1),3);
            
            horizontaledges = (z1==z2);
            nonhorizontaledges = 1 - horizontaledges;
            ivhor = find(horizontaledges);
            
            tedgecrossing = (zray - z1)./(z2-z1);

            talmostzero = abs(tedgecrossing)<geomacc;
            talmostone  = abs(tedgecrossing-1)<geomacc;
            tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
            tendpoint_countablehit = (talmostzero==1 & z2 < 0) + (talmostone==1 & z1 < 0);
    
            y1 = corners(planecorners(planelist(yzsubset),ii),2);
            y2 = corners(planecorners(planelist(yzsubset),ii+1),2);
            yedge = y1 + tedgecrossing.*(y2-y1);
            yedgeveryclosetostart = abs(yedge-ystart)<geomacc;

            edgecrossings = edgecrossings + ...
                (tinside & yedge > ystart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                (tendpoint_countablehit & yedge > ystart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;

            addto_closetocorner = (yedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck);
            edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
            closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));

            addto_closetoedge = (yedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck);
            closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
            if any(addto_closetoedge)
               edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii; 
            end
            if ~isempty(ivhor)
                insidehorizontaledge = insidehorizontaledge*0;
                smallvertdistance = smallvertdistance*0;
                smallvertdistance(ivhor) = abs(zray(ivhor)-z1(ivhor))<geomacc;               
                insidehorizontaledge(ivhor) = ...
                    ( ( ystart(ivhor)-y1(ivhor)>geomacc & y2(ivhor)-ystart(ivhor)>geomacc ) | ...
                      ( ystart(ivhor)-y2(ivhor)>geomacc & y1(ivhor)-ystart(ivhor)>geomacc ) ).* ...
                      smallvertdistance(ivhor);
                addto_closetoedge = (horizontaledges == 1 & insidehorizontaledge == 1 );
                closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                if any(addto_closetoedge)
                    edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                end
                closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                closetocornerofhorizontaledge(ivhor) = ...
                    (abs( ystart(ivhor)-y1(ivhor) ) < geomacc | ...
                     abs( ystart(ivhor)-y2(ivhor) ) < geomacc ).* ...
                      smallvertdistance(ivhor);
                addto_closetocorner = (horizontaledges == 1 & closetocornerofhorizontaledge == 1 );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
            end                        
       end
       
       hitvec(yzsubset) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
        
        cornerhit(yzsubset) = sign(closetocorner);
        iv_replace_highest_edgenumber_with_zero = ...
            find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                 edgeswithcorners_that_were_hit(:,1)== 1);
        edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
        cornerhitnumbers(yzsubset) = max(edgeswithcorners_that_were_hit.').';

        edgehit(yzsubset) = (closetoedge>0 );
        edgehitnumbers(yzsubset) = edgenumbers_that_were_hit;
        
        if showtext >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the yzplane projections test:'])      
        end
    end
    
end
