function [hitvec,edgehit,edgehitnumbers,cornerhit,cornerhitnumbers] = ...
    EDpoinpla(xpoints,planelist,minvals,maxvals,planecorners,corners,...
    ncornersperplanevec,planenvecs,geomacc,showtext)
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
%   planelist       List, [nplanes,1], of plane numbers to check the N points
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
%                    simply 1,2,3,... for each plane, given by the order of 
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
% Peter Svensson (peter.svensson@ntnu.no) 31 Oct. 2023
% 
% [hitvec,edgehit,edgehitnumbers,cornerhit] = EDpoinpla(xpoints,planelist,...
% minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs,geomacc,showtext);

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
% 16 Aug 2021 Converted two logical arrays (addto_closetorcorner and 
% addto_closetoedge) to double since Matlab 2018 protested, but Matlab 2020 did not.
% 8 June 2022: Fixed a bug that made the function crash when planes with different numbers
% of corners ended up in the same xysubset etc.
% 31 Oct. 2023 Added the possibility that this function is called with some
% 3-corner planes and some 4-corner planes etc.
% 31 Oct. 2023 Found a bug which would happen rarely: if the ray goes
% exactly through a corner, then the hit was micounted.

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

% For some uses of this function, the rows in planecorners might have
% different numbers of corners. That is fixed here (31 Oct. 2023).

columnswithzeros = sum(planecorners==0);

if any(columnswithzeros)
    iv = find(columnswithzeros);
    for ii = 1:length(iv)
        colwzero = iv(ii);
        rowswzero = find(planecorners(:,colwzero)==0);
        planecorners(rowswzero,colwzero) = planecorners(rowswzero,1);
    end
end

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
    % Fixed a bug 8 June 2022
    % If one subset contained planes with different numbers of corners,
    % this function crashed, because the for-loop stepping through the
    % number of edges used the highest number of edges in the list and
    % therefore hit zeros. In each subset, we have to first check if all
    % planes have the same number of edges. If they do, the old code can be
    % used. If they don't we create subsubsets, using cells (since each
    % subsubset will have different numbers of entries and we don't know
    % how many subsubsets will be needed.

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

        % 8 June 2022 We create subsubsets as cells, one for each value in
        % edgenumbers.

        xysubsubsets = cell(length(edgenumbers),1);
        for ii = 1:length(edgenumbers)
            xysubsubsets{ii} = xysubset( find(numberofedgestocheck==edgenumbers(ii))  );
        end

        for jj = 1:length(edgenumbers)
            numberofedgestocheck = ncornersperplanevec(planelist(xysubsubsets{jj}));  

            yray = xpoints(xysubsubsets{jj},2);
            xstart = xpoints(xysubsubsets{jj},1);
            edgecrossings = zeros(size(xysubsubsets{jj}));
            
            closetoedge = zeros(size(xysubsubsets{jj}));
            addto_closetoedge = zeros(size(xysubsubsets{jj}));
            edgenumbers_that_were_hit = zeros(size(xysubsubsets{jj}));
            
            closetocorner = zeros(size(xysubsubsets{jj}));
            addto_closetocorner = zeros(size(xysubsubsets{jj})); 
            edgeswithcorners_that_were_hit = zeros(length(xysubsubsets{jj}),2);
            insidehorizontaledge = zeros(size(xysubsubsets{jj}));
            smallvertdistance = zeros(size(xysubsubsets{jj}));
            closetocornerofhorizontaledge = zeros(size(xysubsubsets{jj}));
            
            % Use a parametric representation for each edge:
            % x_edge = x_1 + t*(x_2 - x_1)
            % y_edge = y_1 + t*(y_2 - y_1)
            % Find t by setting y_ray = y_edge

            for ii = 1:double((edgenumbers(jj)))
    
                y1 = corners(planecorners(planelist(xysubsubsets{jj}),ii),2);
                y2 = corners(planecorners(planelist(xysubsubsets{jj}),ii+1),2);
                
                horizontaledges = (y1==y2);
                nonhorizontaledges = 1 - horizontaledges;
                ivhor = find(horizontaledges);
                
                tedgecrossing = (yray - y1)./(y2-y1);
    
                talmostzero = abs(tedgecrossing)<geomacc;
                talmostone  = abs(tedgecrossing-1)<geomacc;
                tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
                % Error found 31 Oct. 2023: a corner hit should be
                % countable if the other edge endpoint is below **the ray**, 
                % not below zero!
                % tendpoint_countablehit = (talmostzero==1 & y2 < 0) + (talmostone==1 & y1 < 0);
                tendpoint_countablehit = (talmostzero==1 & y2 < yray) + (talmostone==1 & y1 < yray);
    
                x1 = corners(planecorners(planelist(xysubsubsets{jj}),ii),1);
                x2 = corners(planecorners(planelist(xysubsubsets{jj}),ii+1),1);
                xedge = x1 + tedgecrossing.*(x2-x1);
                xedgeveryclosetostart = abs(xedge-xstart)<geomacc;

                edgecrossings = edgecrossings + ...
                    (tinside & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                    (tendpoint_countablehit & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;
                % 16 Aug 2021 Converted the logical addto_closetocorner to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetocorner = double( (xedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck) );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
    
                % 16 Aug 2021 Converted the logical addto_closetoedge to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetoedge = double( (xedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck) );
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
                    % 16 Aug 2021 Converted the logical addto_closetoedge to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetoedge = double( (horizontaledges == 1 & insidehorizontaledge == 1 ) );
                    closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                    if any(addto_closetoedge)
                        edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                    end
                    closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                    closetocornerofhorizontaledge(ivhor) = ...
                        (abs( xstart(ivhor)-x1(ivhor) ) < geomacc | ...
                         abs( xstart(ivhor)-x2(ivhor) ) < geomacc ).* ...
                          smallvertdistance(ivhor);
                    % 16 Aug 2021 Converted the logical addto_closetocorner to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetocorner = double( (horizontaledges == 1 & closetocornerofhorizontaledge == 1 ) );
                    edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
    
                end 
    
            end       

            hitvec(xysubsubsets{jj}) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
            
            cornerhit(xysubsubsets{jj}) = sign(closetocorner);
            iv_replace_highest_edgenumber_with_zero = ...
                find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                     edgeswithcorners_that_were_hit(:,1)== 1);
            edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
            cornerhitnumbers(xysubsubsets{jj}) = max(edgeswithcorners_that_were_hit.').';
    
            edgehit(xysubsubsets{jj}) = (closetoedge>0 );
            edgehitnumbers(xysubsubsets{jj}) = edgenumbers_that_were_hit;
        end
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

        % 8 June 2022 We create subsubsets as cells, one for each value in
        % edgenumbers.

        xzsubsubsets = cell(length(edgenumbers),1);
        for ii = 1:length(edgenumbers)
            xzsubsubsets{ii} = xzsubset( find(numberofedgestocheck==edgenumbers(ii))  );
        end

       for jj = 1:length(edgenumbers)

            numberofedgestocheck = ncornersperplanevec(planelist(xzsubsubsets{jj}));    
        
            zray = xpoints(xzsubsubsets{jj},3);
            xstart = xpoints(xzsubsubsets{jj},1);
            edgecrossings = zeros(size(xzsubsubsets{jj}));
            closetoedge = zeros(size(xzsubsubsets{jj}));
            addto_closetoedge = zeros(size(xzsubsubsets{jj}));
            edgenumbers_that_were_hit = zeros(size(xzsubsubsets{jj}));
            
            closetocorner = zeros(size(xzsubsubsets{jj}));
            addto_closetocorner = zeros(size(xzsubsubsets{jj})); 
            edgeswithcorners_that_were_hit = zeros(length(xzsubsubsets{jj}),2);
            insidehorizontaledge = zeros(size(xzsubsubsets{jj}));
            smallvertdistance = zeros(size(xzsubsubsets{jj}));
            closetocornerofhorizontaledge = zeros(size(xzsubsubsets{jj}));
            
            % Use a parametric representation for each edge:
            % x_edge = x_1 + t*(x_2 - x_1)
            % z_edge = z_1 + t*(z_2 - z_1)
            % Find t by setting z_ray = z_edge
          
            for ii = 1:double((edgenumbers(jj)))
    
                z1 = corners(planecorners(planelist(xzsubsubsets{jj}),ii),3);
                z2 = corners(planecorners(planelist(xzsubsubsets{jj}),ii+1),3);
                
                horizontaledges = (z1==z2);
                nonhorizontaledges = 1 - horizontaledges;
                ivhor = find(horizontaledges);
                
                tedgecrossing = (zray - z1)./(z2-z1);
    
                talmostzero = abs(tedgecrossing)<geomacc;
                talmostone  = abs(tedgecrossing-1)<geomacc;
                tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
                % Error found 31 Oct. 2023: a corner hit should be
                % countable if the other edge endpoint is below **the ray**, 
                % not below zero!
                % tendpoint_countablehit = (talmostzero==1 & z2 < 0) + (talmostone==1 & z1 < 0);
                tendpoint_countablehit = (talmostzero==1 & z2 < zray) + (talmostone==1 & z1 < zray);
        
                x1 = corners(planecorners(planelist(xzsubsubsets{jj}),ii),1);
                x2 = corners(planecorners(planelist(xzsubsubsets{jj}),ii+1),1);
                xedge = x1 + tedgecrossing.*(x2-x1);
                xedgeveryclosetostart = abs(xedge-xstart)<geomacc;
    
                edgecrossings = edgecrossings + ...
                    (tinside & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                    (tendpoint_countablehit & xedge > xstart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;
    
                % 16 Aug 2021 Converted the logical addto_closetocorner to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetocorner = double( (xedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck) );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
    
                % 16 Aug 2021 Converted the logical addto_closetoedge to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetoedge = double( (xedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck) );
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
    
                    % 16 Aug 2021 Converted the logical addto_closetoedge to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetoedge = double( (horizontaledges == 1 & insidehorizontaledge == 1 ) );
                    closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                    if any(addto_closetoedge)
                        edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                    end
                    closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                    closetocornerofhorizontaledge(ivhor) = ...
                        (abs( xstart(ivhor)-x1(ivhor) ) < geomacc | ...
                         abs( xstart(ivhor)-x2(ivhor) ) < geomacc ).* ...
                          smallvertdistance(ivhor);
                      
                    % 16 Aug 2021 Converted the logical addto_closetocorner to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetocorner = double( (horizontaledges == 1 & closetocornerofhorizontaledge == 1 ) );
                    edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
                end                        
           end
       
           hitvec(xzsubsubsets{jj}) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
            
            cornerhit(xzsubsubsets{jj}) = sign(closetocorner);
            iv_replace_highest_edgenumber_with_zero = ...
                find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                     edgeswithcorners_that_were_hit(:,1)== 1);
            edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
            cornerhitnumbers(xzsubsubsets{jj}) = max(edgeswithcorners_that_were_hit.').';
    
            edgehit(xzsubsubsets{jj}) = (closetoedge>0 );
            edgehitnumbers(xzsubsubsets{jj}) = edgenumbers_that_were_hit;
       end

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

        % 8 June 2022 We create subsubsets as cells, one for each value in
        % edgenumbers.

        yzsubsubsets = cell(length(edgenumbers),1);
        for ii = 1:length(edgenumbers)
            yzsubsubsets{ii} = yzsubset( find(numberofedgestocheck==edgenumbers(ii))  );
        end

       for jj = 1:length(edgenumbers)
            zray = xpoints(yzsubsubsets{jj},3);
            ystart = xpoints(yzsubsubsets{jj},2);
    
            edgecrossings = zeros(size(yzsubsubsets{jj}));
            closetoedge = zeros(size(yzsubsubsets{jj}));
            addto_closetoedge = zeros(size(yzsubsubsets{jj}));
            edgenumbers_that_were_hit = zeros(size(yzsubsubsets{jj}));
            
            closetocorner = zeros(size(yzsubsubsets{jj}));
            addto_closetocorner = zeros(size(yzsubsubsets{jj})); 
            edgeswithcorners_that_were_hit = zeros(length(yzsubsubsets{jj}),2);
            insidehorizontaledge = zeros(size(yzsubsubsets{jj}));
            smallvertdistance = zeros(size(yzsubsubsets{jj}));
            closetocornerofhorizontaledge = zeros(size(yzsubsubsets{jj}));
            
            % Use a parametric representation for each edge:
            % y_edge = y_1 + t*(y_2 - y_1)
            % z_edge = z_1 + t*(z_2 - z_1)
            % Find t by setting z_ray = z_edge        
          
            for ii = 1:double((edgenumbers(jj)))
                z1 = corners(planecorners(planelist(yzsubsubsets{jj}),ii),3);
                z2 = corners(planecorners(planelist(yzsubsubsets{jj}),ii+1),3);
                
                horizontaledges = (z1==z2);
                nonhorizontaledges = 1 - horizontaledges;
                ivhor = find(horizontaledges);
                
                tedgecrossing = (zray - z1)./(z2-z1);
    
                talmostzero = abs(tedgecrossing)<geomacc;
                talmostone  = abs(tedgecrossing-1)<geomacc;
                tinside = tedgecrossing > 0 & tedgecrossing < 1 & talmostzero == 0 & talmostone == 0;
                % Error found 31 Oct. 2023: a corner hit should be
                % countable if the other edge endpoint is below **the ray**, 
                % not below zero!
                % tendpoint_countablehit = (talmostzero==1 & z2 < 0) + (talmostone==1 & z1 < 0);
                tendpoint_countablehit = (talmostzero==1 & z2 < zray) + (talmostone==1 & z1 < zray);
        
                y1 = corners(planecorners(planelist(yzsubsubsets{jj}),ii),2);
                y2 = corners(planecorners(planelist(yzsubsubsets{jj}),ii+1),2);
                yedge = y1 + tedgecrossing.*(y2-y1);
                yedgeveryclosetostart = abs(yedge-ystart)<geomacc;
    
                edgecrossings = edgecrossings + ...
                    (tinside & yedge > ystart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) + ...
                    (tendpoint_countablehit & yedge > ystart).*(ii <= numberofedgestocheck).*(nonhorizontaledges) ;
    
                % 16 Aug 2021 Converted the logical addto_closetocorner to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetocorner = double( (yedgeveryclosetostart==1 & (talmostzero+talmostone) > 0).*(ii <= numberofedgestocheck) );
                edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
    
                % 16 Aug 2021 Converted the logical addto_closetoedge to
                % double since some Matlab versions protested against the
                % addition/multiplication on the lines below.
                addto_closetoedge = double( (yedgeveryclosetostart==1 & tinside & (talmostzero+talmostone) == 0).*(ii <= numberofedgestocheck) );
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
                      
                    % 16 Aug 2021 Converted the logical addto_closetoedge to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetoedge = double( (horizontaledges == 1 & insidehorizontaledge == 1 ) );
                    closetoedge = closetoedge + addto_closetoedge.*(1-sign(closetoedge));
                    if any(addto_closetoedge)
                        edgenumbers_that_were_hit = edgenumbers_that_were_hit + sign(addto_closetoedge)*ii;                     
                    end
                    closetocornerofhorizontaledge = closetocornerofhorizontaledge*0;
                    closetocornerofhorizontaledge(ivhor) = ...
                        (abs( ystart(ivhor)-y1(ivhor) ) < geomacc | ...
                         abs( ystart(ivhor)-y2(ivhor) ) < geomacc ).* ...
                          smallvertdistance(ivhor);
                    % 16 Aug 2021 Converted the logical addto_closetocorner to
                    % double since some Matlab versions protested against the
                    % addition/multiplication on the lines below.
                    addto_closetocorner = double( (horizontaledges == 1 & closetocornerofhorizontaledge == 1 ) );
                    edgeswithcorners_that_were_hit(:,2) = edgeswithcorners_that_were_hit(:,2) + sign(addto_closetocorner).*(sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    edgeswithcorners_that_were_hit(:,1) = edgeswithcorners_that_were_hit(:,1) + sign(addto_closetocorner).*(1-sign(edgeswithcorners_that_were_hit(:,1)))*ii;
                    closetocorner = closetocorner + addto_closetocorner.*(1-sign(closetocorner));
                end                        
           end
           
           hitvec(yzsubsubsets{jj}) = (edgecrossings==1 & closetocorner==0 & closetoedge==0);
            
            cornerhit(yzsubsubsets{jj}) = sign(closetocorner);
            iv_replace_highest_edgenumber_with_zero = ...
                find(edgeswithcorners_that_were_hit(:,2)== max(edgenumbers)  & ...
                     edgeswithcorners_that_were_hit(:,1)== 1);
            edgeswithcorners_that_were_hit(iv_replace_highest_edgenumber_with_zero,2) = 0;
            cornerhitnumbers(yzsubsubsets{jj}) = max(edgeswithcorners_that_were_hit.').';
    
            edgehit(yzsubsubsets{jj}) = (closetoedge>0 );
            edgehitnumbers(yzsubsubsets{jj}) = edgenumbers_that_were_hit;

       end

        if showtext >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the yzplane projections test:'])      
        end
    end
    
end
