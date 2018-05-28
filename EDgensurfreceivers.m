function [surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
    surfacegaussorder,EDversionnumber)
% EDgensurfreceivers distributes receivers along all the faces/polygons of
% a polyhedral object. A Gauss-Legendre distribution is used. The function
% can handle general quadrilaterals and triangles.
% The quadrature order is given as an input parameter, and this value will
% be used for the largest polygon, and then scaled down for smaller
% polygons. For triangles, the same quadrature order is used for all.
% 
% Input parameters:
%   planedata           Struct
%   surfacegaussorder   The quadratureorder that will be used for the largest
%                       plane/face/polygon of the model. A value of 5
%                       implies that 5*5 receivers will be generated across
%                       the largest surface.
%   EDversionnumber
%   showtext            (optional)
% 
% Output parameters:
%   surfacerecs         Matrix, [nrec,3], of the coordinates
%   surfacerecnvecs     Matrix, [nrec,3], with the plane normal vector for
%                       each coordinate
%   surfacerecweights   Vector, [nrec,1], with the weight factor for each
%                       receiver point. These weight factors acts as dS in
%                       the Helmholtz integral, that is, they have the unit
%                       m2.
% References for triangle quadrature:
% 1. D.A. Dunavant, "High degree efficient symmetrical Gaussian quadrature
%    rules for the triangle," Int. J. Num. Meth. Eng. 21, pp. 1129-1148
%    (1985).
% 2. P. M. Juhl, "The boundary element method for sound field calculations,"
%    Report nr. 55, The acoustics laboratory, Technical University of Denmark, 1993.
% 
% Uses the function lgwt from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 10 Apr 2018

% 1 Mar 2018 First version
% 8 Mar 2018 Changed the displacement of the surface receivers to 0.1mm
% 4 Apr 2018 Introduced the triangle element quadrature, using gaussTRINEW
% from OpenBEM
% 5 Apr 2018 Small error if no triangular elements. Implemented quadrature
% from the Duvanant paper instead of the gaussTRINEW function.
% 10 Apr 2018 Added some quadrature values

% distancefromsurf = 1e-4;
distancefromsurf = 1e-3;

nplanes = size(planedata.planecorners,1);
if any(planedata.ncornersperplanevec>4)
    error('ERROR: Unfortunately, only triangular and quadrilateral surfaces ',...
          'can be handled by this version of EDgensurfreceivers');
end

iv3 = find(planedata.ncornersperplanevec == 3);
iv4 = find(planedata.ncornersperplanevec == 4);

planeextensions = planedata.maxvals - planedata.minvals;
planeextensions = max(planeextensions.').';
planerelsize = planeextensions/max(planeextensions);

gaussorderperplane = ceil(planerelsize*surfacegaussorder);
if ~isempty(iv3)
    if surfacegaussorder >8 
       disp('   WARNING: surfacegaussorder is decreased to 8 (for the triangular elements)') 
       surfacegaussorder = 8;
    end
    gaussorderperplane(iv3) = 0;
end
    
surfacerecs = [];
surfacerecnvecs = [];
surfacerecweights = [];

% First all quadrilateral elements 

for ii = max(gaussorderperplane):-1:1
   n2 = ii^2;
   iv = find(gaussorderperplane == ii);
   [x,w] = lgwt(ii,0,1);
   
   x = x(end:-1:1);
   xy = [repmat(x,ii,1) reshape(x(:,ones(1,ii)).',n2,1)];
   ww = prod([repmat(w,ii,1) reshape(w(:,ones(1,ii)).',n2,1)],2);   
   
   for jj = 1:length(iv)
       planenumber = iv(jj);
       c1 = planedata.corners(planedata.planecorners(planenumber,1),:);
       c2 = planedata.corners(planedata.planecorners(planenumber,2),:);
       c3 = planedata.corners(planedata.planecorners(planenumber,3),:);
       c4 = planedata.corners(planedata.planecorners(planenumber,4),:);
       xvec1 = c2 - c1;
       alen = norm(xvec1);
       xvec2 = c3 - c4;
       clen = norm(xvec2);
       yvec1 = c4 - c1;
       dlen = norm(yvec1);
       yvec2 = c3 - c2;
       blen = norm(yvec2);
       crossvec1 = c3 - c1;
       plen = norm(crossvec1);
       crossvec2 = c4 - c2;
       qlen = norm(crossvec2);
%        if norm(crossvec1) ~= norm(crossvec2)
%           error('ERROR: general quadrilaterals have not been implemented yet') 
%        end
% % % % % %        length1 = norm(xvec);
% % % % % %        length2 = norm(yvec);
% % % % % %        planearea = length1*length2;

        planearea = sqrt( 4*plen^2*qlen^2 - (blen^2 + dlen^2 - alen^2 - clen^2)^2 )/4;

       recnvecs_oneplane = planedata.planeeqs(planenumber,1:3);
       recnvecs_oneplane = recnvecs_oneplane(ones(n2,1),:);
       recweights_oneplane = ww*planearea;

       startpointsx = c1(ones(n2,1),:);
       startpointsx = startpointsx + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);

       endpointsx = c4(ones(n2,1),:);
       endpointsx = endpointsx + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
       
       startpointsy = c1(ones(n2,1),:);
       startpointsy = startpointsy + yvec1(ones(n2,1),:).*xy(:,[2 2 2]);

       endpointsy = c2(ones(n2,1),:);
       endpointsy = endpointsy + yvec2(ones(n2,1),:).*xy(:,[2 2 2]);

       recs_oneplane = startpointsx + (endpointsx - startpointsx).*xy(:,[2 2 2]);
%        recs_oneplane = startpointsx + (endpointsx - startpointsx).*xy(:,[1 1 1]) + (endpointsy - startpointsy).*xy(:,[2 2 2]);
       recs_oneplane = recs_oneplane + recnvecs_oneplane*distancefromsurf;       
%        startpoint = startpoint(ones(n2,1),:);
%        recs_oneplane = startpoint + xvec(ones(n2,1),:).*xy(:,[1 1 1]) + yvec(ones(n2,1),:).*xy(:,[2 2 2]);
%        recs_oneplane = recs_oneplane + recnvecs_oneplane*distancefromsurf;
       
       surfacerecs = [surfacerecs;recs_oneplane];
       surfacerecnvecs = [surfacerecnvecs;recnvecs_oneplane];
       surfacerecweights = [surfacerecweights;recweights_oneplane];
   end
end

% Now all the triangular elements

if ~isempty(iv3)
%     gaussvalues = gaussTRINEW(surfacegaussorder);
    switch surfacegaussorder
        case {0,1}
           gaussvalues = [1/3 1/3 1]; 
        case {2}
            gaussvalues = [2/3 1/6 1/3;...
                           1/6 2/3 1/3;
                           1/6 1/6 1/3];
        case {3}
            gaussvalues = [1/3 1/3 -9/32; ...
                0.2 0.6 25/96; ...
                0.2 0.2 25/96; ...
                0.6 0.2 25/96];
        case {4}
            gaussvalues = [0.108103018168070  0.445948490915965 0.223381589678011; ...
                           0.445948490915965  0.108103018168070 0.223381589678011; ...
                           0.445948490915965  0.445948490915965 0.223381589678011; ...
                           0.816847572980459  0.091576213509771 0.109951743655322; ...
                           0.091576213509771  0.816847572980459 0.109951743655322; ...
                           0.091576213509771  0.091576213509771 0.109951743655322];
        case {5}   % 7 gauss points
            gaussvalues = [1/3  1/3  0.225;...
                           0.059715871789770 0.470142064105115 0.132394152788506;...
                           0.470142064105115 0.059715871789770 0.132394152788506;...
                           0.470142064105115 0.470142064105115 0.132394152788506;...
                           0.797426985353087 0.101286507323456 0.125939180544827;...
                           0.101286507323456 0.797426985353087 0.125939180544827;...
                           0.101286507323456 0.101286507323456 0.125939180544827];
        case{6}  % 12 gauss points
            gaussvalues = [0.501426509658179 0.249286745170910 0.116786275726379;...
                           0.249286745170910 0.501426509658179 0.116786275726379;...
                           0.249286745170910 0.249286745170910 0.116786275726379;...
                           0.873821971016996 0.063089014491502 0.050844906370207;...
                           0.063089014491502 0.873821971016996 0.050844906370207;...
                           0.063089014491502 0.063089014491502 0.050844906370207;...                           
                           0.053145049844817 0.310352451033784 0.082851075618374;...
                           0.310352451033784 0.053145049844817 0.082851075618374;...
                           0.053145049844817 0.636502499121399 0.082851075618374;...
                           0.636502499121399 0.053145049844817 0.082851075618374;...
                           0.636502499121399 0.310352451033784 0.082851075618374;...
                           0.310352451033784 0.636502499121399 0.082851075618374];
        case{7,8} % 16 gauss points 
              gaussvalues = [1/3 1/3 0.144315607677787;...
                             0.081414823414554 0.459292588292723 0.095091634267285;...
                             0.459292588292723 0.081414823414554  0.095091634267285;...
                             0.459292588292723 0.459292588292723 0.095091634267285;...
                             0.658861384496480 0.170569307751760 0.103217370534718;...
                             0.170569307751760 0.658861384496480 0.103217370534718;...
                             0.170569307751760 0.170569307751760 0.103217370534718;...
                             0.898905543365938 0.050547228317031 0.032458497623198;...
                             0.050547228317031 0.898905543365938 0.032458497623198;...
                             0.050547228317031 0.050547228317031 0.032458497623198;...
                             0.008394777409958 0.263112829634638 0.027230314174435;...
                             0.728492392955404 0.008394777409958 0.027230314174435;...
                             0.008394777409958 0.728492392955404 0.027230314174435;...
                             0.728492392955404 0.263112829634638 0.027230314174435;...
                             0.263112829634638 0.008394777409958 0.027230314174435;...
                             0.263112829634638 0.728492392955404 0.027230314174435];
                                                       
        otherwise
            error(['ERROR: triangle quadrature not implemented for this number: ',int2str(surfacegaussorder)])
    end
    n2 = size(gaussvalues,1);    
end

for ii = 1:length(iv3)
    planenumber = iv3(ii);
    
    corners_x = planedata.corners(planedata.planecorners(planenumber,1:3),1);
    corners_y = planedata.corners(planedata.planecorners(planenumber,1:3),2);
    corners_z = planedata.corners(planedata.planecorners(planenumber,1:3),3);
    
    sidelength1 = norm( [corners_x(2)-corners_x(1) corners_y(2)-corners_y(1) corners_z(2)-corners_z(1)]); 
    sidelength2 = norm( [corners_x(3)-corners_x(2) corners_y(3)-corners_y(2) corners_z(3)-corners_z(2)]); 
    sidelength3 = norm( [corners_x(1)-corners_x(3) corners_y(1)-corners_y(3) corners_z(1)-corners_z(3)]); 
    semiperimeter = (sidelength1 + sidelength2 + sidelength3)/2;
    
    % Heron's area formula
    
    planearea = sqrt( semiperimeter*(semiperimeter-sidelength1)*(semiperimeter-sidelength2)*(semiperimeter-sidelength3) );

    % Convert the normalized triangle coordinates to the 3D coordinates
    
    psi=[gaussvalues(:,1) gaussvalues(:,2) 1-gaussvalues(:,1)-gaussvalues(:,2)];
    
    xgausspoints = psi*corners_x;
    ygausspoints = psi*corners_y;
    zgausspoints = psi*corners_z;

    recnvecs_oneplane = planedata.planeeqs(planenumber,1:3);
    recnvecs_oneplane = recnvecs_oneplane(ones(n2,1),:);
    surfacerecnvecs = [surfacerecnvecs;recnvecs_oneplane];
    
    recs_oneplane = [xgausspoints ygausspoints zgausspoints];
    recs_oneplane = recs_oneplane + recnvecs_oneplane*distancefromsurf;
    surfacerecs = [surfacerecs;recs_oneplane];

    recweights_oneplane = gaussvalues(:,3)*planearea;
    surfacerecweights = [surfacerecweights;recweights_oneplane];
    
end









