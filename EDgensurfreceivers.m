function [surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
    surfacegaussorder,EDversionnumber,showtext);
% EDgensurfreceivers distributes receivers along all the faces/polygons of
% a polyhedral object. A Gauss-Legendre distribution is used. The function
% can handle general quadrilaterals and triangles.
% The quadrature order is given as an input parameter, and this value will
% be used for the largest polygon, and then scaled down for smaller
% polygons.
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
%
% Uses the function lgwt from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 8 March 2018

% 1 Mar 2018 First version
% 8 Mar 2018 Changed the displacement of the surface receivers to 0.1mm

distancefromsurf = 1e-4;

nplanes = size(planedata.planecorners,1);
if any(planedata.ncornersperplanevec~=4)
    error('ERROR: Unfortunately, only quadrilateral surfaces ',...
          'can be handled by this version of EDgensurfreceivers');
end

planeextensions = planedata.maxvals - planedata.minvals;
planeextensions = max(planeextensions.').';
planerelsize = planeextensions/max(planeextensions);

gaussorderperplane = ceil(planerelsize*surfacegaussorder);

surfacerecs = [];
surfacerecnvecs = [];
surfacerecweights = [];

for ii = surfacegaussorder:-1:1
   n2 = ii^2;
   iv = find(gaussorderperplane == ii);
   [x,w] = lgwt(ii,0,1);
   x = x(end:-1:1);
   xy = [repmat(x,ii,1) reshape(x(:,ones(1,ii)).',n2,1)];
   ww = prod([repmat(w,ii,1) reshape(w(:,ones(1,ii)).',n2,1)],2);   
   for jj = 1:length(iv)
       planenumber = iv(jj);
       xvec = planedata.corners(planedata.planecorners(planenumber,2),:) - planedata.corners(planedata.planecorners(planenumber,1),:);
       yvec = planedata.corners(planedata.planecorners(planenumber,4),:) - planedata.corners(planedata.planecorners(planenumber,1),:);
       crossvec1 = planedata.corners(planedata.planecorners(planenumber,3),:) - planedata.corners(planedata.planecorners(planenumber,1),:);
       crossvec2 = planedata.corners(planedata.planecorners(planenumber,4),:) - planedata.corners(planedata.planecorners(planenumber,2),:);
       if norm(crossvec1) ~= norm(crossvec2)
          error('ERROR: general quadrilaterals have not been implemented yet') 
       end
       length1 = norm(xvec);
       length2 = norm(yvec);
       planearea = length1*length2;
       
       recnvecs_oneplane = planedata.planeeqs(planenumber,1:3);
       recnvecs_oneplane = recnvecs_oneplane(ones(n2,1),:);
       recweights_oneplane = ww*planearea;

       startpoint = planedata.corners(planedata.planecorners(planenumber,1),:);
       startpoint = startpoint(ones(n2,1),:);
       recs_oneplane = startpoint + xvec(ones(n2,1),:).*xy(:,[1 1 1]) + yvec(ones(n2,1),:).*xy(:,[2 2 2]);
       recs_oneplane = recs_oneplane + recnvecs_oneplane*distancefromsurf;
       
       surfacerecs = [surfacerecs;recs_oneplane];
       surfacerecnvecs = [surfacerecnvecs;recnvecs_oneplane];
       surfacerecweights = [surfacerecweights;recweights_oneplane];
   end
end










