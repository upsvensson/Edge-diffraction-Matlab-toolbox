function [surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
    surfacegaussorder,EDversionnumber)
% EDgensurfreceivers - Distributes receivers along all the faces/polygons of
% a polyhedral object. The quadrature points and weights for each face are
% given by EDpolygonquadrature: a Gauss-Legendre distribution for
% quadrilaterals, a symmetric triangle rule for triangles, and the triangle
% rule on a division into triangles for faces with more corners (which
% must be convex).
% The quadrature order is given as an input parameter, and this value will
% be used for the largest quadrilateral, and then scaled down for smaller
% ones. For the other faces, the same quadrature order is used for all
% (at most 8).
% 
% Input parameters:
%   planedata           Struct
%   surfacegaussorder   The quadratureorder that will be used for the largest
%                       plane/face/polygon of the model. A value of 5
%                       implies that 5*5 receivers will be generated across
%                       the largest surface.
%   EDversionnumber
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
% See also EDpolygonquadrature
%
% Peter Svensson (peter.svensson@ntnu.no) 18 Oct. 2026
%
% [surfacerecs,surfacerecnvecs,...
%    surfacerecweights] = EDgensurfreceivers(planedata,surfacegaussorder,...
%    EDversionnumber);

% 1 Mar 2018 First version
% 8 Mar 2018 Changed the displacement of the surface receivers to 0.1mm
% 4 Apr 2018 Introduced the triangle element quadrature, using gaussTRINEW
% from OpenBEM
% 5 Apr 2018 Small error if no triangular elements. Implemented quadrature
% from the Duvanant paper instead of the gaussTRINEW function.
% 10 Apr 2018 Added some quadrature values
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).
% 18 Oct. 2026 The quadrature was moved to the new function
% EDpolygonquadrature (by Claude, prompted by PS). The receivers and
% weights are unchanged, except for triangles with surfacegaussorder = 3,
% where the weights were half as large as they should be. Faces with more than 4 corners are now handled
% (they come last in the output).

% distancefromsurf = 1e-4;
distancefromsurf = 1e-3;

nplanes = size(planedata.planecorners,1);

iv3 = find(planedata.ncornersperplanevec == 3);
iv5 = find(planedata.ncornersperplanevec > 4);

planeextensions = planedata.maxvals - planedata.minvals;
planeextensions = max(planeextensions.').';
planerelsize = planeextensions/max(planeextensions);

gaussorderperplane = ceil(planerelsize*surfacegaussorder);
if ~isempty(iv3) || ~isempty(iv5)
    if surfacegaussorder >8 
       disp('   WARNING: surfacegaussorder is decreased to 8 (for the triangular elements)') 
       surfacegaussorder = 8;
    end
    gaussorderperplane([iv3(:);iv5(:)]) = 0;
end
    
surfacerecs = [];
surfacerecnvecs = [];
surfacerecweights = [];

% First all quadrilateral elements, from the highest quadrature order to
% the lowest, then all triangles, and last the faces with more corners.

for ii = max(gaussorderperplane):-1:1
   iv = find(gaussorderperplane == ii);
   for jj = 1:length(iv)
       [surfacerecs,surfacerecnvecs,surfacerecweights] = edaddoneplane(planedata,...
           iv(jj),ii,distancefromsurf,surfacerecs,surfacerecnvecs,surfacerecweights);
   end
end

for ii = 1:length(iv3)
    [surfacerecs,surfacerecnvecs,surfacerecweights] = edaddoneplane(planedata,...
        iv3(ii),surfacegaussorder,distancefromsurf,surfacerecs,surfacerecnvecs,surfacerecweights);
end

for ii = 1:length(iv5)
    [surfacerecs,surfacerecnvecs,surfacerecweights] = edaddoneplane(planedata,...
        iv5(ii),surfacegaussorder,distancefromsurf,surfacerecs,surfacerecnvecs,surfacerecweights);
end

%--------------------------------------------------------------------------
function [surfacerecs,surfacerecnvecs,surfacerecweights] = edaddoneplane(planedata,...
    planenumber,gaussorder,distancefromsurf,surfacerecs,surfacerecnvecs,surfacerecweights)
% Adds the receivers of one plane, moved the distance distancefromsurf
% out from the surface.

ncorners = planedata.ncornersperplanevec(planenumber);
corners = planedata.corners(planedata.planecorners(planenumber,1:ncorners),:);
[recs_oneplane,ww,planearea] = EDpolygonquadrature(corners,gaussorder);
n2 = size(recs_oneplane,1);

recnvecs_oneplane = planedata.planeeqs(planenumber,1:3);
recnvecs_oneplane = recnvecs_oneplane(ones(n2,1),:);
recs_oneplane = recs_oneplane + recnvecs_oneplane*distancefromsurf;

surfacerecs = [surfacerecs;recs_oneplane];
surfacerecnvecs = [surfacerecnvecs;recnvecs_oneplane];
surfacerecweights = [surfacerecweights;ww*planearea];
