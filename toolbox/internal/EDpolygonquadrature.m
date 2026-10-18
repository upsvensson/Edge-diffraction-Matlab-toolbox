function [points,weights,area] = EDpolygonquadrature(corners,gaussorder)
% EDpolygonquadrature - Quadrature points and weights for a plane polygon.
% Gives the quadrature points across one plane, convex polygon, and the
% weights, which sum to 1. An integral over the polygon is then
%   int f dS  =  area * sum( weights .* f(points) )
% The rule depends on the number of corners:
%   3 corners   A symmetric triangle rule (Dunavant), of degree gaussorder
%               (0 to 8). The number of points is 1,3,4,6,7,12,16,16 for
%               gaussorder = 0/1,2,3,4,5,6,7,8.
%   4 corners   A tensor Gauss-Legendre rule, gaussorder*gaussorder
%               points, mapped bilinearly onto the quadrilateral. The
%               weights are the products of the 1D weights, which is exact
%               for parallelograms but does not include the varying
%               Jacobian of a general quadrilateral.
%   >4 corners  The polygon is divided into triangles, with one corner
%               at the center point (the mean of the corners), and the
%               triangle rule is used on each triangle. The polygon must
%               be convex.
% The area is computed with Heron's formula for triangles, Bretschneider's
% formula for quadrilaterals, and as the sum of the triangle areas for
% polygons with more corners.
%
% Input parameters:
%   corners         Matrix, [ncorners,3], of the corner coordinates, in
%                   order around the polygon (either direction).
%   gaussorder      The quadrature order, an integer >= 0. For 4 corners,
%                   the value 0 is treated as 1.
%
% Output parameters:
%   points          Matrix, [npoints,3], of the quadrature points.
%   weights         Vector, [npoints,1], of the weights. They sum to 1.
%   area            The area of the polygon, in m2.
%
% References for the triangle quadrature:
% 1. D.A. Dunavant, "High degree efficient symmetrical Gaussian quadrature
%    rules for the triangle," Int. J. Num. Meth. Eng. 21, pp. 1129-1148
%    (1985).
% 2. P. M. Juhl, "The boundary element method for sound field calculations,"
%    Report nr. 55, The acoustics laboratory, Technical University of
%    Denmark, 1993.
%
% See also EDcheckinputstructs, EDgensurfreceivers, lgwt
%
% Peter Svensson (peter.svensson@ntnu.no) 18 Oct. 2026
%
% [points,weights,area] = EDpolygonquadrature(corners,gaussorder);

% 18 Oct. 2026 First version (by Claude, prompted by PS). The code was
% moved from EDcheckinputstructs (piston sources) and EDgensurfreceivers
% (surface receivers), which had two copies of it. The points and weights
% are the same as before, except for one bug fix: the triangle rule of
% order 3 had weights which summed to 0.5 instead of 1. New: polygons
% with more than 4 corners.

[ncorners,ncolumns] = size(corners);
if ncorners < 3 || ncolumns ~= 3
    error('ERROR: EDpolygonquadrature needs a matrix [ncorners,3] with at least 3 corners')
end

if ncorners == 3
    gaussvalues = edtrianglerule(gaussorder);
    psi = [gaussvalues(:,1) gaussvalues(:,2) 1-gaussvalues(:,1)-gaussvalues(:,2)];
    points = [psi*corners(:,1) psi*corners(:,2) psi*corners(:,3)];
    weights = gaussvalues(:,3);
    area = edtrianglearea(corners);

elseif ncorners == 4
    nGL = max(gaussorder,1);
    n2 = nGL^2;
    [x,w] = lgwt(nGL,0,1);
    x = x(end:-1:1);
    xy = [repmat(x,nGL,1) reshape(x(:,ones(1,nGL)).',n2,1)];
    weights = prod([repmat(w,nGL,1) reshape(w(:,ones(1,nGL)).',n2,1)],2);

    c1 = corners(1,:);
    c2 = corners(2,:);
    c3 = corners(3,:);
    c4 = corners(4,:);
    xvec1 = c2 - c1;
    xvec2 = c3 - c4;
    startpointsx = c1(ones(n2,1),:) + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);
    endpointsx   = c4(ones(n2,1),:) + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
    points = startpointsx + (endpointsx - startpointsx).*xy(:,[2 2 2]);

    % Bretschneider's formula, with the side lengths a,b,c,d and the
    % diagonals p,q
    alen = norm(c2 - c1);
    blen = norm(c3 - c2);
    clen = norm(c3 - c4);
    dlen = norm(c4 - c1);
    plen = norm(c3 - c1);
    qlen = norm(c4 - c2);
    area = sqrt( 4*plen^2*qlen^2 - (blen^2 + dlen^2 - alen^2 - clen^2)^2 )/4;

else
    % Check that the polygon is plane and convex. The normal vector is
    % computed from the vector area (Newell's method).
    nextcorners = corners([2:ncorners 1],:);
    vectorarea = sum(cross(corners,nextcorners,2),1)/2;
    nvec = vectorarea/norm(vectorarea);
    polygonsize = max(max(corners,[],1) - min(corners,[],1));
    distfromplane = (corners - corners(ones(ncorners,1),:))*nvec.';
    if any(abs(distfromplane) > 1e-9*polygonsize)
        error('ERROR: The corners of a polygon with more than 4 corners must be in one plane')
    end
    prevcorners = corners([ncorners 1:ncorners-1],:);
    turnvecs = cross(corners - prevcorners,nextcorners - corners,2);
    if any(turnvecs*nvec.' < -1e-9*polygonsize^2)
        error('ERROR: A polygon with more than 4 corners must be convex')
    end

    midpoint = mean(corners,1);
    gaussvalues = edtrianglerule(gaussorder);
    psi = [gaussvalues(:,1) gaussvalues(:,2) 1-gaussvalues(:,1)-gaussvalues(:,2)];
    npertriangle = size(gaussvalues,1);
    points = zeros(npertriangle*ncorners,3);
    weights = zeros(npertriangle*ncorners,1);
    triangleareas = zeros(ncorners,1);
    for ii = 1:ncorners
        trianglecorners = [midpoint;corners(ii,:);nextcorners(ii,:)];
        iv = (ii-1)*npertriangle + (1:npertriangle);
        points(iv,:) = psi*trianglecorners;
        triangleareas(ii) = edtrianglearea(trianglecorners);
        weights(iv) = gaussvalues(:,3)*triangleareas(ii);
    end
    area = sum(triangleareas);
    weights = weights/area;
end

%--------------------------------------------------------------------------
function area = edtrianglearea(corners)
% Heron's formula, with the semiperimeter

sidelength1 = norm(corners(2,:) - corners(1,:));
sidelength2 = norm(corners(3,:) - corners(2,:));
sidelength3 = norm(corners(1,:) - corners(3,:));
semiperimeter = (sidelength1 + sidelength2 + sidelength3)/2;
area = sqrt( semiperimeter*(semiperimeter-sidelength1)*(semiperimeter-sidelength2)*(semiperimeter-sidelength3) );

%--------------------------------------------------------------------------
function gaussvalues = edtrianglerule(gaussorder)
% Symmetric triangle rules (Dunavant 1985). Each row: two of the three
% barycentric coordinates, and the weight. The weights sum to 1.

switch gaussorder
    case {0,1}
        gaussvalues = [1/3 1/3 1];
    case {2}
        gaussvalues = [2/3 1/6 1/3;...
                       1/6 2/3 1/3;
                       1/6 1/6 1/3];
    case {3}
        % NB! Until 18 Oct. 2026, the weights were -9/32 and 25/96, which
        % sum to 0.5 instead of 1.
        gaussvalues = [1/3 1/3 -27/48; ...
                       0.2 0.6 25/48; ...
                       0.2 0.2 25/48; ...
                       0.6 0.2 25/48];
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
    case {6}  % 12 gauss points
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
    case {7,8} % 16 gauss points
        gaussvalues = [1/3 1/3 0.144315607677787;...
                       0.081414823414554 0.459292588292723 0.095091634267285;...
                       0.459292588292723 0.081414823414554 0.095091634267285;...
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
        error(['ERROR: triangle quadrature not implemented for this number: ',int2str(gaussorder)])
end
