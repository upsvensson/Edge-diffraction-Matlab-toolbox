function [rs,thetas,zs,A] = EDcoordtrans1(xsou,xwedge,nvec1,Bmatrix)
% EDcoordtrans1 - Transforms one set of cartesian coordinates to edge-related cylindrical coordinates.
% The cyl. coord. system is defined so that:
%     A z-axis is placed along the edge, from the given endpoint 1 to the given
%     endpoint 2.
%     The origo of the cyl. syst. will be edge endpoint 1.
%     The theta-angles of the cyl. coord. syst. will refer to the
%     reference plane of the edge.
% The ref. plane of the edge is described by its normal vector.
% NB! The order of the edge points is important!!
% The vector going from xwedge(1,:) to xwedge(2,:) must be
% oriented so that if the RH thumb is along this vector, the tips
% of the fingers "come out of" the open face of plane1, i.e. where nvec1
% is the normal vector.
%
% Input parameters:
%	xsou		Matrix, [n1,3] of cartesian coordinates of n1 points. 
%	xwedge 		Matrix, [2,3], with the cartesian coordinates of the two 
%				wedge end points: [xw1 yw1 zw1;xw2 yw2 zw2].
%   nvec1       List, [1,3], with the normal vector of the reference plane
%               of the edge.
%   Bmatrix     (optional) The row-shuffled version of the A matrix.
%
% Output parameters:
%	rs, thetas, zs		cyl. coord. of the points in xsou
%   A           The A-matrix used for computing the coordinate
%               transformation.
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2017
%
% [rs,thetas,zs,A] = EDcoordtrans1(xsou,xwedge,nvec1,Bmatrix);

% 18 Jan 2006   Functioning version
% 1 Dec. 2014   Added the output parameter A and optional input parameter B
% 27 Nov. 2017 Copied from ESIE2toolbox. Change: wrote out the code for
%              ESIE2cross.

if nargin == 3
    Bmatrixgiven = 0;
else
    Bmatrixgiven = 1;
end

npoints = size(xsou,1);
if npoints == 1
    xneworigo = xwedge(1,:);
    if Bmatrixgiven == 0
        xknown1 = xwedge(2,:) - xneworigo;
        xknown1 = xknown1 / norm(xknown1);

        % % % % % A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' ESIE2cross(nvec1.',xknown1.') nvec1.']);
        A = [nvec1(2)*xknown1(3)-nvec1(3)*xknown1(2) ; nvec1(3)*xknown1(1)-nvec1(1)*xknown1(3) ; nvec1(1)*xknown1(2)-nvec1(2)*xknown1(1)];
        A = inv([xknown1.' A nvec1.']);
        % % % % % A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' A nvec1.']);

        % % % % % xsou =    (A*( xsou.' - xneworigo(ones(npoints,1),:).' )).';
        xsou =    (A([2 3 1],:)*( xsou.' - xneworigo.' )).';
    else
        xsou =    (Bmatrix*( xsou.' - xneworigo.' )).';        
    end
    
    % % % % % rs = sqrt( sum(xsou(:,1:2).'.^2) ).';
    rs = norm(xsou(1:2));
    zs = xsou(:,3);
    thetas = 0;
    if rs > 0
        thetas = real( acos( xsou(1)./rs ).*( xsou(2) ~= 0) );
        thetas = thetas + pi*( (xsou(2)==0) & xsou(1) < 0 );
        thetas = thetas.*( xsou(2) >=0 ) + (2*pi - thetas).*( xsou(2) < 0 );
    end
    
else
    
    xneworigo = xwedge(1,:);
    if Bmatrixgiven == 0
        xknown1 = xwedge(2,:) - xneworigo;
        xknown1 = xknown1 / sqrt( sum( xknown1.^2 ));

        a = nvec1.';
        b = xknown1.';
        c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
            a(3,:).*b(1,:)-a(1,:).*b(3,:)
            a(1,:).*b(2,:)-a(2,:).*b(1,:)];
        c = reshape(c,size(a));

    %     A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' A nvec1.']);
%         A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' ESIE2cross(nvec1.',xknown1.') nvec1.']);
        A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' c nvec1.']);
    else
       A = Bmatrix; 
    end
    npoints = size(xsou,1);
    xsou =    (A*( xsou.' - xneworigo(ones(npoints,1),:).' )).';

    rs = sqrt( sum(xsou(:,1:2).'.^2) ).';
    zs = xsou(:,3);
    thetas = zeros(npoints,1);
    iv = find(rs>0);
    if ~isempty(iv)
        %%%thetas = 0*( xsou(:,2) == 0) + acos( xsou(:,1)./(rs+eps) ).*( xsou(:,2) ~= 0)
        thetas(iv) = real( acos( xsou(iv,1)./rs(iv) ).*( xsou(iv,2) ~= 0) );
        thetas(iv) = thetas(iv) + pi*( (xsou(iv,2)==0) & xsou(iv,1) < 0 );
        thetas(iv) = thetas(iv).*( xsou(iv,2) >=0 ) + (2*pi - thetas(iv)).*( xsou(iv,2) < 0 );
    end
end