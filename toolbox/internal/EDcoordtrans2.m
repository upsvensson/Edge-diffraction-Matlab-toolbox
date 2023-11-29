function [rs,thetas,zs,rr,thetar,zr] = EDcoordtrans2(xsou,xrec,xwedge,nvec1)
% EDcoordtrans2 - Transforms two sets of cartesian coordinates to edge-related cylindrical coordinates.
% The cyl. coord. system is defined so that:
%   - A z-axis is placed along the edge, from the given endpoint 1 to the given
%     endpoint 2.
%   - The origo of the cyl. syst. will be edge endpoint 1.
%   - The theta-angles of the cyl. coord. syst. will refer to the
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
%	xrec		Matrix, [n2,3] of cartesian coordinates of n2 other points. 
%	xwedge 		Matrix, [2,3], with the cartesian coordinates of the two 
%				wedge end points: [xw1 yw1 zw1;xw2 yw2 zw2].
%   nvec1       List, [1,3], with the normal vector of the reference plane
%               of the edge.
%
% Output parameters:
%	rs, thetas, zs		cyl. coord. of the points in xsou
%	rr, thetar, zr		cyl. coord. of the points in xrec
%
% Uses the subroutine EDcross
%
% Peter Svensson (svensson@iet.ntnu.no) 27 Nov. 2017
%
% [rs,thetas,zs,rr,thetar,zr] = EDcoordtrans2(xsou,xrec,xwedge,nvec1);

% 19 Jul 2010 Functioning version
% 27 Nov. 2017 Copied from ESIE2toolbox

xneworigo = xwedge(1,:);

xknown1 = xwedge(2,:) - xneworigo;
xknown1 = xknown1 / sqrt( sum( xknown1.^2 ));

%xknown3 = nvec1;

%xknown2 = cross(xknown3,xknown1);
%xknown2 = ESIE2cross(xknown3.',xknown1.').';

%xknown = [xknown1.' xknown2.' xknown3.'];
%xnew = [0 1 0;0 0 1;1 0 0];
% A = xnew*inv(xknown);

% A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' ESIE2cross(nvec1.',xknown1.') nvec1.']);
A = [0 1 0;0 0 1;1 0 0]/[xknown1.' EDcross(nvec1.',xknown1.') nvec1.'];

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

npoints = size(xrec,1);
if npoints >0
	xrec =    (A*( xrec.' - xneworigo(ones(npoints,1),:).' )).';
	rr = sqrt( sum(xrec(:,1:2).'.^2) ).';
	zr = xrec(:,3);
	thetar = zeros(npoints,1);
	iv = find(rr>0);
	if ~isempty(iv)	
		%%%thetar = 0*( xrec(:,2) == 0) + acos( xrec(:,1)./(rr+eps) ).*( xrec(:,2) ~= 0);
		thetar(iv) = real( acos( xrec(iv,1)./rr(iv) ).*( xrec(iv,2) ~= 0) );
		thetar(iv) = thetar(iv) + pi*( (xrec(iv,2)==0) & xrec(iv,1) < 0 );
		thetar(iv) = thetar(iv).*( xrec(iv,2) >=0 ) + (2*pi - thetar(iv)).*( xrec(iv,2) < 0 );
	end
else
	rr = []; thetar =[]; zr = [];
end
