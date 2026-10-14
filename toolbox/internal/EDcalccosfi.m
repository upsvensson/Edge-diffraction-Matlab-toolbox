function [Distances,Cosfi] = EDcalccosfi(coordreceivers,coordsurfacesources,nvecsurfacesources)
% EDcalccosfi - Computes the distances and cos(phi) from surface points to receivers.
% Gives a distance matrix, and a matrix with the cosine of the angle
% between the normal vector of each surface point and the direction to
% each receiver. Used by EDpropagateESIEBEM.
%
% Input parameters:
%   coordreceivers       [nreceivers,3]
%   coordsurfacesources  [nsurfacesources,3] the surface points
%   nvecsurfacesources   [nsurfacesources,3] the normal vectors at the
%                        surface points
%
% Output parameters:
%   Distances            [nsurfacesources,nreceivers]
%   Cosfi                [nsurfacesources,nreceivers]
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [Distances,Cosfi] = EDcalccosfi(coordreceivers,coordsurfacesources,...
%    nvecsurfacesources);

% 970901 Functioning version
% 971002 Included the input parameter dipolevec, giving the possibility
%        to use dipole sources.
% 971110 Took care of the distance = 0 possibility.
% 971111 Speeded up calculations.
% 040606 Updated to Matlab 6.5
% 13 Dec. 2017 Copied from RADtoolbox to GEOtoolbox
% 18 Jan 2018 Copied to EDtoolbox
% 2 Oct. 2023 Rewrote it to clarify the code and the use.
%
% Uses the function EDcalcdist
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

nsurfacesources = size(coordsurfacesources,1);

nreceivers = size(coordreceivers,1);

% Distances will get size [nsurfacesources,nreceivers]:
Distances = EDcalcdist(coordsurfacesources,coordreceivers);

Cosfi = zeros(nsurfacesources,nreceivers);

for ii = 1:nsurfacesources
	coordsinglesource = coordsurfacesources(ii,:);
	nvecsinglesource = nvecsurfacesources(ii,:);
	xdiff = coordreceivers - coordsinglesource(ones(nreceivers,1),:);
	Cosfi(ii,:) = sum( (xdiff.*nvecsinglesource(ones(nreceivers,1),:)).' )./(Distances(ii,:) + eps);
end
