function xis = EDfindis(xsou,ivec,planeeqs)
% EDfindis - Returns the image source coordinates
% Returns the image source coordinates
% via mirroring one or N sources in a list of planes.
%
% Input parameters:
%   xsou		Matrix, [N,3], of the source coordinates.
%   ivec		List, [nplanes,1], of plane numbers that the source
%               should be mirrored in. If N ~= 1, then nplanes must
%               be = N.
%   planeeqs	Matrix [nplanes,4] of the plane equations.
%
% Output paramaters:
%   xis			Matrix [N,3] of all image sources to the source(s)
%               via all the planes specified in ivec and planeeqs.
%
% Uses no special subroutines
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2017
%
% xis = EDfindis(xsou,ivec,planeeqs);

% 6 July 2003 Functioning version
% 27 Nov. 2017 Copied to EDtoolbox from ESIE2toolbox. Removed the input
% parameters nsou and onesvec.

nsou = size(xsou,1);
onesvec = ones(1,3);

if nsou == 1
    Sbig = xsou(ones(length(ivec),1),:);
    t = planeeqs(ivec,4) - sum((planeeqs(ivec,1:3).').*(Sbig.')).';

    xis = Sbig + 2*t(:,onesvec).*planeeqs(ivec,1:3);
else
    t = planeeqs(ivec,4) - sum((planeeqs(ivec,1:3).').*(xsou.')).';

    xis = xsou + 2*t(:,onesvec).*planeeqs(ivec,1:3);
end

