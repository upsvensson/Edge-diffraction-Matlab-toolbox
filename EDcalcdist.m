function distances = EDcalcdist(x1,x2);
% EDcalcdist - Gives distances between sets of points in 3D space.
% EDcalcdist returns a matrix of distances from all n1 points
% in x1 to all n2 points in x2.
%
% Input parameters:
%   x1      Matrix, [n1,3], with the coordinates of n1 points.
%   x2      Matrix, [n2,3], with the coordinates of n2 points.
%
% Output parameters:
%   distances   Matrix, [n1,n2], with the distances from all
%               points in x1 to all points in x2.
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
% Peter Svensson (peter.svensson@ntnu.no) 28 Nov. 2017
%
% distances = EDcalcdist(x1,x2);

% 10 Oct. 2000 First version
% 28 Nov. 2017 Copied to EDtoolbox

[n11,n12] = size(x1);
[n21,n22] = size(x2);

if n12 ~= 3 || n22 ~= 3
	disp('Coordinates must be specified as [x y z]')
else
	if n11 == 1
		x1 = x1(ones(n21,1),:);
		d1 = x1-x2;
		distances = sqrt(sum( (d1.').^2 ));
	elseif n21 == 1
		x2 = x2(ones(n11,1),:);
		d1 = x1-x2;
		distances = sqrt(sum( (d1.').^2 )).';
	else
		distances = zeros(n11,n21);
		for ii = 1:n21
			x2point = x2(ii,:);
			x2point = x2point(ones(n11,1),:);
			d1 = x1 - x2point;
			distances(:,ii) = sqrt(sum( (d1.').^2 )).';
		end
	end
end
