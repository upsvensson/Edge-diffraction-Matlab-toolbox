function distances = EDcalcdist(x1,x2)
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
% Peter Svensson (peter.svensson@ntnu.no) 3 July 2019
%
% distances = EDcalcdist(x1,x2);

% 10 Oct. 2000 First version
% 28 Nov. 2017 Copied to EDtoolbox
% 3 July 2019 Implemented a speed-up if x2 has many more points than x1.

[n11,n12] = size(x1);
[n21,n22] = size(x2);

if n12 ~= 3 || n22 ~= 3
	disp('Coordinates must be specified as [x y z]')
else
	if n11 == 1
		distances = sqrt(sum( ((x1(ones(n21,1),:)-x2).').^2 ));
	elseif n21 == 1
		distances = sqrt(sum( ((x1-x2(ones(n11,1),:)).').^2 )).';
	else
		distances = zeros(n11,n21);
        if n21 < n11
            onesvec = ones(n11,1);
       		for ii = 1:n21
                x2point = x2(ii,:);
                distances(:,ii) = sqrt(sum( ((x1 - x2point(onesvec,:)).').^2 )).';
            end
        else
            onesvec = ones(n21,1);
            for ii = 1:n11
                x1point = x1(ii,:);
                distances(ii,:) = sqrt(sum( ((x2 - x1point(onesvec,:)).').^2 ));                
            end
        end
	end
end
