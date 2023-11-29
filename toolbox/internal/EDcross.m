function c = EDcross(a,b)
% EDcross - Stripped down version of Matlab's built in function cross.
%
% Input parameters:
%   a,b     Two vectors, of size (3,1). Size (1,3) will give an error
%           message.
% 
% Output parameter:
%   c       One vector, of size (3,1)
%
%  c = EDcross(a,b)

%CROSS  Vector cross product.
%   C = CROSS(A,B) returns the cross product of the vectors
%   A and B.  That is, C = A x B.  A and B must be 3 element
%   vectors.
%
%   C = CROSS(A,B) returns the cross product of A and B along the
%   first dimension of length 3.
%

%   Clay M. Thompson
%   updated 12-21-94, Denise Chen
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 5.18 $  $Date: 2002/04/09 00:29:46 $

% Special case: A and B are vectors
%rowvec = 0;
%if ndims(a)==2 & ndims(b)==2
%    if size(a,1)==1, a = a(:); rowvec = 1; end
%    if size(b,1)==1, b = b(:); rowvec = 1; end
%end;

% Calculate cross product
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
c = reshape(c,size(a));
