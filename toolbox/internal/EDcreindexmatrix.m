function ivmatrix = EDcreindexmatrix(ndivvec)
% EDcreindeixmatrix creates a matrix with index numbers.
%
% Input parameters:
%   ndivvec     A matrix, [1,specorder], with the maximum counter
%               number for each dimension.
%
% Output parameters:
%   ivmatrix    A matrix, [prod(ndivvec),specorder] of all possible
%               combinations of the integers
%
% Peter Svensson (peter.svensson@ntnu.no) 16 Mar 2018
%
% ivmatrix = EDcreindexmatrix(ndivvec)

% 5 May 2005 Functioning version in EDB1toolbox
% 16 Mar 2018 Copied to EDtoolbox

n = length(ndivvec);

maxval = max(ndivvec);
if maxval > 2^32
    error('ERROR: This version of EDcreindexmatrix can not create such large matrics')    
end

if n == 2
	iv1 = uint32([1:ndivvec(1)].');
	iv2 = uint32([1:ndivvec(2)]);
	iv1 = iv1(:,ones(1,ndivvec(2)));
	iv2 = iv2(ones(ndivvec(1),1),:);

    ivmatrix = [reshape(iv1.',prod(ndivvec),1) reshape(iv2.',prod(ndivvec),1)];
elseif n >= 3
    ivmatrix = EDcreindexmatrix(ndivvec(2:n));
    ivmatrix = repmat(ivmatrix,[ndivvec(1),1]);
    ivfirstcol = uint32([1:ndivvec(1)].');
    ivfirstcol = ivfirstcol(:,ones(1,prod(ndivvec(2:n))));
    ivfirstcol = reshape(ivfirstcol.',prod(ndivvec),1);
    ivmatrix = [ivfirstcol ivmatrix];
end
