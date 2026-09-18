function Hq = EDmatvecH(Q,Hsubmatrixdata,Hsubdatalists,ivuselists)
% EDmatvecH mutliplies a sparse matrix H with a vertical vector Q (H*Q).
% The sparse matrix is stored with the format specified with Hsubmatrixdata.
% This is used both in the iterative solution of the integral equation and
% in the GMRES version.
%
% Input parameters:
%   Q               Vertical vector of size [N,1];
%   Hsubmatrixdata  Struct with the submatrix structure information
%   Hsubdatalists   Cell array with the actual submatrixdata
%   ivuselists      A vector which points at the elements that should be used for each submatrix.
%
% Output parameters:
%   Hq              Vertical vector of size [N,1]
%
% Peter Svensson  11 Sep. 2026

% 11 Sep. 2026 First version, by Claude, after prompting by PS.

N = double(Hsubmatrixdata.bigmatrixendnums(end));
nsub = size(Hsubmatrixdata.listofsubmatrices,1);
reftoshortlist = Hsubmatrixdata.reftoshortlist(:);

groupstart = [1; find(diff(reftoshortlist)~=0)+1];
groupend   = [groupstart(2:end)-1; nsub];
ngroups = length(groupstart);

Hq = zeros(N,1,'like',1i);

for g = 1:ngroups
    nn0 = groupstart(g);
    nnrange = groupstart(g):groupend(g);
    m = numel(nnrange);

    ii_list = Hsubmatrixdata.listofsubmatrices(nnrange,2);
    jj_list = Hsubmatrixdata.listofsubmatrices(nnrange,3);
    outstart = double(Hsubmatrixdata.bigmatrixstartnums(ii_list)) - 1;   % [m,1]
    instart  = double(Hsubmatrixdata.bigmatrixstartnums(jj_list)) - 1;   % [m,1]

    ivuse    = double(ivuselists{nn0});   % cast: avoid uint64+double-vector broadcast error
    Hsubdata = Hsubdatalists{nn0};
    nHsub1   = double(Hsubmatrixdata.bigmatrixendnums(ii_list(1))) - ...
               double(Hsubmatrixdata.bigmatrixstartnums(ii_list(1))) + 1;
    row_Hsub = rem(ivuse-1,nHsub1) + 1;
    col_Hsub = (ivuse-row_Hsub)/nHsub1 + 1;

    % Expand this group's indices/values fresh -- these temporaries are
    % discarded after the accumarray call, never accumulated across groups.
    rowidx = reshape(row_Hsub(:) + outstart(:).', [], 1);
    colidx = reshape(col_Hsub(:) + instart(:).', [], 1);
    vals   = repmat(Hsubdata(:), m, 1);

    contrib = vals .* Q(colidx);
    Hq = Hq + accumarray(rowidx, contrib, [N,1]);
end
