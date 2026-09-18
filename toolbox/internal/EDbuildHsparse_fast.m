function Hsp = EDbuildHsparse_fast(Hsubmatrixdata,Hsubdatalists,ivuselists)
% EDbuildHsparse_fast - faster/lower-memory version of EDbuildHsparse.
%
% Two changes vs. the original EDbuildHsparse:
%
%  1) Preallocates the triplet arrays (rows,cols,vals) up front instead of
%     growing cell arrays and concatenating with cell2mat. Avoids the
%     transient ~2x peak-memory spike that cell2mat causes when the
%     cell contents and the freshly-concatenated array briefly coexist.
%
%  2) Submatrices are already sorted by shortlist number in
%     Hsubmatrixdata.listofsubmatrices (see EDinteg_submatrixstructure,
%     "[~,reftosortvec] = sort(reftoshortlist)"). All submatrices sharing
%     one shortlist number have IDENTICAL row_Hsub/col_Hsub/Hsubdata and
%     differ only by their (outstart,instart) block offset. Rather than
%     looping submatrix by submatrix (nsub iterations), this loops once
%     per shortlist GROUP (nuniquesubmatrices iterations) and fills all
%     of that group's triplets in one vectorized step.
%
% Sanity check before trusting this in production: compare against
% EDbuildHsparse on a small case, e.g.
%   norm(Hsp_fast - Hsp_simple, 'fro') / norm(Hsp_simple,'fro')
% should be exactly (or to roundoff) zero, since no values are summed
% differently -- only the construction order changes.
%
% Peter Svensson 3 Sep. 2026
%
% Hsp = EDbuildHsparse_fast(Hsubmatrixdata,Hsubdatalists,ivuselists)

% 3 Sep. 2026 First version, by Claude, prompted by PS

N = double(Hsubmatrixdata.bigmatrixendnums(end));
nsub = size(Hsubmatrixdata.listofsubmatrices,1);
reftoshortlist = Hsubmatrixdata.reftoshortlist(:);

% ---- Group boundaries. nn's are already sorted by shortlist number, so
% groups are contiguous runs of equal reftoshortlist values.
groupstart = [1; find(diff(reftoshortlist)~=0)+1];
groupend   = [groupstart(2:end)-1; nsub];
ngroups = length(groupstart);

% ---- Pass 1 (cheap): find each group's per-submatrix triplet count, to
% preallocate exactly. Only touches ivuselists (index lists), not the
% (larger) Hsubdatalists value arrays.
grouplen = zeros(ngroups,1);
for g = 1:ngroups
    grouplen(g) = numel(ivuselists{groupstart(g)});
end
groupcount = groupend - groupstart + 1;      % submatrices per group
totalnnz = sum(grouplen.*groupcount);

rows = zeros(totalnnz,1);
cols = zeros(totalnnz,1);
vals = zeros(totalnnz,1,'like',1i);          % complex double

% ---- Pass 2: fill, one vectorized block per shortlist group
ptr = 0;
for g = 1:ngroups
    nn0 = groupstart(g);
    nnrange = groupstart(g):groupend(g);
    m = numel(nnrange);

    ii_list = Hsubmatrixdata.listofsubmatrices(nnrange,2);
    jj_list = Hsubmatrixdata.listofsubmatrices(nnrange,3);
    outstart = double(Hsubmatrixdata.bigmatrixstartnums(ii_list)) - 1;   % [m,1]
    instart  = double(Hsubmatrixdata.bigmatrixstartnums(jj_list)) - 1;   % [m,1]

    ivuse    = double(ivuselists{nn0});   % cast: uint64 can't broadcast against a double vector below
    Hsubdata = Hsubdatalists{nn0};
    nHsub1   = double(Hsubmatrixdata.bigmatrixendnums(ii_list(1))) - ...
               double(Hsubmatrixdata.bigmatrixstartnums(ii_list(1))) + 1;
    row_Hsub = rem(ivuse-1,nHsub1) + 1;
    col_Hsub = (ivuse-row_Hsub)/nHsub1 + 1;

    len = numel(row_Hsub);
    idx = ptr + (1:len*m);

    % Broadcasting: [len,1] + [1,m] -> [len,m], then column-major reshape
    % to [len*m,1]. repmat(...,m,1) below uses the same len-fastest
    % ordering, so rows/cols/vals stay aligned.
    rows(idx) = reshape(row_Hsub(:) + outstart(:).', len*m, 1);
    cols(idx) = reshape(col_Hsub(:) + instart(:).', len*m, 1);
    vals(idx) = repmat(Hsubdata(:), m, 1);

    ptr = ptr + len*m;
end

%Hsp = sparse(rows,cols,vals,N,N);
Hsp = fsparse(rows,cols,vals,[N,N]);
Hsp = fsparse(rows,cols,real(vals),[N,N]);
Hspi = fsparse(rows,cols,imag(vals),[N,N]);
