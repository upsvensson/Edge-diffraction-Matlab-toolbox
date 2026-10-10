function Hsp = EDbuildHsparse(Hsubmatrixdata,Hsubdatalists,ivuselists)
% EDbuildHsparse - Assemble the full sparse edge-source H-matrix from the
% submatrix lists produced by EDcreateHsubmatrices.
%
% This reproduces, index-for-index, the block bookkeeping used in the
% accumarray-based mat-vec loop in EDintegralequation_convex_tf.m
% (the "shortlistprev / Hsubdata / ivuse" reuse logic included), but
% instead of multiplying by Q it just records global (row,col,val)
% triplets. Because each (to-edgepair,from-edgepair) combination maps to
% exactly one submatrix, no global (row,col) pair is ever produced twice,
% so sparse() reproduces the accumarray result exactly.
%
% Input parameters:
%   Hsubmatrixdata     Struct from EDinteg_submatrixstructure
%   Hsubdatalists      Cell array of submatrix values, from EDcreateHsubmatrices
%   ivuselists         Cell array of submatrix index maps, from EDcreateHsubmatrices
%
% Output parameters:
%   Hsp                 Sparse matrix, size [N,N], N = Qvecendnums(end)
%
% Peter Svensson 9 Oct. 2026
%
% Hsp = EDbuildHsparse(Hsubmatrixdata,Hsubdatalists,ivuselists);

% 2 Sep. 2026 FIrst version, created by Claude, prompted by Peter Svensson
% 9 Oct. 2026 Code cleanup, changed bigmatrix... to Qvec...

nQ = double(Hsubmatrixdata.Qvecendnums(end));
nsub = size(Hsubmatrixdata.listofsubmatrices,1);

rowcell = cell(nsub,1);
colcell = cell(nsub,1);
valcell = cell(nsub,1);

shortlistprev = 0;
Hsubdata = [];
row_Hsub = [];
col_Hsub = [];

for nn = 1:nsub
    ii = Hsubmatrixdata.listofsubmatrices(nn,2);
    jj = Hsubmatrixdata.listofsubmatrices(nn,3);
    ivoutvec = Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii);
    ivinvec  = Hsubmatrixdata.Qvecstartnums(jj):Hsubmatrixdata.Qvecendnums(jj);

    shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);
    if shortlistnumber ~= shortlistprev
        Hsubdata = Hsubdatalists{nn};
        ivuse    = ivuselists{nn};
        nHsub1   = length(ivoutvec);
        row_Hsub = rem(ivuse-1,nHsub1) + 1;
        col_Hsub = (ivuse-row_Hsub)/nHsub1 + 1;
        shortlistprev = shortlistnumber;
    end

    rowcell{nn} = double(ivoutvec(1) - 1 + row_Hsub(:));
    colcell{nn} = double(ivinvec(1)  - 1 + col_Hsub(:));
    valcell{nn} = Hsubdata(:);
end

rows = cell2mat(rowcell);
cols = cell2mat(colcell);
vals = cell2mat(valcell);

Hsp = sparse(rows,cols,vals,nQ,nQ);
