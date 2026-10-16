function [Fgeom,pathlengths,propdiagnostics] = EDcalcpropagatematrix_block(edgedata,edgetoedgedata,...
    Hsubmatrixdata,vispartedgesfromr,rRmat,thetaRmat,zRmat,doesQsegmenthavevalues)
% EDcalcpropagatematrix_block - Propagation rows from the edge sources to a block of receivers.
% The same computation as EDcalcpropagatematrix with frequency = 0 and
% ESIEpropagationimplementation = 'accurate' (isthinhole = 0), but for
% many receivers at once: the loop is over the edge pairs, and for each
% edge pair all the receivers of the block are handled with array
% operations. The receiver-independent quantities (node positions,
% weights, mdist) are computed once per edge pair, and the near-singular
% (Lorentzian) corrections of all receivers are computed in one call.
% This is much faster than one call of EDcalcpropagatematrix per receiver,
% where Matlab's overhead per operation dominates (arrays of only n1*n2
% elements).
%
% For each receiver r and frequency (wavenumber k), the propagation row
% vector is Fgeom(r,:).*exp(-1i*k*(pathlengths(r,:)-Rstart)).
%
% Input parameters:
%   edgedata, edgetoedgedata, Hsubmatrixdata   Structs (see EDdatastructs)
%   vispartedgesfromr   [nedges,nrec] the visibility of the edges from
%                       the receivers of the block (Rdata.vispartedgesfromr)
%   rRmat, thetaRmat, zRmat   [nedges,nrec] the cylindrical coordinates of
%                       the receivers of the block relative to each edge
%   doesQsegmenthavevalues    [nedgepairs,1] 1 for the edge pairs where the
%                       edge source amplitudes are not all zero
%
% Output parameters:
%   Fgeom           [nrec,nQ] the frequency-independent part of the
%                   propagation row vectors
%   pathlengths     [nrec,nQ] the path lengths l + m (0 where Fgeom is 0)
%   propdiagnostics struct with totals over the receivers of the block:
%                   .npairs, .npairs_eligible, .npairs_corrected (counted
%                   per receiver and edge pair), .nrows, .nrows_corrected,
%                   and .correctedpairs, [edge1 edge2 corrected_rows
%                   all_rows] summed over the receivers, one row per edge
%                   pair with corrections
%
% See also EDcalcpropagatematrix, EDintegralequation_convex_tf
%
% Peter Svensson (peter.svensson@ntnu.no) 16 Oct. 2026
%
% [Fgeom,pathlengths,propdiagnostics] = EDcalcpropagatematrix_block(edgedata,...
%    edgetoedgedata,Hsubmatrixdata,vispartedgesfromr,rRmat,thetaRmat,zRmat,...
%    doesQsegmenthavevalues);

% 16 Oct. 2026 First version (by Claude, prompted by PS). The local
% functions gaussLegendreNodes, legendreVandermonde, glStandardCached,
% neumannQupward and modifiedWeightsLorentzianVec are copies of those in
% EDcalcpropagatematrix.

nQ = double(Hsubmatrixdata.Qvecendnums(end));
nrec = size(vispartedgesfromr,2);
Fgeom = zeros(nrec,nQ);
pathlengths = zeros(nrec,nQ);
refto = edgetoedgedata.reftoshortlistE;

alledgepairsthin = prod(double(Hsubmatrixdata.isthinplaneedgepair));
noedgepairsthin  = prod(1-double(Hsubmatrixdata.isthinplaneedgepair));

propdiagnostics = struct('npairs',0,'npairs_eligible',0,'npairs_corrected',0, ...
    'nrows',0,'nrows_corrected',0,'correctedpairs',zeros(0,4));
correctedpairs = zeros(size(Hsubmatrixdata.edgepairlist,1),4);
ncorrectedpairs = 0;

for ii = 1:size(Hsubmatrixdata.edgepairlist,1)
    if doesQsegmenthavevalues(ii) ~= 1
        continue
    end
    edge2 = Hsubmatrixdata.edgepairlist(ii,1);
    edge1 = Hsubmatrixdata.edgepairlist(ii,2);
    % The receivers of the block that see edge2
    rr = find(vispartedgesfromr(edge2,:) > 0);
    if isempty(rr)
        continue
    end
    nr = numel(rr);

    n1 = double(Hsubmatrixdata.nedgeelems(edge1));
    n2 = double(Hsubmatrixdata.nedgeelems(edge2));
    len2 = edgedata.edgelengthvec(edge2);
    if alledgepairsthin == 1
        thinplaneboost = 2;
    elseif noedgepairsthin == 1
        thinplaneboost = 1;
    else
        thinplaneboost = double(Hsubmatrixdata.isthinplaneedgepair(ii)+1);
    end
    ny = pi/(2*pi-edgedata.closwedangvec(edge2));

    % Receiver-independent quantities (computed for every edge pair; see
    % the bug fix of 15 Oct. 2026 in EDcalcpropagatematrix)
    n2vec = full(Hsubmatrixdata.quadraturematrix_pos(1:n2,n2));
    weightvec2 = full(Hsubmatrixdata.quadraturematrix_weights(1:n2,n2));
    dzvec2 = weightvec2*edgedata.edgelengthvec(edge2);
    ze2_re2 = (0 + n2vec*len2);
    n1vec = full(Hsubmatrixdata.quadraturematrix_pos(1:n1,n1));
    weightvec1 = full(Hsubmatrixdata.quadraturematrix_weights(1:n1,n1));
    dzvec1 = weightvec1*edgedata.edgelengthvec(edge1);

    ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );
    re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );

    thetain = edgetoedgedata.thetae1sho(refto(edge1,edge2));
    thetain_end = edgetoedgedata.thetae2sho(refto(edge1,edge2));
    if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_end) < 1e-10
        acrossface_in = 1;
    else
        acrossface_in = 0;
        xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
        xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
        yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
        yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
        xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
        ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);
        [re1_re2,thetae1_re2,ze1_re2] = EDcoordtrans1([xe1_re2 ye1_re2 ze1_re2],[0 0 0;0 0 len2],[0 1 0]);
        re1_re2 = re1_re2(:); thetae1_re2 = thetae1_re2(:); ze1_re2 = ze1_re2(:);
    end

    % Index vectors: element (i,j) of the n1-by-n2 block is number
    % i + (j-1)*n1, as in EDcalcpropagatematrix.
    n2vertmat = reshape(repmat(1:n2,n1,1),n1*n2,1);
    n1vertmat = repmat((1:n1).',n2,1);

    mdist = sqrt( (ze1_re2(n1vertmat) - ze2_re2(n2vertmat)).^2 + re1_re2(n1vertmat).^2 );   % [n1n2 x 1]

    % Receiver-dependent quantities, [n1n2 x nr]
    rR = rRmat(edge2,rr);           % [1 x nr]
    thetaout = thetaRmat(edge2,rr);
    zR = zRmat(edge2,rr);
    ldist = sqrt( (zR - ze2_re2(n2vertmat)).^2 + rR.^2 );
    pathlengths(rr,Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) = (ldist + mdist).';

    ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
    ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
    ch = ( ch + 1./ch)/2;

    % 1/(l*m), the frequency-independent part of the kernel
    kernelfactor = 1./ldist./mdist;

    if acrossface_in == 1
        % As EDcalcpropagatematrix_accuratebeta in EDcalcpropagatematrix
        s1 = [1 1];  s2 = [1 -1];  termboost = 2;
        betaterms = {sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))), ...
                     sin(ny*(pi + thetain - thetaout  ))./(ch - cos(ny*(pi + thetain - thetaout   )))};

        R0vec = sqrt( (ze1_re2 - zR).^2 + (re1_re2 + rR).^2 );          % [n1 x nr]
        zavec = ( zR.*re1_re2 + ze1_re2.*rR ) ./ ( re1_re2 + rR );        % [n1 x nr]
        kappavec = ny^2*(re1_re2+rR).^4 ./ (2*re1_re2.^2.*rR.^2.*R0vec.^2); % [n1 x nr]

        if n2 > 1
            localspacing = max(diff(sort(ze2_re2)));
        else
            localspacing = len2;
        end
        zrelthreshold = 5;
        insiderange = (zavec >= -0.05*len2) & (zavec <= 1.05*len2);   % [n1 x nr]

        dzvec2_row = dzvec2(:).';
        ze2_row = ze2_re2(:).';
        patched = zeros(n1*nr,n2);          % rows: (i,r) with i fastest
        anyflag = false(n1,nr);
        for kk = 1:2
            phi_k = pi + thetain*s1(kk) + thetaout*s2(kk);          % [1 x nr]
            B1_k  = 2*sin(ny*phi_k/2).^2 ./ kappavec;                 % [n1 x nr]
            num_k = sin(ny*phi_k) ./ kappavec;                        % [n1 x nr]
            flag_k = insiderange & (sqrt(max(B1_k,0)) < zrelthreshold*localspacing);
            anyflag = anyflag | flag_k;

            % The exact term, [n1 x n2 x nr] -> [n1*nr x n2] (rows (i,r))
            termk = reshape(betaterms{kk},n1,n2,nr);
            contrib = reshape(permute(termk,[1 3 2]),n1*nr,n2);
            contrib = termboost*contrib.*dzvec2_row;

            flaggedrows = find(flag_k(:));
            if ~isempty(flaggedrows)
                [~,wnodes_std,Sstd,tfine_std,wfine_std,Pint_std] = glStandardCached(n2);
                zaf  = zavec(flaggedrows);
                B1f  = max(B1_k(flaggedrows),0);
                numf = num_k(flaggedrows);
                xa   = 2*zaf/len2 - 1;
                Bhat = 4*B1f/len2^2;
                WT = modifiedWeightsLorentzianVec(Sstd,wnodes_std,tfine_std,wfine_std,Pint_std,xa,Bhat);   % [n2 x nf]
                wt_physical = termboost*(2/len2) * (numf .* WT.');           % [nf x n2]
                lorentzian_at_nodes = numf ./ ((ze2_row - zaf).^2 + B1f);    % [nf x n2]
                plain_lorentzian_est = termboost * lorentzian_at_nodes .* dzvec2_row;
                contrib(flaggedrows,:) = contrib(flaggedrows,:) + (wt_physical - plain_lorentzian_est);
            end
            patched = patched + contrib;
        end
        % Back to [n1 x n2 x nr], then to [n1n2 x nr]
        patched = permute(reshape(patched,n1,nr,n2),[1 3 2]);
        Fsub = -ny/4/2/pi*thinplaneboost * (dzvec1(:) .* patched) .* reshape(kernelfactor,n1,n2,nr);
        Fsub = reshape(Fsub,n1*n2,nr);

        ncorrected = sum(anyflag,1);       % [1 x nr]
        propdiagnostics.npairs_eligible = propdiagnostics.npairs_eligible + nr;
        propdiagnostics.nrows           = propdiagnostics.nrows + n1*nr;
        propdiagnostics.nrows_corrected = propdiagnostics.nrows_corrected + sum(ncorrected);
        propdiagnostics.npairs_corrected = propdiagnostics.npairs_corrected + nnz(ncorrected);
        if any(ncorrected)
            ncorrectedpairs = ncorrectedpairs + 1;
            correctedpairs(ncorrectedpairs,:) = [edge1 edge2 sum(ncorrected) n1*nr];
        end
    else
        beta =  sin(ny*(pi + thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) + thetaout   ))) + ...
                sin(ny*(pi + thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) - thetaout   ))) + ...
                sin(ny*(pi - thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) + thetaout   ))) + ...
                sin(ny*(pi - thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) - thetaout   )));
        Fsub = -ny/4/2/pi*beta.*kernelfactor.*(dzvec2(n2vertmat).*dzvec1(n1vertmat))*thinplaneboost;
    end

    propdiagnostics.npairs = propdiagnostics.npairs + nr;
    Fgeom(rr,Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) = Fsub.';
end
propdiagnostics.correctedpairs = correctedpairs(1:ncorrectedpairs,:);

end

%--------------------------------------------------------------------------
% Copies of the local functions in EDcalcpropagatematrix
function [x,w] = gaussLegendreNodes(n)
    % Golub-Welsch algorithm: standard n-point Gauss-Legendre nodes and
    % weights on [-1,1]. Verified to match EDdistelements(n,2)'s hardcoded
    % tables (mapped to [0,1]) to ~1e-16.
    if n == 1
        x = 0; w = 2; return
    end
    ivec = (1:n-1);
    betavec = ivec ./ sqrt(4*ivec.^2 - 1);
    J = diag(betavec,1) + diag(betavec,-1);
    [V,D] = eig(J);
    x = diag(D);
    [x,idx] = sort(x);
    V = V(:,idx);
    w = 2*(V(1,:).^2)';
    x = x(:); w = w(:);
end

function P = legendreVandermonde(x,Nmax)
    % P(j,n+1) = P_n(x(j)), n=0..Nmax-1
    x = x(:);
    M = length(x);
    P = zeros(M,Nmax);
    P(:,1) = 1;
    if Nmax > 1
        P(:,2) = x;
    end
    for n = 2:Nmax-1
        P(:,n+1) = ((2*n-1)*x.*P(:,n) - (n-1)*P(:,n-1))/n;
    end
end
function [x,w,S,tfine,wfine,Pint] = glStandardCached(n)
    % Standard n-point Gauss-Legendre nodes x and weights w on [-1,1]
    % (ascending), and S(j,m+1) = (2m+1)*P_m(x_j), m = 0..n-1.
    % Also a finer 4n-point rule (tfine,wfine) and the matrix Pint
    % [4n x n] which interpolates nodal values at x to the fine nodes
    % (exact for polynomials of degree <= n-1).
    % These depend only on n, so they are cached between calls.
    persistent cachedn cachedx cachedw cachedS cachedt cachedwf cachedP
    if isempty(cachedn)
        cachedn = [];  cachedx = {};  cachedw = {};  cachedS = {};
        cachedt = {};  cachedwf = {};  cachedP = {};
    end
    idx = find(cachedn == n,1);
    if isempty(idx)
        [x,w] = gaussLegendreNodes(n);
        P = legendreVandermonde(x,n);
        S = P .* (2*(0:n-1)+1);
        [tfine,wfine] = gaussLegendreNodes(4*n);
        Pfine = legendreVandermonde(tfine,n);          % [4n x n]
        Pint = (Pfine*S.') .* (w(:).'/2);              % [4n x n]
        cachedn(end+1) = n;
        cachedx{end+1} = x;   cachedw{end+1} = w;   cachedS{end+1} = S;
        cachedt{end+1} = tfine;  cachedwf{end+1} = wfine;  cachedP{end+1} = Pint;
    else
        x = cachedx{idx};  w = cachedw{idx};  S = cachedS{idx};
        tfine = cachedt{idx};  wfine = cachedwf{idx};  Pint = cachedP{idx};
    end
end

function Q = neumannQupward(z,Nmax)
    % Legendre functions of the second kind Q_0..Q_{Nmax-1} at the complex
    % values z (off the interval [-1,1]), by upward recurrence from the
    % closed forms of Q_0 and Q_1. Returns [nz x Nmax]. Vectorised over z.
    % The upward recurrence is stable only for z close to [-1,1]; the
    % caller uses it only where rho^(2*Nmax) < 1e4, with
    % rho = |z + sqrt(z^2-1)|.
    zcol = z(:);
    Q = zeros(length(zcol),Nmax);
    Q(:,1) = 0.5*log((zcol+1)./(zcol-1));
    if Nmax > 1
        Q(:,2) = zcol.*Q(:,1) - 1;
    end
    for n = 1:Nmax-2
        Q(:,n+2) = ( (2*n+1)*zcol.*Q(:,n+1) - n*Q(:,n) )/(n+1);
    end
end

function WT = modifiedWeightsLorentzianVec(S,wnodes,tfine,wfine,Pint,za,B)
    % Modified weights for g(t) = 1/((t-za)^2 + B) on the standard GL
    % nodes, for a whole vector of (za,B) pairs at once:
    %   WT(j,q) = integral over [-1,1] of g_q(t)*L_j(t) dt
    % where L_j is the Lagrange polynomial of node j.
    %   S,wnodes,tfine,wfine,Pint   from glStandardCached
    %   za,B    [nz x 1], real, B >= 0
    % Returns WT [N x nz], one column per (za,B) pair.
    %
    % Two cases, decided by rho = |p + sqrt(p^2-1)|, p = za + i*sqrt(B):
    %
    % 1. Narrow peak (pole close to the interval, rho^(2N) < 1e4):
    %      1/((t-za)^2+B) = Im( 1/(t-p) ) / sqrt(B)
    %    and the single-pole weights follow from Neumann's formula
    %      int P_m(t)/(t-z) dt = -2 Q_m(z):
    %      w_j(z) = -w_j * sum_m (2m+1) P_m(x_j) Q_m(z)
    %    The conjugate pole gives the complex conjugate, so one
    %    evaluation is enough.
    %
    % 2. Wide peak (rho^(2N) >= 1e4): the 4N-point rule integrates
    %    g(t)*L_j(t) to rounding accuracy (error ~ rho^(-7N) <= 1e-14),
    %    so the weights are obtained directly from one matrix product.
    za = za(:);
    sB = sqrt(max(B(:),eps));
    N = size(S,1);
    nz = length(za);
    WT = zeros(N,nz);
    p = za + 1i*sB;
    rho = abs(p + sqrt(p.^2-1));
    rho = max(rho,1./rho);
    isnarrow = (2*N)*log(rho) < log(1e4);

    iv = find(isnarrow);
    if ~isempty(iv)
        Q = neumannQupward(p(iv),N);                              % [nnarrow x N]
        WT(:,iv) = -(S*imag(Q).') .* wnodes(:) ./ (sB(iv).');
    end

    iv = find(~isnarrow);
    if ~isempty(iv)
        gfine = 1./( (tfine(:) - za(iv).').^2 + (sB(iv).').^2 );    % [4N x nwide]
        WT(:,iv) = Pint.'*(wfine(:).*gfine);
    end
end
