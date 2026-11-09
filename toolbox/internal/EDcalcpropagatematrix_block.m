function [Fgeom,pathlengths,propdiagnostics] = EDcalcpropagatematrix_block(edgedata,edgetoedgedata,...
    Hsubmatrixdata,vispartedgesfromr,rRmat,thetaRmat,zRmat,doesQsegmenthavevalues,polemethod)
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
% NEAR-SINGULAR ROWS
% For an edge pair across a face, each of the two singular beta terms,
%   beta(z2) = sin(ny*phi)/(cosh(ny*eta(z2)) - cos(ny*phi)),
% has a pair of complex-conjugate poles close to the apex za when the
% receiver is close to a zone boundary (the plane of a face). A row (an
% edge-1 node and a receiver) is flagged when the peak is narrow compared
% with the node spacing. Two treatments of the flagged rows (polemethod):
%   'taylor'  (ESIEpropagationimplementation = 'accurate') The poles are
%             taken from the Taylor expansion of Svensson & Calamia (2006),
%             za +- 1i*sqrt(B1). The Lorentzian model is integrated by
%             product integration and the remainder beta - model with the
%             plain Gauss-Legendre rule. The remainder is not smooth on the
%             scale of the peak (the true poles are shifted, and its odd
%             part is nearly singular), which gives an irregular, slowly
%             converging error.
%   'exact'   (ESIEpropagationimplementation = 'exactpole') The exact pole
%             p = zp + 1i*sqrt(Bp) is found by Newton iteration, started from
%             the Taylor values, with a cancellation-free form of
%               F(z) = cosh(ny*eta) - cos(ny*phi) = 2*sinh(ny*g)^2 + 2*sin(ny*phi/2)^2,
%               g = asinh(sqrt(u/2)),
%               u = ch0 - 1 = (rS+rR)^2 (z-za)^2 / (rS rR (m l + rS rR - (z-zS)(z-zR))).
%             Then beta(z)*((z-zp)^2+Bp) is smooth, and the whole row is
%             integrated by product integration with the weight function
%             1/((z-zp)^2+Bp): row element j = termboost*W_j*beta_j*((z_j-zp)^2+Bp),
%             where the smooth factor 1/(l*m)*q is represented by the
%             polynomial through the nodes (singularity swap quadrature,
%             af Klinteberg & Barnett, BIT Numer. Math. 2021). For receivers
%             very close to the zone boundary (sin(ny*phi/2)^2 < 1e-6), beta
%             at the nodes is recomputed with the cancellation-free form. If
%             the Newton iteration does not converge for a row, the 'taylor'
%             treatment is used for that row.
%
% For each receiver r and frequency (wavenumber k), the propagation row
% vector is Fgeom(r,:).*exp(-1i*k*(pathlengths(r,:)-Rstart)).
%
% EDGE NUMBERING (the same in all the integral equation functions)
% Each kernel is one first-order diffraction step, from an incoming point
% via a diffracting edge to an outgoing point, and the diffracting edge is
% always edge 2 (position z2):
%   source term:      q0(z3 <- z2) = G(z3 <- z2 <- S)
%   integral eq.:     q(z3 <- z2)  = q0(z3 <- z2) + int H(z3 <- z2 <- z1) q(z2 <- z1) dz1
%   propagation:      p(R)         = int int F(R <- z2 <- z1) q(z2 <- z1) dz1 dz2
% q(zb <- za) is the edge source amplitude at za (on the edge that the q
% segment is "from"), radiating towards zb (on the "to" edge). The wedge
% index ny, the beta terms and the cylindrical coordinates (the suffix
% _re2, e.g. ze1_re2) are those of edge 2. S or R takes the place of
% edge 1 or edge 3 where the kernel starts or ends at a point.
% In Hsubmatrixdata, edgepairlist(ii,:) = [to-edge from-edge] of the q
% segment ii, and edgetripletlist(nn,:) = [edge3 edge2 edge1]. Inside a q
% segment, the node index of the from-edge runs fastest: element
% ifrom + (ito-1)*nfrom.
% In this function: edge1 -> edge2 -> R. The q segment ii has
% edgepairlist(ii,:) = [edge2 edge1].
%
% Input parameters:
%   edgedata, edgetoedgedata, Hsubmatrixdata   Structs (see EDdatastructs)
%   vispartedgesfromr   [nedges,nrec] the visibility of the edges from
%                       the receivers of the block (Rdata.vispartedgesfromr)
%   rRmat, thetaRmat, zRmat   [nedges,nrec] the cylindrical coordinates of
%                       the receivers of the block relative to each edge
%   doesQsegmenthavevalues    [nedgepairs,1] 1 for the edge pairs where the
%                       edge source amplitudes are not all zero
%   polemethod          (optional) 'taylor' (default) or 'exact', the
%                       treatment of the near-singular rows (see above)
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
%                   pair with corrections. With polemethod = 'exact' also
%                   .nrows_exactpole and .nrows_fallback: the flagged rows
%                   (per beta term) with the exact pole, and those where the
%                   Newton iteration failed ('taylor' treatment).
%
% See also EDcalcpropagatematrix, EDintegralequation_convex_tf
%
% Peter Svensson (peter.svensson@ntnu.no) 9 Nov. 2026
%
% [Fgeom,pathlengths,propdiagnostics] = EDcalcpropagatematrix_block(edgedata,...
%    edgetoedgedata,Hsubmatrixdata,vispartedgesfromr,rRmat,thetaRmat,zRmat,...
%    doesQsegmenthavevalues,polemethod);

% 16 Oct. 2026 First version (by Claude, prompted by PS). The local
% functions gaussLegendreNodes, legendreVandermonde, glStandardCached,
% neumannQupward and modifiedWeightsLorentzianVec are copies of those in
% EDcalcpropagatematrix.
% 18 Oct. 2026 (by Claude, prompted by PS) The edge numbering (edge 2 is
% always the diffracting edge) is described in the header.
% 9 Nov. 2026 (by Claude, prompted by PS) New optional input polemethod:
% 'exact' finds the exact poles of the near-singular beta terms by Newton
% iteration and integrates the flagged rows by product integration with
% the exact pole pair as the weight function. 'taylor' (the default) is
% the earlier treatment, unchanged. New propdiagnostics fields
% .nrows_exactpole and .nrows_fallback (only for 'exact').

if nargin < 9 || isempty(polemethod)
    polemethod = 'taylor';
end
if ~ischar(polemethod) || ~any(strcmp(polemethod,{'taylor','exact'}))
    error('EDcalcpropagatematrix_block:polemethod','polemethod must be ''taylor'' or ''exact''')
end
useexactpole = strcmp(polemethod,'exact');

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
nrows_exactpole = 0;
nrows_fallback = 0;

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

            if useexactpole
                flaggedrows = find(flag_k(:));
                if ~isempty(flaggedrows)
                    [~,wnodes_std,Sstd,tfine_std,wfine_std,Pint_std] = glStandardCached(n2);
                    zaf  = zavec(flaggedrows);
                    B1f  = max(B1_k(flaggedrows),0);
                    numf = num_k(flaggedrows);
                    % The geometry of each flagged row (i,r)
                    [iif,rrf] = ind2sub([n1 nr],flaggedrows);
                    zSf = ze1_re2(iif);  zSf = zSf(:);
                    rSf = re1_re2(iif);  rSf = rSf(:);
                    zRf = zR(rrf);       zRf = zRf(:);
                    rRf = rR(rrf);       rRf = rRf(:);
                    phif = phi_k(rrf);   phif = phif(:);
                    sinhalf2 = sin(ny*phif/2).^2;
                    snphi = sin(ny*phif);
                    [zpf,Bpf,okf] = edexactpoles(zSf,rSf,zRf,rRf,zaf(:),B1f(:),sinhalf2,ny,len2);
                    iok = find(okf);
                    inok = find(~okf);
                    nrows_exactpole = nrows_exactpole + numel(iok);
                    nrows_fallback = nrows_fallback + numel(inok);
                    if ~isempty(iok)
                        WT = modifiedWeightsLorentzianVec(Sstd,wnodes_std,tfine_std,wfine_std,Pint_std,...
                            2*zpf(iok)/len2-1,4*Bpf(iok)/len2^2);                       % [n2 x nok]
                        Wg = (2/len2)*WT.';                                               % [nok x n2]
                        % beta at the nodes: the values already computed, except for
                        % receivers extremely close to the zone boundary, where
                        % cosh - cos loses digits (relative error ~ eps/sin^2(ny*phi/2))
                        betaacc = contrib(flaggedrows(iok),:)./(termboost*dzvec2_row);  % [nok x n2]
                        irec = find(sinhalf2(iok) < 1e-6);
                        if ~isempty(irec)
                            ik = iok(irec);
                            betaacc(irec,:) = edbetaaccurate(ze2_row,zSf(ik),rSf(ik),zRf(ik),rRf(ik),...
                                zaf(ik),sinhalf2(ik),snphi(ik),ny);
                        end
                        contrib(flaggedrows(iok),:) = termboost*Wg.*betaacc.*((ze2_row - zpf(iok)).^2 + Bpf(iok));
                    end
                    if ~isempty(inok)
                        % The original correction (Taylor model of the pole)
                        fr = flaggedrows(inok);
                        xa   = 2*zaf(inok)/len2 - 1;
                        Bhat = 4*B1f(inok)/len2^2;
                        WT = modifiedWeightsLorentzianVec(Sstd,wnodes_std,tfine_std,wfine_std,Pint_std,xa,Bhat);
                        wt_physical = termboost*(2/len2) * (numf(inok) .* WT.');
                        lorentzian_at_nodes = numf(inok) ./ ((ze2_row - zaf(inok)).^2 + B1f(inok));
                        plain_lorentzian_est = termboost * lorentzian_at_nodes .* dzvec2_row;
                        contrib(fr,:) = contrib(fr,:) + (wt_physical - plain_lorentzian_est);
                    end
                end
            else
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
if useexactpole
    propdiagnostics.nrows_exactpole = nrows_exactpole;
    propdiagnostics.nrows_fallback = nrows_fallback;
end

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

%--------------------------------------------------------------------------
function [F,dF] = edFexact(z,zS,rS,zR,rR,za,sinhalf2,ny)
    % F(z) = cosh(ny*eta(z)) - cos(ny*phi) without cancellation, and dF/dz.
    % z can be complex. All inputs are column vectors (or broadcast).
    a = z - zS;  b = z - zR;
    m = sqrt(a.^2 + rS.^2);  l = sqrt(b.^2 + rR.^2);
    A = (rS + rR).^2;
    D = rS.*rR.*(m.*l + rS.*rR - a.*b);
    u = A.*(z - za).^2./D;                  % = ch0 - 1
    s = sqrt(u/2);
    g = asinh(s);                           % = eta/2
    F = 2*sinh(ny*g).^2 + 2*sinhalf2;
    if nargout > 1
        Dp = rS.*rR.*(a./m.*l + m.*b./l - (a + b));
        du = A.*(2*(z - za).*D - (z - za).^2.*Dp)./D.^2;
        small = abs(s) < 1e-8;
        ss = s;  ss(small) = 1;
        dFdu = ny*sinh(2*ny*g)./(2*ss.*sqrt(1 + s.^2));
        dFdu(small) = ny^2;
        dF = dFdu.*du;
    end
end

function [zp,Bp,ok] = edexactpoles(zS,rS,zR,rR,za,B1,sinhalf2,ny,len2)
    % The zero p (upper half-plane) of F(z) = cosh(ny*eta) - cos(ny*phi)
    % near the apex za, by Newton iteration from the Taylor estimate
    % za + 1i*sqrt(B1). ok = 0 where the iteration fails (then the caller
    % uses the Taylor model).
    p = za + 1i*sqrt(B1);
    active = true(size(p));
    for it = 1:12
        ia = find(active);
        if isempty(ia)
            break
        end
        [F,dF] = edFexact(p(ia),zS(ia),rS(ia),zR(ia),rR(ia),za(ia),sinhalf2(ia),ny);
        step = F./dF;
        step(~isfinite(step)) = 0;
        p(ia) = p(ia) - step;
        % Rows that have converged (quadratic convergence) drop out
        active(ia) = abs(step) > 1e-14*abs(imag(p(ia)));
    end
    F = edFexact(p,zS,rS,zR,rR,za,sinhalf2,ny);
    ok = isfinite(p) & imag(p) > 0 & sinhalf2 > 1e-30 & ...
         abs(F) <= 1e-10*2*sinhalf2 & abs(real(p) - za) < 0.25*len2 + 4*sqrt(B1);
    zp = real(p);
    Bp = imag(p).^2;
end

function beta = edbetaaccurate(zrow,zS,rS,zR,rR,za,sinhalf2,snphi,ny)
    % One beta term at the real nodes zrow [1 x n2], for the rows given by
    % the column vectors zS ... snphi [nrows x 1]; cancellation-free.
    a = zrow - zS;  b = zrow - zR;
    m = sqrt(a.^2 + rS.^2);  l = sqrt(b.^2 + rR.^2);
    u = (rS + rR).^2.*(zrow - za).^2./(rS.*rR.*(m.*l + rS.*rR - a.*b));
    g = asinh(sqrt(u/2));
    beta = snphi./(2*sinh(ny*g).^2 + 2*sinhalf2);
end
