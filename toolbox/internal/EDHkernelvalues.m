function [Hvalues,mdist] = EDHkernelvalues(edgedata,edgetoedgedata,edge1,edge2,edge3,...
    n1vec,n2vec,n3vec,dzvec,n1ind,n2ind,n3ind,thinplaneboost)
% EDHkernelvalues - The frequency-independent part of the H-matrix of the edge source integral equation.
% For the edge triplet edge1 -> edge2 -> edge3 (edge 2 is the diffracting
% edge), and the node positions n1vec, n2vec, n3vec on the three edges,
% the kernel values
%   (2-acrossface_out)*Hsubconstant*ny*beta/m*dz1*thinplaneboost
% are computed for a list of node combinations (n1ind(i), n2ind(i),
% n3ind(i)). The frequency-dependent factor exp(-1i*k*m) is not included;
% m is returned, so that the caller can apply it.
%
% Used by EDintegralequation_convex_tf (PART 1, with the index lists of
% the H-submatrix layout) and EDrefineQvec (with full grids of new nodes).
%
% Input parameters:
%   edgedata, edgetoedgedata   structs, see help EDdatastructs
%   edge1, edge2, edge3        the edge triplet (edgetripletlist(nn,:) is
%                       [edge3 edge2 edge1])
%   n1vec, n2vec, n3vec        the normalized node positions (0 to 1) on
%                       the three edges
%   dzvec               [n1,1] the quadrature weights on edge 1, in m
%   n1ind, n2ind, n3ind        index vectors into n1vec, n2vec, n3vec, one
%                       element per kernel value. They must be column
%                       vectors (or arrays of the same size).
%   thinplaneboost      1, or 2 for a thin-plane triplet
%
% Output parameters:
%   Hvalues             the kernel values, one per element of the index
%                       vectors
%   mdist               the distance m from the node on edge 1 to the
%                       node on edge 2, one per element of the index vectors
%
% See also EDintegralequation_convex_tf, EDrefineQvec
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Oct. 2026
%
% [Hvalues,mdist] = EDHkernelvalues(edgedata,edgetoedgedata,edge1,edge2,edge3,...
%    n1vec,n2vec,n3vec,dzvec,n1ind,n2ind,n3ind,thinplaneboost);

% 28 Oct. 2026 First version (by Claude, prompted by PS): the kernel
% computation of PART 1 of EDintegralequation_convex_tf, moved here
% unchanged, so that EDrefineQvec can use the same code.

Hsubconstant = -1/4/pi/2;
refto = edgetoedgedata.reftoshortlistE;
len2 = edgedata.edgelengthvec(edge2);

n2vec = n2vec(:);
n3vec = n3vec(:);

% Bug found and fixed 22 Sep. 2014: lines below did not always detect if another
% edge was skewed relative to the first one, because only the
% thetae1sho was used. Necessary to check thetae1sho and thetae2sho.

ny = pi/(2*pi-edgedata.closwedangvec(edge2));
thetain  = edgetoedgedata.thetae1sho(refto(edge1,edge2));
thetain_other_end_of_edge  = edgetoedgedata.thetae2sho(refto(edge1,edge2));
if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_other_end_of_edge) < 1e-10
   acrossface_in = 1;
else
   acrossface_in = 0;
end
thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
thetaout_other_end_of_edge = edgetoedgedata.thetae2sho(refto(edge3,edge2));
if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetaout-thetaout_other_end_of_edge) < 1e-10
   acrossface_out = 1;
else
   acrossface_out = 0;
end

ze2_re2 = n2vec*len2;
ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2))   );

n1vec = n1vec(:).';

ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );

if acrossface_in*acrossface_out == 0 % This case is a preparation for non-convex objects.
    % Reconvert the cylindrical coordinates of
    % edge1-re2 into cartesian! Then we can use the EDcoordtrans.
    % We define our own cartesian coord syst such that the reference
    % edge (edge2) has its starting point in [0 0 0], the z-axis along its
    % own edge, and the x-axis is along (onto) the reference plane.

    xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
    xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
    yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
    yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
    xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
    ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);

    [re1_re2,thetae1_re2,~] = EDcoordtrans1([xe1_re2.' ye1_re2.' ze1_re2.'],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
    re1_re2 = re1_re2(:).';
    thetae1_re2 = thetae1_re2(:).';

    xedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
    xedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
    yedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
    yedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
    xe3_re2 = xedge3start + n3vec.*(xedge3end-xedge3start);
    ye3_re2 = yedge3start + n3vec.*(yedge3end-yedge3start);

    [re3_re2,thetae3_re2,~] = EDcoordtrans1([xe3_re2 ye3_re2 ze3_re2],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
    re3_re2 = re3_re2(:);
    thetae3_re2 = thetae3_re2(:);

else

    re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );
    re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

end

ldist = (ze3_re2(n3ind) - ze2_re2(n2ind)).^2   + re3_re2(n3ind).^2;
ldist = sqrt(ldist);

mdist = (ze1_re2(n1ind).' - ze2_re2(n2ind)).^2   + (re1_re2(n1ind).').^2;
mdist = sqrt(mdist);

% The factor exp(-1i*k*mdist) is not included here: it is
% applied per frequency (see EDintegralequation_convex_tf).
dz_overm = dzvec(n1ind)./mdist;

ch = ((ze2_re2(n2ind) - ze1_re2(n1ind).').*(ze2_re2(n2ind) - ze3_re2(n3ind)) + mdist.*ldist)./re1_re2(n1ind).'./re3_re2(n3ind);
ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
ch = ( ch + 1./ch)/2;

% PS 28 March 2011: The 4 beta terms are identical!
if acrossface_in == 1 && acrossface_out == 1
   beta = 4*( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))));
else
   beta = ( sin(ny*(pi + thetae1_re2(n1ind).' + thetae3_re2(n3ind)  ))./(ch - cos(ny*(pi + thetae1_re2(n1ind).' + thetae3_re2(n3ind)   ))+eps*10)) + ...
          ( sin(ny*(pi + thetae1_re2(n1ind).' - thetae3_re2(n3ind)  ))./(ch - cos(ny*(pi + thetae1_re2(n1ind).' - thetae3_re2(n3ind)   ))+eps*10)) + ...
          ( sin(ny*(pi - thetae1_re2(n1ind).' + thetae3_re2(n3ind)  ))./(ch - cos(ny*(pi - thetae1_re2(n1ind).' + thetae3_re2(n3ind)   ))+eps*10)) + ...
          ( sin(ny*(pi - thetae1_re2(n1ind).' - thetae3_re2(n3ind)  ))./(ch - cos(ny*(pi - thetae1_re2(n1ind).' - thetae3_re2(n3ind)   ))+eps*10));
end

if ismember(edge1,edgedata.offedges)==1 || ismember(edge2,edgedata.offedges)==1 || ismember(edge3,edgedata.offedges)==1
   beta = beta*0;
end

Hvalues = (2-acrossface_out)*Hsubconstant*ny*beta.*dz_overm*thinplaneboost;

end
