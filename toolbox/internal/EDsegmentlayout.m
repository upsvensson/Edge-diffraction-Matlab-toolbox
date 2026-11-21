function [segdata,Sdata_regular] = EDsegmentlayout(edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,planedata,controlparameters)
% EDsegmentlayout - Edge segments for sources that see only part of an edge (non-convex models).
% When a source point (a source, or an image source from EDimagepointdata)
% sees only a part [c1,c2] of an edge, a source term on the full-edge
% nodes would have to be zero outside [c1,c2]. That step inside the
% integration interval limits the Gauss-Legendre quadrature (in H*q and in
% the propagation, which integrate over the from-edge) to an erratic
% error of order 1/n. Instead, the part [c1,c2] gets its own Gauss-Legendre
% rule, as a "segment": a q segment (X <- segment) for every edge X that
% the edge sees. These q segments are known (the source term only), and
% are not unknowns of the integral equation. They enter in two ways:
%   1. One step of the H-matrix takes them to the ordinary unknowns,
%      b(Y <- X) = int_segment H(Y <- X <- z1) q0seg(X <- z1) dz1,
%      which is added to the right-hand side of the integral equation.
%   2. They are propagated to the receivers directly (the second-order
%      diffraction source -> segment -> X -> receiver).
% Every integral over the segment has a smooth integrand, so the
% quadrature converges regularly.
%
% This function finds the segments (groups), removes them from the source
% visibility used for the ordinary source term, and builds the layouts of
% the segment q vectors, in the same form as Hsubmatrixdata (the fields
% that EDinteg_souterm, EDcalcpropagatematrix_block and EDpropagateQvec
% use), with the extra fields .fromnodelist (the nodes and weights of the
% segment), .segmentlimits ([c1 c2]) and .regpairindex (the ordinary q
% segment with the same edge pair).
%
% Input parameters:
%   edgedata, edgetoedgedata, Hsubmatrixdata, Sdata, planedata
%                   structs, see help EDdatastructs. Sdata is the one used
%                   by the integral equation (with image sources, from
%                   EDimagepointdata).
%   controlparameters   struct; not used now (for future options)
%
% Output parameters:
%   segdata         struct with the fields
%     .ngroups      the number of segments
%     .groups       [ngroups,4] = [source point, edge, c1, c2]
%     .grouplayout  {ngroups,1} the layout of each segment's q vector
%     .groupoffset  [ngroups,1] the position of each segment's q vector in
%                   the combined vector (0-based)
%     .layout       the layout of the combined vector (all segments)
%   Sdata_regular   Sdata, where the segments' edges are marked as not
%                   visible for their source points (for the ordinary
%                   source term)
%
% Number of nodes: a segment gets as many nodes as its full edge
% (Hsubmatrixdata.nedgeelems), so that a convergence study with ngauss
% refines the segments too.
%
% See also EDintegralequation_convex_tf, EDinteg_souterm,
% EDcalcpropagatematrix_block, EDimagepointdata, EDedgepairvisibility
%
% Peter Svensson 21 Nov. 2026
%
% [segdata,Sdata_regular] = EDsegmentlayout(edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,planedata,controlparameters);

% 21 Nov. 2026 First version (by Claude, prompted by PS).

tol = 1e-9;
nedges = size(edgedata.edgecorners,1);
nsourcepoints = size(Sdata.coordinates,1);
Sdata_regular = Sdata;

segdata = struct('ngroups',0,'groups',zeros(0,4),'grouplayout',{{}},'groupoffset',zeros(0,1),'layout',[]);

if ~isfield(Sdata,'vispartedgesfroms_start') || ~isfield(Sdata,'vispartedgesfroms_end')
    return
end
visstart = double(Sdata.vispartedgesfroms_start);
visend = double(Sdata.vispartedgesfroms_end);
visflag = double(Sdata.vispartedgesfroms) ~= 0;
partial = visflag & (visstart > tol | visend < 1-tol);
if ~any(partial(:))
    return
end

[iedge,ipoint] = find(partial);
groups = [ipoint iedge visstart(partial) visend(partial)];
groups = sortrows(groups,[1 2 3]);
ngroups = size(groups,1);

% The ordinary source term: these edges are not visible for these points
vp = double(Sdata_regular.vispartedgesfroms);
vp(sub2ind(size(vp),groups(:,2),groups(:,1))) = 0;
Sdata_regular.vispartedgesfroms = vp;

grouplayout = cell(ngroups,1);
groupoffset = zeros(ngroups,1);
offset = 0;
allpairs = zeros(0,2);
allnodes = cell(0,1);
alllimits = zeros(0,2);
allreg = zeros(0,1);
allthin = zeros(0,1);
allvis = zeros(0,1);
anyvis = false;

for g = 1:ngroups
    e = groups(g,2);
    c1 = groups(g,3);
    c2 = groups(g,4);
    na = double(Hsubmatrixdata.nedgeelems(e));
    [x01,w01] = EDgetGLrule(na);
    nodes = [c1 + (c2-c1)*x01(:), (c2-c1)*w01(:)];

    reg = find(Hsubmatrixdata.edgepairlist(:,2) == e);
    npairs = numel(reg);
    L = struct;
    L.edgepairlist = Hsubmatrixdata.edgepairlist(reg,:);
    L.nedgeelems = Hsubmatrixdata.nedgeelems;
    nto = double(Hsubmatrixdata.nedgeelems(L.edgepairlist(:,1)));
    L.Qvecendnums = cumsum(na*nto(:));
    L.Qvecstartnums = [1;L.Qvecendnums(1:end-1)+1];
    L.isthinplaneedgepair = Hsubmatrixdata.isthinplaneedgepair(reg);
    L.quadraturematrix_pos = Hsubmatrixdata.quadraturematrix_pos;
    L.quadraturematrix_weights = Hsubmatrixdata.quadraturematrix_weights;
    L.fromnodelist = repmat({nodes},npairs,1);
    L.segmentlimits = repmat([c1 c2],npairs,1);
    L.regpairindex = reg(:);

    % Partly visible edge pairs (EDedgepairvisibility): the node pairs of
    % the segment and the to-edge
    L.visrowmask = [];
    if isfield(edgetoedgedata,'partialvispairs') && ~isempty(edgetoedgedata.partialvispairs)
        pv = sort(double(edgetoedgedata.partialvispairs),2);
        vis = ones(L.Qvecendnums(end),1);
        for p = 1:npairs
            X = double(L.edgepairlist(p,1));
            [ispv,locpv] = ismember(sort([X e]),pv,'rows');
            if ispv
                nX = double(Hsubmatrixdata.nedgeelems(X));
                uX = full(Hsubmatrixdata.quadraturematrix_pos(1:nX,nX));
                v = EDedgepairvisibility('nodepairs',planedata,edgedata,e,X,nodes(:,1),uX,...
                    edgetoedgedata.partialvisoccluders{locpv});
                vis(L.Qvecstartnums(p):L.Qvecendnums(p)) = v(:);
            end
        end
        if any(vis ~= 1)
            L.visrowmask = vis;
            anyvis = true;
        end
    end

    grouplayout{g} = L;
    groupoffset(g) = offset;
    nQg = L.Qvecendnums(end);
    offset = offset + nQg;

    allpairs = [allpairs;L.edgepairlist]; %#ok<AGROW>
    allnodes = [allnodes;L.fromnodelist]; %#ok<AGROW>
    alllimits = [alllimits;L.segmentlimits]; %#ok<AGROW>
    allreg = [allreg;L.regpairindex]; %#ok<AGROW>
    allthin = [allthin;double(L.isthinplaneedgepair(:))]; %#ok<AGROW>
    if isempty(L.visrowmask)
        allvis = [allvis;ones(nQg,1)]; %#ok<AGROW>
    else
        allvis = [allvis;L.visrowmask]; %#ok<AGROW>
    end
end

% The combined layout (all segments), for the propagation
C = struct;
C.edgepairlist = allpairs;
C.nedgeelems = Hsubmatrixdata.nedgeelems;
nfrom = cellfun(@(x) size(x,1),allnodes);
nto = double(Hsubmatrixdata.nedgeelems(allpairs(:,1)));
C.Qvecendnums = cumsum(nfrom(:).*nto(:));
C.Qvecstartnums = [1;C.Qvecendnums(1:end-1)+1];
C.isthinplaneedgepair = allthin;
C.quadraturematrix_pos = Hsubmatrixdata.quadraturematrix_pos;
C.quadraturematrix_weights = Hsubmatrixdata.quadraturematrix_weights;
C.fromnodelist = allnodes;
C.segmentlimits = alllimits;
C.regpairindex = allreg;
if anyvis
    C.visrowmask = allvis;
else
    C.visrowmask = [];
end

segdata.ngroups = ngroups;
segdata.groups = groups;
segdata.grouplayout = grouplayout;
segdata.groupoffset = groupoffset;
segdata.layout = C;

if nedges == 0 || nsourcepoints == 0
    segdata.ngroups = 0;
end

end
