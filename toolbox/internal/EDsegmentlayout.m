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
% 22 Nov. 2026: to-edge splits (group type 2). Where a zone boundary of
% the first edge e, for a source point, crosses the second edge X (with
% the apex point inside the part of e that the source point sees), the
% source term q0(X <- e) jumps along X. The integral equation needs only
% its values at the nodes of X, but the propagation integrates over X,
% and a Gauss-Legendre rule over a jump converges erratically, as O(1/n).
% A type-2 group has two q segments (X <- e): the source term at nodes of
% X that are split at the jump points (.tonodelist, a Gauss-Legendre rule
% with as many nodes as the edge on each part, .tobreaks the part
% limits), and, with the opposite sign, the source term at the ordinary
% nodes of X. Propagated together with the segments, this replaces the
% propagation of the ordinary source term of the pair by the split one.
% It is not added to the right-hand side of the integral equation (no
% H-step). Found with a folded plate and an oblique plane wave, where the
% reflection boundary of a side edge of one plate crossed the far edge of
% the other plate, and the result scattered by up to 3% with ngauss.
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
%     .ngroups      the number of groups (segments and to-edge splits)
%     .groups       [ngroups,4] = [source point, edge, c1, c2]
%     .grouptype    [ngroups,1] 1 = segment, 2 = to-edge split
%     .groupX       [ngroups,1] the edge X of a to-edge split (0 for a
%                   segment)
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
% Peter Svensson 22 Nov. 2026
%
% [segdata,Sdata_regular] = EDsegmentlayout(edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,planedata,controlparameters);

% 21 Nov. 2026 First version (by Claude, prompted by PS).
% 22 Nov. 2026 To-edge splits (group type 2), the fields .grouptype,
% .groupX, and .tonodelist and .tobreaks of the layouts (by Claude,
% prompted by PS). Version 2.

tol = 1e-9;
nedges = size(edgedata.edgecorners,1);
nsourcepoints = size(Sdata.coordinates,1);
Sdata_regular = Sdata;

segdata = struct('ngroups',0,'groups',zeros(0,4),'grouptype',zeros(0,1),'groupX',zeros(0,1),...
    'grouplayout',{{}},'groupoffset',zeros(0,1),'layout',[]);

if ~isfield(Sdata,'vispartedgesfroms_start') || ~isfield(Sdata,'vispartedgesfroms_end')
    return
end
visstart = double(Sdata.vispartedgesfroms_start);
visend = double(Sdata.vispartedgesfroms_end);
visflag = double(Sdata.vispartedgesfroms) ~= 0;
partial = visflag & (visstart > tol | visend < 1-tol);

[iedge,ipoint] = find(partial);
groups = [ipoint iedge visstart(partial) visend(partial)];
groups = sortrows(groups,[1 2 3]);
ngroups1 = size(groups,1);

% 22 Nov. 2026: the to-edge splits (group type 2). For each source point,
% each edge e that it sees (the part [c1,c2]) and each edge X of the pairs
% (X <- e): the points of X where X is on a zone boundary of e for that
% source point, with the apex point inside [c1,c2]. The source term
% q0(X <- e) jumps there, along X.
tosplits = zeros(0,5);      % [source point, e, c1, c2, X]
tosplitpoints = cell(0,1);
if nedges > 0 && nsourcepoints > 0 && isfield(Sdata,'sourcetype') && strcmp(Sdata.sourcetype,'monopole') ...
        && isfield(Sdata,'rSsho') && isfield(Sdata,'reftoshortlistS')
    for col = 1:nsourcepoints
        for e = find(visflag(:,col)).'
            c1 = visstart(e,col);
            c2 = visend(e,col);
            ir = Sdata.reftoshortlistS(e,col);
            rS = double(Sdata.rSsho(ir));
            thetaS = double(Sdata.thetaSsho(ir));
            zS = double(Sdata.zSsho(ir));
            reg = find(Hsubmatrixdata.edgepairlist(:,2) == e);
            for X = unique(double(Hsubmatrixdata.edgepairlist(reg,1))).'
                tlist = edtosplitpoints(edgedata,e,X,rS,thetaS,zS,c1,c2);
                if ~isempty(tlist)
                    tosplits = [tosplits; col e c1 c2 X]; %#ok<AGROW>
                    tosplitpoints{end+1,1} = tlist; %#ok<AGROW>
                end
            end
        end
    end
end
ngroups2 = size(tosplits,1);
ngroups = ngroups1 + ngroups2;
if ngroups == 0
    return
end

% The ordinary source term: these edges are not visible for these points
vp = double(Sdata_regular.vispartedgesfroms);
if ngroups1 > 0
    vp(sub2ind(size(vp),groups(:,2),groups(:,1))) = 0;
end
Sdata_regular.vispartedgesfroms = vp;
grouptype = [ones(ngroups1,1); 2*ones(ngroups2,1)];
groupX = [zeros(ngroups1,1); tosplits(:,5)];
groups = [groups; tosplits(:,1:4)];

grouplayout = cell(ngroups,1);
groupoffset = zeros(ngroups,1);
offset = 0;
allpairs = zeros(0,2);
allnodes = cell(0,1);
alltonodes = cell(0,1);
alltobreaks = cell(0,1);
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
    if grouptype(g) == 2
        % 22 Nov. 2026: a to-edge split. Two q segments (X <- e): the
        % first with the nodes of X split at the jump points, the second
        % with the ordinary nodes of X (with the opposite sign in
        % edsegmentterms). The from-edge nodes are those of the segment
        % [c1,c2] for a partly visible e, otherwise the ordinary ones.
        X = groupX(g);
        reg = find(Hsubmatrixdata.edgepairlist(:,1) == X & Hsubmatrixdata.edgepairlist(:,2) == e);
        reg = reg(1);
        ispartial = c1 > tol || c2 < 1-tol;
        if ispartial
            fromnodes = nodes;
            unodes = nodes(:,1);
        else
            fromnodes = [];
            unodes = full(Hsubmatrixdata.quadraturematrix_pos(1:na,na));
        end
        nX = double(Hsubmatrixdata.nedgeelems(X));
        bp = [0 tosplitpoints{g-ngroups1}(:).' 1];
        [y01,v01] = EDgetGLrule(nX);
        tonodes = zeros(0,2);
        for ib = 1:numel(bp)-1
            tonodes = [tonodes; bp(ib) + (bp(ib+1)-bp(ib))*y01(:), (bp(ib+1)-bp(ib))*v01(:)]; %#ok<AGROW>
        end
        uX = full(Hsubmatrixdata.quadraturematrix_pos(1:nX,nX));
        L = struct;
        L.edgepairlist = Hsubmatrixdata.edgepairlist([reg reg],:);
        L.nedgeelems = Hsubmatrixdata.nedgeelems;
        L.Qvecendnums = cumsum([na*size(tonodes,1); na*nX]);
        L.Qvecstartnums = [1;L.Qvecendnums(1)+1];
        L.isthinplaneedgepair = Hsubmatrixdata.isthinplaneedgepair([reg reg]);
        L.quadraturematrix_pos = Hsubmatrixdata.quadraturematrix_pos;
        L.quadraturematrix_weights = Hsubmatrixdata.quadraturematrix_weights;
        L.fromnodelist = {fromnodes; fromnodes};
        L.segmentlimits = [c1 c2; c1 c2];
        L.tonodelist = {tonodes; []};
        L.tobreaks = {bp; []};
        L.regpairindex = [reg; reg];
        L.visrowmask = [];
        vis = ones(L.Qvecendnums(end),1);
        if isfield(edgetoedgedata,'partialvispairs') && ~isempty(edgetoedgedata.partialvispairs)
            pv = sort(double(edgetoedgedata.partialvispairs),2);
            [ispv,locpv] = ismember(sort([X e]),pv,'rows');
            if ispv
                v = EDedgepairvisibility('nodepairs',planedata,edgedata,e,X,unodes,tonodes(:,1),...
                    edgetoedgedata.partialvisoccluders{locpv});
                vis(L.Qvecstartnums(1):L.Qvecendnums(1)) = v(:);
                if ispartial
                    v = EDedgepairvisibility('nodepairs',planedata,edgedata,e,X,unodes,uX,...
                        edgetoedgedata.partialvisoccluders{locpv});
                    vis(L.Qvecstartnums(2):L.Qvecendnums(2)) = v(:);
                end
            end
        end
        if ~ispartial && isfield(Hsubmatrixdata,'visrowmask') && ~isempty(Hsubmatrixdata.visrowmask)
            vis(L.Qvecstartnums(2):L.Qvecendnums(2)) = ...
                Hsubmatrixdata.visrowmask(Hsubmatrixdata.Qvecstartnums(reg):Hsubmatrixdata.Qvecendnums(reg));
        end
        if any(vis ~= 1)
            L.visrowmask = vis;
            anyvis = true;
        end
        grouplayout{g} = L;
        groupoffset(g) = offset;
        nQg = L.Qvecendnums(end);
        offset = offset + nQg;
        allpairs = [allpairs;L.edgepairlist]; %#ok<AGROW>
        allnodes = [allnodes;L.fromnodelist]; %#ok<AGROW>
        alltonodes = [alltonodes;L.tonodelist]; %#ok<AGROW>
        alltobreaks = [alltobreaks;L.tobreaks]; %#ok<AGROW>
        alllimits = [alllimits;L.segmentlimits]; %#ok<AGROW>
        allreg = [allreg;L.regpairindex]; %#ok<AGROW>
        allthin = [allthin;double(L.isthinplaneedgepair(:))]; %#ok<AGROW>
        allvis = [allvis;vis]; %#ok<AGROW>
        continue
    end

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
    alltonodes = [alltonodes;cell(npairs,1)]; %#ok<AGROW>
    alltobreaks = [alltobreaks;cell(npairs,1)]; %#ok<AGROW>
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
iregfrom = find(nfrom == 0);
nfrom(iregfrom) = double(Hsubmatrixdata.nedgeelems(allpairs(iregfrom,2)));
nto = double(Hsubmatrixdata.nedgeelems(allpairs(:,1)));
ntosp = cellfun(@(x) size(x,1),alltonodes);
nto(ntosp > 0) = ntosp(ntosp > 0);
C.Qvecendnums = cumsum(nfrom(:).*nto(:));
C.Qvecstartnums = [1;C.Qvecendnums(1:end-1)+1];
C.isthinplaneedgepair = allthin;
C.quadraturematrix_pos = Hsubmatrixdata.quadraturematrix_pos;
C.quadraturematrix_weights = Hsubmatrixdata.quadraturematrix_weights;
C.fromnodelist = allnodes;
C.tonodelist = alltonodes;
C.tobreaks = alltobreaks;
C.segmentlimits = alllimits;
C.regpairindex = allreg;
if anyvis
    C.visrowmask = allvis;
else
    C.visrowmask = [];
end

segdata.ngroups = ngroups;
segdata.groups = groups;
segdata.grouptype = grouptype;
segdata.groupX = groupX;
segdata.grouplayout = grouplayout;
segdata.groupoffset = groupoffset;
segdata.layout = C;

if nedges == 0 || nsourcepoints == 0
    segdata.ngroups = 0;
end

end

%--------------------------------------------------------------------------
function tlist = edtosplitpoints(edgedata,e,X,rS,thetaS,zS,c1,c2)
% 22 Nov. 2026: the points t (0..1) along edge X where one of the four
% beta terms of edge e, for the source point (rS,thetaS,zS) (cylindrical
% coordinates relative to e), has its pole (sin(ny*phi/2) = 0, phi = pi
% +- thetaS +- thetaX), with the apex point inside the part [c1,c2] of e.
% Found on 201 sample points, then refined by the secant method (as
% edsplitpoints in EDcalcpropagatematrix_block).
tlist = [];
lene = edgedata.edgelengthvec(e);
ny = pi/(2*pi-edgedata.closwedangvec(e));
nsamp = 201;
tf = linspace(0,1,nsamp).';
[rf,thf] = edXgeom(edgedata,e,X,tf);
% An edge X in a plane through e (constant angle) never crosses a zone
% boundary of e
ok = rf > 1e-9;
if ~any(ok)
    return
end
th0 = thf(find(ok,1));
dth = abs(thf(ok) - th0);
dth = min(dth,abs(dth-2*pi));
if max(dth) < 1e-6
    return
end
s1 = [1 1 -1 -1];
s2 = [1 -1 1 -1];
for kk = 1:4
    gk = sin(ny*(pi + s1(kk)*thetaS + s2(kk)*thf)/2);
    ic = find(gk(1:end-1).*gk(2:end) < 0 & rf(1:end-1) > 1e-9 & rf(2:end) > 1e-9);
    for ic1 = ic(:).'
        ta = tf(ic1); tb = tf(ic1+1);
        ga = gk(ic1); gb = gk(ic1+1);
        for it = 1:8
            tn = tb - gb*(tb-ta)/(gb-ga);
            tn = min(max(tn,tf(ic1)),tf(ic1+1));
            [~,thn] = edXgeom(edgedata,e,X,tn);
            gn = sin(ny*(pi + s1(kk)*thetaS + s2(kk)*thn)/2);
            ta = tb; ga = gb; tb = tn; gb = gn;
            if abs(gn) < 1e-14 || ga == gb
                break
            end
        end
        [rn,~,zn] = edXgeom(edgedata,e,X,tb);
        za = (zn*rS + zS*rn)/(rS + rn);
        if abs(gb) < 1e-8 && za > c1*lene - 1e-6*lene && za < c2*lene + 1e-6*lene && tb > 1e-6 && tb < 1-1e-6
            tlist = [tlist tb]; %#ok<AGROW>
        end
    end
end
if ~isempty(tlist)
    tlist = sort(tlist);
    tlist = tlist([true diff(tlist) > 1e-9]);
end
end

%--------------------------------------------------------------------------
function [r,th,z] = edXgeom(edgedata,e,X,t)
% 22 Nov. 2026: the cylindrical coordinates, relative to edge e, of the
% points t (0..1) along edge X.
t = t(:);
P = edgedata.edgestartcoords(X,:) + t*(edgedata.edgeendcoords(X,:) - edgedata.edgestartcoords(X,:));
[r,th,z] = EDcoordtrans1(P,[edgedata.edgestartcoords(e,:);edgedata.edgeendcoords(e,:)],...
    edgedata.edgenvecs(e,:),reshape(edgedata.edgerelatedcoordsysmatrices(e,:),3,3));
r = r(:); th = th(:); z = z(:);
if edgedata.closwedangvec(e) > 0
    th(th > 2*pi-1e-6) = th(th > 2*pi-1e-6) - 2*pi;
end
end
