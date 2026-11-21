function [out1,out2] = EDedgepairvisibility(mode,planedata,edgedata,varargin)
% EDedgepairvisibility - Visibility between pairs of edges, for non-convex objects.
% Two modes:
%
% 'classify': a conservative classification of edge pairs as fully
% visible, hidden or partly visible from each other.
%   [visclass,occluders] = EDedgepairvisibility('classify',planedata,edgedata,pairs)
% All the straight lines between edge A and edge B lie inside the
% tetrahedron spanned by the four end points. A plane can only obstruct
% some of these lines if it has end points strictly on both sides. If A
% is on one side and B on the other, the cross-section of the tetrahedron
% with the plane is a quadrilateral (the crossing points of the lines
% A1-B1, A1-B2, A2-B2, A2-B1). The plane hides the pair completely if the
% quadrilateral is inside the plane's polygon, and does not obstruct any
% line if the two do not overlap. In all other cases (overlap, or an edge
% that crosses the plane, or an end point in the plane), the plane is
% kept as a possible occluder and the pair is "partly visible".
% The classification is exact when a single plane hides a pair; when
% several planes together hide a pair, it is classified as partly visible,
% and the node pairs are then tested one by one (mode 'nodepairs'). So no
% visible part is ever missed.
% Edge pairs that belong to the same plane (the path runs in that plane)
% get the earlier test: the line between the two edge midpoints. If it
% passes through a plane, or through an edge or corner of a plane (across
% an indent), the pair is hidden, otherwise fully visible. But if one of
% the edges crosses the line of the other (21 Nov. 2026), it is partly on
% the face side and partly across an indent: the pair is partly visible,
% with all the planes (but the edges' own) as possible occluders.
%
% 'nodepairs': the visibility of all the node pairs of two edges.
%   vis = EDedgepairvisibility('nodepairs',planedata,edgedata,edgeA,edgeB,uA,uB,occluderplanes)
% For each pair of points (uA(i) on edge A, uB(j) on edge B), vis(i,j) is
% 1 if the straight line between them is not obstructed by any of the
% planes in occluderplanes, 0 if it is, and 1/2 (2/3) if it passes
% exactly through an edge (a corner) of a plane (the same rule as in
% EDfindGApaths).
%
% Input parameters:
%   mode            'classify' or 'nodepairs'
%   planedata, edgedata   structs, see EDdatastructs
%   pairs           [npairs,2] edge numbers (mode 'classify')
%   edgeA, edgeB    edge numbers (mode 'nodepairs')
%   uA, uB          the node positions on the two edges, as fractions
%                   (0-1) of the edge lengths, from the start point
%   occluderplanes  the planes to test (from mode 'classify')
%
% Output parameters:
%   visclass        [npairs,1] 0: hidden, 1: partly visible, 2: fully
%                   visible (mode 'classify')
%   occluders       {npairs,1} the possible occluders (plane numbers) of
%                   the partly visible pairs; empty for the others
%   vis             [numel(uA),numel(uB)] (mode 'nodepairs')
%
% See also EDed2geo, EDinteg_submatrixstructure, EDfindGApaths, EDpoinpla
%
% Peter Svensson 21 Nov. 2026
%
% [visclass,occluders] = EDedgepairvisibility('classify',planedata,edgedata,pairs);
% vis = EDedgepairvisibility('nodepairs',planedata,edgedata,edgeA,edgeB,uA,uB,occluderplanes);

% 21 Nov. 2026 (by Claude, prompted by PS) Edge pairs in the same plane,
% where one edge crosses the line of the other, are partly visible.
% 20 Nov. 2026 First version (by Claude, prompted by PS). The idea of a
% conservative classification before the point-to-point tests is from
% Antani et al., Applied Acoustics 73 (2012) 218-233.

geo = edgeomdata(planedata,edgedata);

switch mode
    case 'classify'
        pairs = double(varargin{1});
        npairs = size(pairs,1);
        visclass = 2*ones(npairs,1);
        occluders = cell(npairs,1);
        for ii = 1:npairs
            [visclass(ii),occluders{ii}] = edclassifypair(pairs(ii,1),pairs(ii,2),geo);
        end
        out1 = visclass;
        out2 = occluders;
    case 'nodepairs'
        edgeA = double(varargin{1});
        edgeB = double(varargin{2});
        uA = varargin{3}(:);
        uB = varargin{4}(:);
        occluderplanes = varargin{5};
        PA = geo.edgestartcoords(edgeA,:) + uA*(geo.edgeendcoords(edgeA,:) - geo.edgestartcoords(edgeA,:));
        PB = geo.edgestartcoords(edgeB,:) + uB*(geo.edgeendcoords(edgeB,:) - geo.edgestartcoords(edgeB,:));
        nA = numel(uA);
        nB = numel(uB);
        % All pairs, with the index of edge A running fastest
        iA = repmat((1:nA).',nB,1);
        iB = reshape(repmat(1:nB,nA,1),nA*nB,1);
        fac = edobstructionfactor(PA(iA,:),PB(iB,:),geo,occluderplanes(:).');
        out1 = reshape(fac,nA,nB);
        out2 = [];
    otherwise
        error('EDedgepairvisibility:mode','mode must be ''classify'' or ''nodepairs''')
end

end

%------------------------------------------------------------------------
function geo = edgeomdata(planedata,edgedata)
geo = struct;
geo.planeeqs = planedata.planeeqs;
geo.nvecs = planedata.planeeqs(:,1:3);
geo.minvals = planedata.minvals;
geo.maxvals = planedata.maxvals;
geo.planecorners = planedata.planecorners;
geo.corners = planedata.corners;
geo.ncornersperplanevec = planedata.ncornersperplanevec;
geo.rearsideplane = planedata.rearsideplane(:);
geo.obstructplanes = find(planedata.canplaneobstruct(:) ~= 0);
% A thin plate has two planes with the same polygon (each other's
% rearsideplane): only one of them is tested, so that a path exactly
% through an edge of the plate gets the factor 1/2 once (20 Nov. 2026)
rs = planedata.rearsideplane(geo.obstructplanes);
rs = rs(:);
geo.obstructplanes = geo.obstructplanes(~(rs > 0 & ismember(rs,geo.obstructplanes) & rs < geo.obstructplanes));
geo.planesatedge = double(edgedata.planesatedge);
geo.edgestartcoords = edgedata.edgestartcoords;
geo.edgeendcoords = edgedata.edgeendcoords;
% Points closer to a plane than geomacc count as in the plane
geo.geomacc = 1e-9;
end

%------------------------------------------------------------------------
function [visclass,occl] = edclassifypair(edgeA,edgeB,geo)
% 0 hidden, 1 partly visible, 2 fully visible; occl = possible occluders

acc = geo.geomacc;
A = [geo.edgestartcoords(edgeA,:);geo.edgeendcoords(edgeA,:)];
B = [geo.edgestartcoords(edgeB,:);geo.edgeendcoords(edgeB,:)];
ownplanes = [geo.planesatedge(edgeA,:) geo.planesatedge(edgeB,:)];
ownplanes = ownplanes(ownplanes > 0);
rear = geo.rearsideplane(ownplanes);
ownplanes = unique([ownplanes(:);rear(rear > 0)]);
testplanes = setdiff(geo.obstructplanes,ownplanes);
testplanes = testplanes(:).';
occl = zeros(1,0);

pA = geo.planesatedge(edgeA,:);
pB = geo.planesatedge(edgeB,:);
commonplanes = intersect(pA(pA > 0),pB(pB > 0));
sameplane = ~isempty(commonplanes);
if sameplane
    % 21 Nov. 2026: if one of the edges crosses the line of the other
    % (within their common plane), the other edge is partly on the face
    % side and partly across an indent of the plane. Then the pair is
    % partly visible, and the node pairs are tested (a path across the
    % indent passes exactly through an edge of a neighbouring plane: 1/2).
    np = geo.nvecs(commonplanes(1),:);
    if edcrossesline(A,B,np,acc) || edcrossesline(B,A,np,acc)
        visclass = 1;
        occl = testplanes;
        return
    end
    % The earlier test: the line between the midpoints. A hit inside a
    % plane, or at an edge or corner of a plane, hides the pair.
    mA = mean(A,1);
    mB = mean(B,1);
    fac = edobstructionfactor(mA,mB,geo,testplanes);
    if fac < 1
        visclass = 0;
    else
        visclass = 2;
    end
    return
end

visclass = 2;
for p = testplanes
    dA = A*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    dB = B*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    d4 = [dA;dB];
    if ~(any(d4 > acc) && any(d4 < -acc))
        continue
    end
    if (all(dA > acc) && all(dB < -acc)) || (all(dA < -acc) && all(dB > acc))
        % The cross-section: crossing points of A1-B1, A1-B2, A2-B2, A2-B1
        ia = [1 1 2 2].';
        ib = [1 2 2 1].';
        t = dA(ia)./(dA(ia) - dB(ib));
        Q = A(ia,:) + t.*(B(ib,:) - A(ia,:));
        ov = edquadpolygonoverlap(Q,p,geo);
        if ov == 2
            visclass = 0;
            occl = zeros(1,0);
            return
        elseif ov == 1
            occl = [occl p]; %#ok<AGROW>
        end
    else
        occl = [occl p]; %#ok<AGROW>
    end
end
if ~isempty(occl)
    visclass = 1;
end

end

%------------------------------------------------------------------------
function c = edcrossesline(A,B,np,acc)
% 21 Nov. 2026: true if the end points of segment B are strictly on
% different sides of the line through segment A, within the plane with
% the normal vector np.
dA = A(2,:) - A(1,:);
s1 = dot(cross(dA,B(1,:)-A(1,:)),np);
s2 = dot(cross(dA,B(2,:)-A(1,:)),np);
tol = acc*max(norm(dA),1);
c = (s1 > tol && s2 < -tol) || (s1 < -tol && s2 > tol);
end

%------------------------------------------------------------------------
function ov = edquadpolygonoverlap(Q,p,geo)
% Q [4,3] a convex quadrilateral in the plane p (corners in order).
% ov = 0: no overlap with the polygon of plane p, 1: overlap or contact
% with a boundary, 2: Q is inside the polygon.

nc = geo.ncornersperplanevec(p);
P = geo.corners(geo.planecorners(p,1:nc),:);
% Project onto the coordinate plane where the plane is largest
[~,idrop] = max(abs(geo.nvecs(p,:)));
keep = setdiff(1:3,idrop);
q = Q(:,keep);
P2 = P(:,keep);
scale = max([max(abs(P2(:)));1]);
tol = 1e-9*scale;

qin = zeros(4,1);
for ii = 1:4
    qin(ii) = edpointinpolygon(q(ii,:),P2,tol);
end
pin = zeros(nc,1);
for ii = 1:nc
    pin(ii) = edpointinpolygon(P2(ii,:),q,tol);
end
crossing = false;
for ii = 1:4
    a1 = q(ii,:);
    a2 = q(mod(ii,4)+1,:);
    for jj = 1:nc
        b1 = P2(jj,:);
        b2 = P2(mod(jj,nc)+1,:);
        if edpropercrossing(a1,a2,b1,b2,tol)
            crossing = true;
            break
        end
    end
    if crossing
        break
    end
end

if all(qin == 1) && ~crossing
    ov = 2;
elseif any(qin ~= 0) || any(pin ~= 0) || crossing
    ov = 1;
else
    ov = 0;
end

end

%------------------------------------------------------------------------
function r = edpointinpolygon(x,P,tol)
% 1: x strictly inside the polygon P [n,2], -1: on its boundary (within
% tol), 0: outside. Crossing-number test.
n = size(P,1);
for ii = 1:n
    a = P(ii,:);
    b = P(mod(ii,n)+1,:);
    ab = b - a;
    L2 = ab*ab.';
    if L2 == 0
        continue
    end
    t = ((x - a)*ab.')/L2;
    t = min(max(t,0),1);
    if norm(x - (a + t*ab)) <= tol
        r = -1;
        return
    end
end
inside = false;
for ii = 1:n
    a = P(ii,:);
    b = P(mod(ii,n)+1,:);
    if (a(2) > x(2)) ~= (b(2) > x(2))
        xc = a(1) + (x(2) - a(2))*(b(1) - a(1))/(b(2) - a(2));
        if x(1) < xc
            inside = ~inside;
        end
    end
end
r = double(inside);
end

%------------------------------------------------------------------------
function c = edpropercrossing(a1,a2,b1,b2,tol)
% true if the segments a1-a2 and b1-b2 cross at a point strictly inside
% both
o1 = edorient(a1,a2,b1);
o2 = edorient(a1,a2,b2);
o3 = edorient(b1,b2,a1);
o4 = edorient(b1,b2,a2);
la = norm(a2-a1);
lb = norm(b2-b1);
c = (o1*o2 < 0) && (o3*o4 < 0) && abs(o1) > tol*la && abs(o2) > tol*la && ...
    abs(o3) > tol*lb && abs(o4) > tol*lb;
end

function o = edorient(a,b,c)
o = (b(1)-a(1))*(c(2)-a(2)) - (b(2)-a(2))*(c(1)-a(1));
end

%------------------------------------------------------------------------
function fac = edobstructionfactor(A,B,geo,planelist)
% For segments A(i,:) -> B(i,:): 0 if a plane in planelist obstructs the
% segment, 1 if not, and 1/2 (2/3) if the segment passes exactly through
% an edge (a corner) of a plane. Only crossings where A and B are on
% different sides of a plane, by more than geomacc, count. (The same rule
% as in EDfindGApaths.)
N = size(A,1);
fac = ones(N,1);
for p = planelist
    sA = A*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    sB = B*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    iv = find( ((sA > geo.geomacc) & (sB < -geo.geomacc)) | ((sA < -geo.geomacc) & (sB > geo.geomacc)) );
    iv = iv(fac(iv) ~= 0);
    if isempty(iv)
        continue
    end
    t = sA(iv)./(sA(iv) - sB(iv));
    X = A(iv,:) + t.*(B(iv,:) - A(iv,:));
    [hitvec,edgehit,~,cornerhit] = EDpoinpla(X,p*ones(numel(iv),1),geo.minvals,geo.maxvals,...
        geo.planecorners,geo.corners,geo.ncornersperplanevec,geo.nvecs);
    fac(iv) = fac(iv).*(1 - double(hitvec) - 0.5*double(edgehit) - (1/3)*double(cornerhit));
end
end
