function slopeplanes = EDslopediffractionplanes(planedata,edgedata,settings)
% EDslopediffractionplanes - Finds the planes that are hit by diffraction waves through free space.
% The edge source formulation is exact when an edge is hit by the wave
% from a monopole, or by a diffraction wave that travels along a rigid
% surface (the edges of one plane). When an edge is hit by a diffraction
% wave through free space, the missing slope diffraction makes the
% second- and higher-order diffraction approximate, and the total field
% does not fulfill the rigid boundary condition on the surfaces that this
% wave hits. A second pass with secondary sources on those surfaces can
% correct this (see EDdatastructs, slope diffraction handling).
%
% A plane P needs such a correction if some part of it is hit by a
% diffraction wave through free space, that is, if some part of it is
% visible from an edge E that
%   - is active (not switched off, and not flat: an open angle of pi does
%     not diffract),
%   - does not lie in the plane of P (or of a plane that is coplanar with
%     P; those edges only send waves along the surface of P),
%   - is, at least partly, in front of P,
% and the path from E to the point of P
%   - leaves E strictly inside its open wedge (a path along one of E's
%     faces travels along a surface), and
%   - is not obstructed by any plane.
% For a convex polyhedron no plane is found: each plane sees only the
% edges in its own plane. For a folded plate (an L-shape) the two inner
% faces are found, and the two outer faces are not.
%
% Each plane is sampled with a grid of points (spacing settings.spacing),
% and each candidate edge with settings.nedgesamples points. A plane
% sample is "hit" if it is visible from at least one edge sample.
%
% Input parameters:
%   planedata, edgedata     structs, see EDdatastructs
%   settings        (optional) struct with the fields
%     .spacing      the spacing of the plane samples, in metres. Default:
%                   1/10 of the largest side of the bounding box of the
%                   model.
%     .nedgesamples the number of samples per edge. Default: 9.
%     .showtext     0 (default) or 1
%
% Output parameters:
%   slopeplanes     struct with the fields
%     .planelist        [nsp,1] the planes that need slope diffraction
%                       handling (a part of them is hit)
%     .illuminatingedges {nsp,1} the edges that hit each of those planes
%     .planeedgehit     [nplanes,nedges] logical, true if edge E hits some
%                       part of plane P through free space
%     .samples          {nplanes,1} the sample points of each plane [ns,3]
%     .samplehitcount   {nplanes,1} [ns,1] the number of edges that hit
%                       each sample point (0 = not hit)
%     .hitfraction      [nplanes,1] the fraction of the sample points of
%                       each plane that are hit
%     .spacing, .nedgesamples   the values that were used
%
% See also EDed2geo, EDedgepairvisibility, EDpoinpla, EDcoordtrans1
%
% Peter Svensson (peter.svensson@ntnu.no) 22 Nov. 2026
%
% slopeplanes = EDslopediffractionplanes(planedata,edgedata,settings);

% 22 Nov. 2026 First version (by Claude, prompted by PS).

if nargin < 3 || isempty(settings)
    settings = struct;
end
corners = planedata.corners;
if ~isfield(settings,'spacing') || isempty(settings.spacing)
    settings.spacing = max(max(corners,[],1) - min(corners,[],1))/10;
end
if ~isfield(settings,'nedgesamples') || isempty(settings.nedgesamples)
    settings.nedgesamples = 9;
end
if ~isfield(settings,'showtext') || isempty(settings.showtext)
    settings.showtext = 0;
end

geo = edgeomdata(planedata);
nplanes = size(planedata.planeeqs,1);
nedges = size(edgedata.edgestartcoords,1);
acc = geo.geomacc;
anglemargin = 1e-6;

% Coplanar planes (the same plane equation, the same side)
pe = planedata.planeeqs;
coplanar = false(nplanes);
for pp = 1:nplanes
    coplanar(:,pp) = all(abs(pe(:,1:3) - pe(pp,1:3)) < 1e-9,2) & abs(pe(:,4) - pe(pp,4)) < 1e-9*max(1,abs(pe(pp,4)));
end
coplanar(1:nplanes+1:end) = true;

% Active edges: not switched off, and not flat
activeedges = true(nedges,1);
if isfield(edgedata,'offedges') && ~isempty(edgedata.offedges)
    activeedges(edgedata.offedges) = false;
end
activeedges(abs(edgedata.closwedangvec(:) - pi) < 1e-9) = false;

% Edge samples
u = ((1:settings.nedgesamples).' - 0.5)/settings.nedgesamples;
edgesamples = cell(nedges,1);
for e = 1:nedges
    edgesamples{e} = edgedata.edgestartcoords(e,:) + u*(edgedata.edgeendcoords(e,:) - edgedata.edgestartcoords(e,:));
end

samples = edplanesamples(geo,settings.spacing);
samplehitcount = cell(nplanes,1);
planeedgehit = false(nplanes,nedges);

for P = 1:nplanes
    S = samples{P};
    ns = size(S,1);
    hitcount = zeros(ns,1);
    nvP = geo.nvecs(P,:);
    for E = find(activeedges).'
        % Not an edge in the plane of P (or of a plane coplanar with P)
        if any(coplanar(double(edgedata.planesatedge(E,:)),P))
            continue
        end
        Q = edgesamples{E};
        dQ = Q*nvP.' - geo.planeeqs(P,4);
        Q = Q(dQ > acc,:);                 % the edge samples in front of P
        if isempty(Q) || ns == 0
            continue
        end
        % The plane samples strictly inside the open wedge of E
        openangle = 2*pi - edgedata.closwedangvec(E);
        [~,thS] = EDcoordtrans1(S,[edgedata.edgestartcoords(E,:);edgedata.edgeendcoords(E,:)],...
            edgedata.edgenvecs(E,:),reshape(edgedata.edgerelatedcoordsysmatrices(E,:),3,3));
        thS = thS(:);
        inwedge = thS > anglemargin & thS < openangle - anglemargin;
        if ~any(inwedge)
            continue
        end
        % Unobstructed paths from the edge samples to the plane samples
        iS = find(inwedge);
        [ia,ib] = ndgrid(1:size(Q,1),iS);
        fac = edobstructionfactor(Q(ia(:),:),S(ib(:),:),geo);
        seen = reshape(fac > 0,size(Q,1),numel(iS));
        hitS = false(ns,1);
        hitS(iS) = any(seen,1).';
        if any(hitS)
            planeedgehit(P,E) = true;
            hitcount = hitcount + double(hitS);
        end
    end
    samplehitcount{P} = hitcount;
end

planelist = find(any(planeedgehit,2));
illuminatingedges = cell(numel(planelist),1);
for ii = 1:numel(planelist)
    illuminatingedges{ii} = find(planeedgehit(planelist(ii),:)).';
end
hitfraction = zeros(nplanes,1);
for P = 1:nplanes
    if ~isempty(samplehitcount{P})
        hitfraction(P) = mean(samplehitcount{P} > 0);
    end
end

slopeplanes = struct('planelist',planelist,'illuminatingedges',{illuminatingedges},...
    'planeedgehit',planeedgehit,'samples',{samples},'samplehitcount',{samplehitcount},...
    'hitfraction',hitfraction,'spacing',settings.spacing,'nedgesamples',settings.nedgesamples);

if settings.showtext >= 1
    if isempty(planelist)
        disp('      Slope diffraction: no plane is hit by diffraction waves through free space')
    else
        for ii = 1:numel(planelist)
            disp(['      Slope diffraction: plane ',int2str(planelist(ii)),' (',...
                num2str(100*hitfraction(planelist(ii)),3),'% of it) is hit by edges ',...
                int2str(illuminatingedges{ii}.')])
        end
    end
end

end

%------------------------------------------------------------------------
function geo = edgeomdata(planedata)
geo = struct;
geo.planeeqs = planedata.planeeqs;
geo.nvecs = planedata.planeeqs(:,1:3);
geo.minvals = planedata.minvals;
geo.maxvals = planedata.maxvals;
geo.planecorners = planedata.planecorners;
geo.corners = planedata.corners;
geo.ncornersperplanevec = planedata.ncornersperplanevec;
obstructplanes = (1:size(planedata.planeeqs,1)).';
if isfield(planedata,'canplaneobstruct')
    obstructplanes = find(planedata.canplaneobstruct(:) ~= 0);
end
% A thin plate has two planes with the same polygon: test only one of them
if isfield(planedata,'rearsideplane')
    rs = planedata.rearsideplane(obstructplanes);
    rs = rs(:);
    obstructplanes = obstructplanes(~(rs > 0 & ismember(rs,obstructplanes) & rs < obstructplanes));
end
geo.obstructplanes = obstructplanes(:).';
geo.geomacc = 1e-9;
end

%------------------------------------------------------------------------
function samples = edplanesamples(geo,spacing)
% Sample points inside each plane polygon: a grid with the given spacing
% (cell centres, at least 3 x 3), moved slightly off the polygon sides.
nplanes = size(geo.planeeqs,1);
samples = cell(nplanes,1);
for p = 1:nplanes
    nc = geo.ncornersperplanevec(p);
    C = geo.corners(geo.planecorners(p,1:nc),:);
    c0 = mean(C,1);
    nv = geo.nvecs(p,:);
    e1 = C(2,:) - C(1,:);
    e1 = e1/norm(e1);
    e2 = cross(nv,e1);
    uv = (C - c0)*[e1.' e2.'];
    umin = min(uv(:,1)); umax = max(uv(:,1));
    vmin = min(uv(:,2)); vmax = max(uv(:,2));
    nu = max(3,ceil((umax-umin)/spacing));
    nw = max(3,ceil((vmax-vmin)/spacing));
    [ug,vg] = meshgrid(umin + ((1:nu)-0.5)/nu*(umax-umin), vmin + ((1:nw)-0.5)/nw*(vmax-vmin));
    X = c0 + ug(:)*e1 + vg(:)*e2;
    hitvec = EDpoinpla(X,p*ones(size(X,1),1),geo.minvals,geo.maxvals,...
        geo.planecorners,geo.corners,geo.ncornersperplanevec,geo.nvecs);
    samples{p} = X(hitvec ~= 0,:);
end
end

%------------------------------------------------------------------------
function fac = edobstructionfactor(A,B,geo)
% For segments A(i,:) -> B(i,:): 0 if a plane obstructs the segment, 1 if
% not, and 1/2 (2/3) if the segment passes exactly through an edge (a
% corner) of a plane. Only crossings where A and B are on different sides
% of a plane, by more than geomacc, count (as in EDfindGApaths), so the
% planes that contain an end point never obstruct.
N = size(A,1);
fac = ones(N,1);
for p = geo.obstructplanes
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
