function [coords,nudgeinfo] = EDnudgepoints(coords,planedata,nudgedistance,geomacc)
% EDnudgepoints - Moves points that lie on the extension of a plane slightly off that plane.
% A source or receiver that lies exactly in the (infinite) plane of a
% face, but outside the face polygon, is on a visibility boundary for the
% edges of that face (they are seen at grazing incidence along the face),
% and EDSorRgeo stops with "too close to the extension of plane". Such a
% point is moved a small distance along the plane's normal (to the side
% that the face faces), and checked again, since it could then lie on the
% extension of another plane (at most 6 moves). The total field is
% continuous across the plane, so the result is the one-sided limit; it
% differs from the value on the other side by the small discontinuity
% that the edge diffraction model has there (0.03 - 0.6% in the tests).
%
% Points that lie on a face itself (inside the polygon), or on an edge or
% a corner of a face, are not moved: EDSorRgeo gives its error for those,
% as before (the two sides of a thin plate have different pressures).
%
% Input parameters:
%   coords          [npoints,3] the source or receiver coordinates
%   planedata       struct, see EDdatastructs
%   nudgedistance   (optional) the distance, in metres, of each move.
%                   Default: max(1e-6*(the size of the model), 1000*geomacc)
%   geomacc         (optional) the distance within which a point belongs to
%                   a plane, as in EDSorRgeo. Default: 1e-9 m.
%
% Output parameters:
%   coords          [npoints,3] the coordinates, with the moved points
%   nudgeinfo       [nmoves,3] one row per move: [point number, plane
%                   number, distance moved]. Empty if no point was moved.
%
% See also EDSorRgeo, EDinfrontofplane, EDpoinpla, EDmain_convex
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2026
%
% [coords,nudgeinfo] = EDnudgepoints(coords,planedata,nudgedistance,geomacc);

% 27 Nov. 2026 First version (by Claude, prompted by PS).

if nargin < 4 || isempty(geomacc)
    geomacc = 1e-9;
end
if nargin < 3 || isempty(nudgedistance)
    modelsize = max(max(planedata.corners,[],1) - min(planedata.corners,[],1));
    nudgedistance = max(1e-6*modelsize,1000*geomacc);
end
nudgeinfo = zeros(0,3);
if isempty(coords)
    return
end

nplanes = size(planedata.planeeqs,1);
nvecs = planedata.planeeqs(:,1:3);
c1 = planedata.corners(planedata.planecorners(:,1),:);
c2 = planedata.corners(planedata.planecorners(:,2),:);
skip = false(size(coords,1),1);       % points on a face, edge or corner

for iter = 1:6
    cand = find(~skip);
    if isempty(cand)
        break
    end
    % The signed distances to all planes, as in EDinfrontofplane (relative
    % to the first corner of each plane, or the second if the point is at
    % the first corner)
    P = coords(cand,:);
    np = numel(cand);
    s = zeros(np,nplanes);
    for p = 1:nplanes
        v = P - c1(p,:);
        atcorner = sqrt(sum(v.^2,2)) < geomacc;
        v(atcorner,:) = P(atcorner,:) - c2(p,:);
        s(:,p) = v*nvecs(p,:).';
    end
    [ipt,ipl] = find(abs(s) <= geomacc);
    ipt = ipt(:);
    ipl = ipl(:);
    if isempty(ipt)
        break
    end
    [hitvec,edgehit,~,cornerhit] = EDpoinpla(P(ipt,:),ipl,planedata.minvals,planedata.maxvals,...
        planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,nvecs,geomacc);
    onface = (hitvec(:) ~= 0) | (edgehit(:) ~= 0) | (cornerhit(:) ~= 0);
    % Points on a face, edge or corner of some plane: not moved
    skip(cand(unique(ipt(onface)))) = true;
    ok = ~onface & ~skip(cand(ipt));
    ipt = ipt(ok);
    ipl = ipl(ok);
    if isempty(ipt)
        break
    end
    % Each point is moved along the normal of its first such plane
    [upt,ifirst] = unique(ipt,'first');
    for ii = 1:numel(upt)
        k = cand(upt(ii));
        p = ipl(ifirst(ii));
        coords(k,:) = coords(k,:) + nudgedistance*nvecs(p,:);
        nudgeinfo = [nudgeinfo; k p nudgedistance]; %#ok<AGROW>
    end
end
end
