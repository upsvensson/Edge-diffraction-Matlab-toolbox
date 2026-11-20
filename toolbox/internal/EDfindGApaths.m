function [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
    EDfindGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
    EDversionnumber,filehandlingparameters)
% EDfindGApaths - Finds the geometrical-acoustics paths for a non-convex object.
% Finds the direct sound paths, the specular reflection paths up to the
% order controlparameters.specorder, and the first-order diffraction paths,
% possibly with specular reflections before and after the diffraction.
% For each source and receiver, the image sources and image receivers are
% found (up to the order specorder), together with the parts of each edge
% that each of them can see. These are used for the first-order
% diffraction (here) and for the source term and the propagation of the
% integral equation (via EDimagepointdata).
%
% This function replaces EDfindconvexGApaths when
% geoinputdata.allownonconvex = 1. For a convex object, the two give the
% same paths.
%
% The visible part of an edge, from a source or receiver point or an image
% of one, is found by sampling the edge (controlparameters.nedgesamples
% intervals, default 32) and refining every change between visible and not
% visible by bisection. Visible pieces that are shorter than the sampling
% interval could be missed.
%
% Input parameters:
%   planedata, edgedata, Sdata, Rdata  structs, see help EDdatastructs
%   controlparameters  struct; the fields .difforder (only whether it is
%                   larger than 0), .directsound, .specorder and
%                   .nedgesamples (optional, default 32) are used.
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem, .showtext,
%                   .suppressresultrecycling and .savepathsfile are used.
%
% Output parameters:
%   firstorderpathdata  struct, see help EDdatastructs. The same fields as
%                   from EDfindconvexGApaths (the specular reflections can
%                   be of any order up to specorder), and in addition:
%                   .specreflplanes [nIS,specorder] the reflection planes
%                   .diffpathlist   [ndiff,7] one row per first-order
%                       diffraction component: [source receiver edge Spoint
%                       Rpoint zstart zend], where Spoint and Rpoint are
%                       rows in .Spoints and .Rpoints, and zstart, zend are
%                       the visible part of the edge (fractions 0-1)
%                   .Spoints, .Rpoints  the sources/receivers and their
%                       images: .coordinates, .parent, .order, .planeseq,
%                       .amplitude (product of the reflection factors) and
%                       .edgeintervals [nint,4] = [point edge zstart zend]
%                       (zstart, zend as fractions 0-1). The first rows are
%                       the sources/receivers themselves.
%                   .specorder
%   elapsedtimefindpaths
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDfindconvexGApaths, EDimagepointdata, EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 20 Nov. 2026
%
% [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
%    EDfindGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
%    EDversionnumber,filehandlingparameters);

% 20 Nov. 2026 First version (by Claude, prompted by PS). The rules for
% which planes and edges can be seen follow EDB1findISEStree,
% EDB1speculISES and EDB1diffISESx of the EDB toolbox.

t00 = clock;

showtext = filehandlingparameters.showtext;
difforder = controlparameters.difforder;
directsound = controlparameters.directsound;
specorder = 1;
if isfield(controlparameters,'specorder') && ~isempty(controlparameters.specorder)
    specorder = controlparameters.specorder;
end
nedgesamples = 32;
if isfield(controlparameters,'nedgesamples') && ~isempty(controlparameters.nedgesamples)
    nedgesamples = controlparameters.nedgesamples;
end

stageversion = EDstageversions('EDfindGApaths');

EDinputdatastruct = struct('stage','EDfindGApaths','stageversion',stageversion,...
    'edgedatahash',EDparenthash(edgedata),'Sdatahash',EDparenthash(Sdata),...
    'Rdatahash',EDparenthash(Rdata),'calcdiffpaths',double(difforder > 0),...
    'directsound',directsound,'doallSRcombinations',Sdata.doallSRcombinations,...
    'specorder',specorder,'nedgesamples',nedgesamples);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Can an existing file be used?

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_paths',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_paths.mat'];

if foundmatch == 1
    filedata = load(existingfilename,'firstorderpathdata','elapsedtimefindpaths');
    firstorderpathdata = filedata.firstorderpathdata;
    elapsedtimefindpaths = filedata.elapsedtimefindpaths;
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    elapsedtimefindpaths_new = etime(clock,t00);
    elapsedtimefindpaths = [elapsedtimefindpaths_new elapsedtimefindpaths];
    return
end

%-----------------------------------------------------------------
if isfield(Sdata,'sourcetype') && ~strcmp(Sdata.sourcetype,'monopole')
    error('ERROR: EDfindGApaths (allownonconvex = 1) can only handle monopole sources')
end

nsources = size(Sdata.coordinates,1);
nreceivers = size(Rdata.coordinates,1);
nedges = size(edgedata.edgecorners,1);

geo = edgeomdata(planedata,edgedata);

%-----------------------------------------------------------------
% The direct sound is found in the same way as for a convex object: the
% obstruction test checks every plane that has the source and the receiver
% on different sides. EDfindconvexGApaths is used for that, with no file
% saving.

if directsound ~= 0
    cpdirect = controlparameters;
    cpdirect.difforder = 0;
    fhdirect = filehandlingparameters;
    fhdirect.suppressresultrecycling = 1;
    fhdirect.savepathsfile = 0;
    fhdirect.showtext = 0;
    directpaths = EDfindconvexGApaths(planedata,edgedata,Sdata,Rdata,cpdirect,...
        EDversionnumber,fhdirect);
    directsoundlist = directpaths.directsoundlist;
else
    directsoundlist = zeros(0,3);
end

%-----------------------------------------------------------------
% The image sources and image receivers, up to order specorder

if showtext >= 2
    disp(['      Finding image sources and image receivers up to order ',int2str(specorder)])
end
Spoints = edimagetree(Sdata.coordinates,geo,specorder);
Rpoints = edimagetree(Rdata.coordinates,geo,specorder);

%-----------------------------------------------------------------
% Specular reflections: every image source (order >= 1) to every receiver

if showtext >= 2
   disp(['      Finding specular reflection paths'])
end
specrefllist = zeros(0,3);
specreflIScoords = zeros(0,3);
specreflplanes = zeros(0,specorder);
ivIS = find(Spoints.order >= 1 & Spoints.amplitude ~= 0);
for jj = 1:numel(ivIS)
    row = ivIS(jj);
    isou = Spoints.parent(row);
    if Sdata.doallSRcombinations == 1
        reclist = (1:nreceivers).';
    elseif isou <= nreceivers
        reclist = isou;
    else
        reclist = [];
    end
    if isempty(reclist)
        continue
    end
    fac = edchainfactor(Spoints,row,Rdata.coordinates(reclist,:),geo);
    ivok = find(fac ~= 0);
    if ~isempty(ivok)
        nok = numel(ivok);
        specrefllist = [specrefllist;[isou*ones(nok,1) reclist(ivok) fac(ivok)*Spoints.amplitude(row)]];
        specreflIScoords = [specreflIScoords;repmat(Spoints.coordinates(row,:),nok,1)];
        specreflplanes = [specreflplanes;repmat(Spoints.planeseq(row,:),nok,1)];
    end
end

%-----------------------------------------------------------------
% First-order diffraction: the visible parts of each edge from every
% source/image source and every receiver/image receiver, and then the
% combinations.

if difforder > 0
    if showtext >= 2
       disp(['      Finding the visible parts of the edges'])
    end
    activeedges = (1:nedges).';
    activeedges(edgedata.offedges) = [];
    Spoints.edgeintervals = ededgeintervals(Spoints,activeedges,geo,nedgesamples);
    Rpoints.edgeintervals = ededgeintervals(Rpoints,activeedges,geo,nedgesamples);

    if showtext >= 2
       disp(['      Finding first-order diffraction paths'])
    end
    diffpathlist = zeros(0,7);
    diffpaths = zeros(nreceivers,nsources,nedges);
    for kk = 1:numel(activeedges)
        edgenumber = activeedges(kk);
        Sint = Spoints.edgeintervals(Spoints.edgeintervals(:,2)==edgenumber,:);
        Rint = Rpoints.edgeintervals(Rpoints.edgeintervals(:,2)==edgenumber,:);
        if isempty(Sint) || isempty(Rint)
            continue
        end
        for is = 1:size(Sint,1)
            isou = Spoints.parent(Sint(is,1));
            for ir = 1:size(Rint,1)
                irec = Rpoints.parent(Rint(ir,1));
                if Sdata.doallSRcombinations == 0 && isou ~= irec
                    continue
                end
                zstart = max(Sint(is,3),Rint(ir,3));
                zend = min(Sint(is,4),Rint(ir,4));
                if zend > zstart
                    diffpathlist = [diffpathlist;[isou irec edgenumber Sint(is,1) Rint(ir,1) zstart zend]];
                    diffpaths(irec,isou,edgenumber) = 1;
                end
            end
        end
    end
    edgeisactive = zeros(nedges,1,'uint8');
    if ~isempty(diffpathlist)
        edgeisactive(unique(diffpathlist(:,3))) = 1;
    end
else
    Spoints.edgeintervals = zeros(0,4);
    Rpoints.edgeintervals = zeros(0,4);
    diffpathlist = zeros(0,7);
    diffpaths = [];
    edgeisactive = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pack the output data in a struct

firstorderpathdata = struct;
firstorderpathdata.specreflIScoords = specreflIScoords;
firstorderpathdata.specrefllist     = specrefllist;
firstorderpathdata.specreflplanes   = specreflplanes;
firstorderpathdata.diffpaths        = diffpaths;
firstorderpathdata.edgeisactive     = edgeisactive;
firstorderpathdata.directsoundlist  = directsoundlist;
firstorderpathdata.diffpathlist     = diffpathlist;
firstorderpathdata.Spoints          = Spoints;
firstorderpathdata.Rpoints          = Rpoints;
firstorderpathdata.specorder        = specorder;
firstorderpathdata.ncomponents      = [size(directsoundlist,1) size(specrefllist,1) size(diffpathlist,1)];

firstorderpathdata.EDinputdatahash = EDinputdatahash;

elapsedtimefindpaths = etime(clock,t00);

if filehandlingparameters.savepathsfile == 1
	save(desiredname,'firstorderpathdata','EDinputdatahash','elapsedtimefindpaths','stageversion','EDversionnumber');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function geo = edgeomdata(planedata,edgedata)
% The geometry data that the local functions need, in one struct.

geo = struct;
geo.planeeqs = planedata.planeeqs;
geo.nvecs = planedata.planeeqs(:,1:3);
geo.minvals = planedata.minvals;
geo.maxvals = planedata.maxvals;
geo.planecorners = planedata.planecorners;
geo.corners = planedata.corners;
geo.ncornersperplanevec = planedata.ncornersperplanevec;
geo.reflfactors = planedata.reflfactors(:);
geo.rearsideplane = planedata.rearsideplane(:);
geo.obstructplanes = find(planedata.canplaneobstruct(:) ~= 0);
geo.planesatedge = double(edgedata.planesatedge);
geo.closwedangvec = edgedata.closwedangvec;
geo.edgestartcoords = edgedata.edgestartcoords;
geo.edgeendcoords = edgedata.edgeendcoords;
% Points closer to a plane than geomacc count as in the plane
geo.geomacc = 1e-9;

end

%------------------------------------------------------------------------
function P = edimagetree(basecoords,geo,specorder)
% The base points (sources or receivers) and all their image points up to
% order specorder. A point is mirrored in a plane if it is in front of
% that plane, the plane is not totally absorbing, and the plane is not the
% one of the previous reflection.

nbase = size(basecoords,1);
coords = basecoords;
parent = (1:nbase).';
order = zeros(nbase,1);
planeseq = zeros(nbase,specorder);
previous = zeros(nbase,1);
amplitude = ones(nbase,1);

activeplanes = find(geo.reflfactors ~= 0).';
lastorder = (1:nbase).';
for norder = 1:specorder
    newrows = [];
    for row = lastorder.'
        lastplane = 0;
        if norder > 1
            lastplane = planeseq(row,norder-1);
        end
        for p = activeplanes
            if p == lastplane
                continue
            end
            s = coords(row,:)*geo.nvecs(p,:).' - geo.planeeqs(p,4);
            if s > geo.geomacc
                newcoords = coords(row,:) - 2*s*geo.nvecs(p,:);
                coords = [coords;newcoords];
                parent = [parent;parent(row)];
                order = [order;norder];
                newseq = planeseq(row,:);
                newseq(norder) = p;
                planeseq = [planeseq;newseq];
                previous = [previous;row];
                amplitude = [amplitude;amplitude(row)*geo.reflfactors(p)];
                newrows = [newrows;size(coords,1)];
            end
        end
    end
    lastorder = newrows;
    if isempty(lastorder)
        break
    end
end

P = struct('coordinates',coords,'parent',parent,'order',order,...
    'planeseq',planeseq,'previous',previous,'amplitude',amplitude);

end

%------------------------------------------------------------------------
function fac = edchainfactor(P,row,T,geo)
% For the image point in row 'row' of P (or the point itself, if its order
% is 0), and a number of target points T [M,3], the factor of the specular
% path from the base point, via the reflection planes, to each target:
% 0 if the path is not valid or obstructed, 1 if it is valid, and 1/2
% (edge hit) or 1/3 (corner hit) when a reflection point is right at the
% boundary of its plane. Obstructions right at an edge give 1/2, at a
% corner 2/3, as for the direct sound in EDfindconvexGApaths.

M = size(T,1);
fac = ones(M,1);

% The rows of the chain, from this image point back to the base point
nrefl = P.order(row);
chainrows = zeros(nrefl+1,1);
chainrows(nrefl+1) = row;
for jj = nrefl:-1:1
    chainrows(jj) = P.previous(chainrows(jj+1));
end
basepoint = P.coordinates(chainrows(1),:);

current = T;
for jj = nrefl:-1:1
    iv = find(fac ~= 0);
    if isempty(iv)
        return
    end
    p = P.planeseq(row,jj);
    I = P.coordinates(chainrows(jj+1),:);
    sC = current(iv,:)*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    sI = I*geo.nvecs(p,:).' - geo.planeeqs(p,4);
    ok = (sC > geo.geomacc) & (sI < -geo.geomacc);
    fac(iv(~ok)) = 0;
    iv = iv(ok);
    if isempty(iv)
        return
    end
    t = sC(ok)./(sC(ok) - sI);
    X = current(iv,:) + t.*(I - current(iv,:));
    [hitvec,edgehit,~,cornerhit] = EDpoinpla(X,p*ones(numel(iv),1),geo.minvals,geo.maxvals,...
        geo.planecorners,geo.corners,geo.ncornersperplanevec,geo.nvecs);
    hitfactor = double(hitvec) + 0.5*double(edgehit) + (1/3)*double(cornerhit);
    fac(iv) = fac(iv).*hitfactor;
    % The path segment from the reflection point to the current point
    ivv = find(fac(iv) ~= 0);
    if ~isempty(ivv)
        fac(iv(ivv)) = fac(iv(ivv)).*edobstructionfactor(X(ivv,:),current(iv(ivv),:),geo);
    end
    current(iv,:) = X;
end

% The path segment from the base point to the first reflection point (or
% to the target, if there is no reflection)
iv = find(fac ~= 0);
if ~isempty(iv)
    fac(iv) = fac(iv).*edobstructionfactor(repmat(basepoint,numel(iv),1),current(iv,:),geo);
end

end

%------------------------------------------------------------------------
function fac = edobstructionfactor(A,B,geo)
% For segments A(i,:) -> B(i,:): 0 if a plane obstructs the segment, 1 if
% not, and 1/2 (2/3) if the segment passes exactly through an edge (a
% corner) of a plane. Only crossings where A and B are on different sides
% of a plane, by more than geomacc, count. So planes that contain one of
% the end points (reflection planes, the planes of an edge) never
% obstruct.

N = size(A,1);
fac = ones(N,1);
for p = geo.obstructplanes.'
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

%------------------------------------------------------------------------
function intervals = ededgeintervals(P,edgelist,geo,nsamples)
% The visible parts of the edges in edgelist, from every point (and image
% point) in P. Output: [nint,4] = [row of P, edge, zstart, zend], with
% zstart and zend as fractions (0-1) of the edge length.
%
% A point sees an edge if it is in front of the edge's planes (in front of
% one of them for an open wedge angle larger than pi, both for one smaller
% than pi), and the path to each edge point is valid and unobstructed. For
% an image point, the last reflection plane can not be one of the edge's
% own planes (that reflection is part of the edge's own diffraction).

intervals = zeros(0,4);
u = linspace(0,1,nsamples+1).';
% Nudge the end points a little bit inwards, to avoid hits exactly at
% corners
u(1) = 1e-6;
u(end) = 1-1e-6;
nbisect = 40;

for row = 1:size(P.coordinates,1)
    if P.amplitude(row) == 0
        continue
    end
    point = P.coordinates(row,:);
    nrefl = P.order(row);
    lastplane = 0;
    if nrefl > 0
        lastplane = P.planeseq(row,nrefl);
    end
    for edgenumber = edgelist(:).'
        p1 = geo.planesatedge(edgenumber,1);
        p2 = geo.planesatedge(edgenumber,2);
        if lastplane ~= 0
            ownplanes = [p1 p2];
            ownplanes = [ownplanes geo.rearsideplane(ownplanes(ownplanes>0)).'];
            if any(ownplanes == lastplane)
                continue
            end
            % The edge must be (partly) in front of the last reflection
            % plane
            sedge = [geo.edgestartcoords(edgenumber,:);geo.edgeendcoords(edgenumber,:)]*geo.nvecs(lastplane,:).' - geo.planeeqs(lastplane,4);
            if max(sedge) <= geo.geomacc
                continue
            end
        end
        s1 = point*geo.nvecs(p1,:).' - geo.planeeqs(p1,4);
        s2 = point*geo.nvecs(p2,:).' - geo.planeeqs(p2,4);
        if geo.closwedangvec(edgenumber) < pi
            infront = (s1 > geo.geomacc) || (s2 > geo.geomacc);
        else
            infront = (s1 > geo.geomacc) && (s2 > geo.geomacc);
        end
        if ~infront
            continue
        end
        e1 = geo.edgestartcoords(edgenumber,:);
        e2 = geo.edgeendcoords(edgenumber,:);
        isvisible = edchainfactor(P,row,e1 + u*(e2-e1),geo) ~= 0;
        if ~any(isvisible)
            continue
        end
        % Refine every change between visible and not visible
        ivchange = find(diff(isvisible));
        ulow = u(ivchange);
        uhigh = u(ivchange+1);
        vislow = isvisible(ivchange);
        for kk = 1:nbisect
            umid = (ulow + uhigh)/2;
            vismid = edchainfactor(P,row,e1 + umid*(e2-e1),geo) ~= 0;
            same = (vismid == vislow);
            ulow(same) = umid(same);
            uhigh(~same) = umid(~same);
        end
        uchange = (ulow + uhigh)/2;
        % Build the intervals
        startpoints = [];
        endpoints = [];
        if isvisible(1)
            startpoints = 0;
        end
        for kk = 1:numel(ivchange)
            if vislow(kk)
                endpoints = [endpoints;uchange(kk)];
            else
                startpoints = [startpoints;uchange(kk)];
            end
        end
        if isvisible(end)
            endpoints = [endpoints;1];
        end
        nint = numel(startpoints);
        intervals = [intervals;[row*ones(nint,1) edgenumber*ones(nint,1) startpoints(:) endpoints(:)]];
    end
end

end
