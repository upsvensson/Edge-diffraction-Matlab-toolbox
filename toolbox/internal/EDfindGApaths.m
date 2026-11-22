function [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
    EDfindGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
    EDversionnumber,filehandlingparameters)
% EDfindGApaths - Finds the geometrical-acoustics paths for a non-convex object.
% Finds the direct sound paths, the specular reflection paths up to the
% order controlparameters.specorder, and the first-order diffraction paths,
% possibly with specular reflections before and after the diffraction.
% For each source and receiver, the image sources and image receivers are
% found (up to the order specorder), together with the parts of each edge
% that each of them can see.
%
% specorder can be a number N or 'auto'. With N, all image points up to
% order N are used (a point is mirrored in every plane that it is in front
% of), and a warning is given if image points of order N+1 would see a
% receiver (source) or a part of an edge: the diffraction from the image
% points of order N assumes the reflections of order N+1, so the results
% could have discontinuities. With 'auto', the order is raised one step at
% a time, and a new image point is kept only if its path reaches a
% receiver (source), a part of an edge or a part of a plane (sample
% points: corners, side midpoints, the centroid and a 6x6 grid). It stops
% when an order adds nothing, or at controlparameters.specordermax
% (default 8), with a warning if order specordermax+1 would still
% contribute (e.g. parallel planes that face each other). These are used for the first-order
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
%                   larger than 0), .directsound, .specorder, .specordermax
%                   (only for 'auto', default 8) and .nedgesamples
%                   (optional, default 32) are used.
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
%                   .specreflplanes [nIS,specorderused] the reflection planes
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
%                   .specorder      as in controlparameters (N or 'auto')
%                   .specorderused  the highest order of the image points
%                   .specordertruncated 1 if image points of a higher order
%                       would still contribute (a warning was given), 0 if
%                       not
%   elapsedtimefindpaths
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDfindconvexGApaths, EDimagepointdata, EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 21 Nov. 2026
%
% [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
%    EDfindGApaths(planedata,edgedata,Sdata,Rdata,controlparameters,...
%    EDversionnumber,filehandlingparameters);

% 20 Nov. 2026 First version (by Claude, prompted by PS). The rules for
% which planes and edges can be seen follow EDB1findISEStree,
% EDB1speculISES and EDB1diffISESx of the EDB toolbox.
% 20 Nov. 2026 (by Claude, prompted by PS) For a thin plate, only one of
% its two planes is tested for obstruction (a path exactly through an
% edge of the plate got the factor 1/4 instead of 1/2). Version 2.
% 21 Nov. 2026 (by Claude, prompted by PS) specorder = 'auto' (the new
% default) and controlparameters.specordermax; for a number N, a warning
% if order N+1 would contribute (found with a 100-degree folded plate and
% specorder = 1, where a sector of receivers got the total field 0 instead
% of 1). New fields .specorderused and .specordertruncated. Version 3.
% 22 Nov. 2026 (by Claude, prompted by PS) The plane sequences of the
% specular paths are padded to specorderused columns (an error when the
% image receivers reached a higher order than the image sources).

t00 = clock;

showtext = filehandlingparameters.showtext;
difforder = controlparameters.difforder;
directsound = controlparameters.directsound;
specorder = 1;
if isfield(controlparameters,'specorder') && ~isempty(controlparameters.specorder)
    specorder = controlparameters.specorder;
end
autoorder = ischar(specorder) && strcmp(specorder,'auto');
specordermax = 8;
if isfield(controlparameters,'specordermax') && ~isempty(controlparameters.specordermax)
    specordermax = controlparameters.specordermax;
end
if ~autoorder
    specordermax = [];
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
    'specorder',specorder,'specordermax',specordermax,'nedgesamples',nedgesamples);
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
% The image sources and image receivers.
% specorder = N: all image points up to order N (a point is mirrored in
% every plane that it is in front of). Afterwards it is checked whether
% image points of order N+1 would contribute; if so, a warning is given.
% specorder = 'auto' (21 Nov. 2026): the order is raised one step at a
% time, and a new image point is kept only if its reflection path reaches
% a receiver (or source, for the image receivers), a part of an edge, or
% a part of a plane (sample points). The tree stops growing when an order
% gives no new image points, or at the order specordermax (with a
% warning).

activeedges = (1:nedges).';
activeedges(edgedata.offedges) = [];
if difforder == 0
    edgesforpruning = zeros(0,1);
else
    edgesforpruning = activeedges;
end
if autoorder
    if showtext >= 2
        disp(['      Finding image sources and image receivers (specorder = ''auto'', at most ',int2str(specordermax),')'])
    end
    planesamples = edplanesamples(geo);
    [Spoints,Sgrowing] = edimagetreeauto(Sdata.coordinates,geo,specordermax,Rdata.coordinates,edgesforpruning,nedgesamples,planesamples);
    [Rpoints,Rgrowing] = edimagetreeauto(Rdata.coordinates,geo,specordermax,Sdata.coordinates,edgesforpruning,nedgesamples,planesamples);
    specorderused = max([0;Spoints.order;Rpoints.order]);
    specordertruncated = Sgrowing || Rgrowing;
    if specordertruncated
        warning('EDfindGApaths:specordermax',['specorder = ''auto'' reached specordermax = ',int2str(specordermax),...
            ', and image points of order ',int2str(specordermax+1),' would still contribute. ',...
            'The results could have discontinuities. Increase controlparameters.specordermax.',edwedgetext(geo)])
    end
else
    if showtext >= 2
        disp(['      Finding image sources and image receivers up to order ',int2str(specorder)])
    end
    Spoints = edimagetree(Sdata.coordinates,geo,specorder);
    Rpoints = edimagetree(Rdata.coordinates,geo,specorder);
    specorderused = specorder;
    planesamples = edplanesamples(geo);
    nextS = edtruncationcheck(Spoints,specorder,geo,Rdata.coordinates,edgesforpruning,nedgesamples,planesamples);
    nextR = edtruncationcheck(Rpoints,specorder,geo,Sdata.coordinates,edgesforpruning,nedgesamples,planesamples);
    specordertruncated = (nextS + nextR) > 0;
    if specordertruncated
        warning('EDfindGApaths:specorder',['With specorder = ',int2str(specorder),', ',int2str(nextS),...
            ' image sources and ',int2str(nextR),' image receivers of order ',int2str(specorder+1),...
            ' would also contribute (they see a receiver/source or a part of an edge). ',...
            'The results could have discontinuities (diffraction from the image points of order ',int2str(specorder),...
            ' assumes the reflections of order ',int2str(specorder+1),'). Use specorder = ''auto'' or a higher value.',edwedgetext(geo)])
    end
end

%-----------------------------------------------------------------
% Specular reflections: every image source (order >= 1) to every receiver

if showtext >= 2
   disp(['      Finding specular reflection paths'])
end
specrefllist = zeros(0,3);
specreflIScoords = zeros(0,3);
specreflplanes = zeros(0,specorderused);
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
        % 22 Nov. 2026: padded to specorderused columns (Spoints.planeseq can
        % have fewer columns, when only the image receivers reach that order)
        ps = zeros(1,specorderused);
        ps(1:size(Spoints.planeseq,2)) = Spoints.planeseq(row,:);
        specreflplanes = [specreflplanes;repmat(ps,nok,1)];
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
firstorderpathdata.specorderused    = specorderused;
firstorderpathdata.specordertruncated = double(specordertruncated);
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
% A thin plate has two planes with the same polygon (each other's
% rearsideplane): only one of them is tested, so that a path exactly
% through an edge of the plate gets the factor 1/2 once (20 Nov. 2026)
rs = planedata.rearsideplane(geo.obstructplanes);
rs = rs(:);
geo.obstructplanes = geo.obstructplanes(~(rs > 0 & ismember(rs,geo.obstructplanes) & rs < geo.obstructplanes));
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

%------------------------------------------------------------------------
function [P,stillgrowing] = edimagetreeauto(basecoords,geo,maxorder,targets,edgelist,nedgesamples,planesamples)
% The base points and their image points, with the order raised one step
% at a time (specorder = 'auto'). A new image point is kept only if it
% sees something (edseessomething). stillgrowing = true if image points
% of order maxorder+1 would still be kept.

nbase = size(basecoords,1);
P = struct('coordinates',basecoords,'parent',(1:nbase).','order',zeros(nbase,1),...
    'planeseq',zeros(nbase,maxorder+1),'previous',zeros(nbase,1),'amplitude',ones(nbase,1));
activeplanes = find(geo.reflfactors ~= 0).';
lastorder = (1:nbase).';
stillgrowing = false;
for norder = 1:maxorder+1
    newrows = [];
    for row = lastorder.'
        lastplane = 0;
        if norder > 1
            lastplane = P.planeseq(row,norder-1);
        end
        for p = activeplanes
            if p == lastplane
                continue
            end
            sdist = P.coordinates(row,:)*geo.nvecs(p,:).' - geo.planeeqs(p,4);
            if sdist <= geo.geomacc
                continue
            end
            newseq = P.planeseq(row,:);
            newseq(norder) = p;
            P.coordinates = [P.coordinates;P.coordinates(row,:) - 2*sdist*geo.nvecs(p,:)];
            P.parent = [P.parent;P.parent(row)];
            P.order = [P.order;norder];
            P.planeseq = [P.planeseq;newseq];
            P.previous = [P.previous;row];
            P.amplitude = [P.amplitude;P.amplitude(row)*geo.reflfactors(p)];
            newrow = size(P.coordinates,1);
            keep = edseessomething(P,newrow,geo,targets,edgelist,nedgesamples,planesamples,true);
            if norder > maxorder && keep
                stillgrowing = true;
            end
            if keep && norder <= maxorder
                newrows = [newrows;newrow];
            else
                P.coordinates(newrow,:) = [];
                P.parent(newrow) = [];
                P.order(newrow) = [];
                P.planeseq(newrow,:) = [];
                P.previous(newrow) = [];
                P.amplitude(newrow) = [];
            end
            if stillgrowing
                break
            end
        end
        if stillgrowing
            break
        end
    end
    lastorder = newrows;
    if isempty(lastorder) || stillgrowing
        break
    end
end
usedorder = max([0;P.order]);
P.planeseq = P.planeseq(:,1:max(usedorder,1));
if usedorder == 0
    P.planeseq = zeros(nbase,0);
end

end

%------------------------------------------------------------------------
function nnext = edtruncationcheck(P,specorder,geo,targets,edgelist,nedgesamples,planesamples)
% For specorder = N: the number of image points of order N+1 that would
% see a target (receiver or source) or a part of an edge. Only the
% children of order-N points that reach a part of some plane are tried.

nnext = 0;
rows = find(P.order == specorder & P.amplitude ~= 0);
if isempty(rows)
    return
end
activeplanes = find(geo.reflfactors ~= 0).';
for row = rows.'
    if specorder > 0 && ~edseessomething(P,row,geo,zeros(0,3),zeros(0,1),nedgesamples,planesamples,true)
        continue
    end
    lastplane = 0;
    if specorder > 0
        lastplane = P.planeseq(row,specorder);
    end
    for p = activeplanes
        if p == lastplane
            continue
        end
        sdist = P.coordinates(row,:)*geo.nvecs(p,:).' - geo.planeeqs(p,4);
        if sdist <= geo.geomacc
            continue
        end
        Q = P;
        newseq = zeros(1,specorder+1);
        newseq(1:specorder) = P.planeseq(row,1:specorder);
        newseq(specorder+1) = p;
        Q.coordinates = [P.coordinates;P.coordinates(row,:) - 2*sdist*geo.nvecs(p,:)];
        Q.parent = [P.parent;P.parent(row)];
        Q.order = [P.order;specorder+1];
        Q.planeseq = [[P.planeseq zeros(size(P.planeseq,1),specorder+1-size(P.planeseq,2))];newseq];
        Q.previous = [P.previous;row];
        Q.amplitude = [P.amplitude;P.amplitude(row)*geo.reflfactors(p)];
        if edseessomething(Q,size(Q.coordinates,1),geo,targets,edgelist,nedgesamples,planesamples,false)
            nnext = nnext + 1;
        end
    end
end

end

%------------------------------------------------------------------------
function yes = edseessomething(P,row,geo,targets,edgelist,nedgesamples,planesamples,useplanes)
% true if the (image) point in row 'row' of P has a valid, unobstructed
% path to at least one of: the target points, a sample point of an edge in
% edgelist, or (if useplanes) a sample point of a plane that it is in
% front of (other than its last reflection plane). The edges are sampled
% at nedgesamples+1 points; the planes at the points in planesamples.

yes = false;
if P.amplitude(row) == 0
    return
end
if ~isempty(targets)
    if any(edchainfactor(P,row,targets,geo) ~= 0)
        yes = true;
        return
    end
end
point = P.coordinates(row,:);
nrefl = P.order(row);
lastplane = 0;
if nrefl > 0
    lastplane = P.planeseq(row,nrefl);
end
if ~isempty(edgelist)
    u = linspace(0,1,nedgesamples+1).';
    u(1) = 1e-6;
    u(end) = 1-1e-6;
    for edgenumber = edgelist(:).'
        p1 = geo.planesatedge(edgenumber,1);
        p2 = geo.planesatedge(edgenumber,2);
        if lastplane ~= 0
            ownplanes = [p1 p2];
            ownplanes = [ownplanes geo.rearsideplane(ownplanes(ownplanes>0)).'];
            if any(ownplanes == lastplane)
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
        if any(edchainfactor(P,row,e1 + u*(e2-e1),geo) ~= 0)
            yes = true;
            return
        end
    end
end
if useplanes
    for p = find(geo.reflfactors ~= 0).'
        if p == lastplane || isempty(planesamples{p})
            continue
        end
        if point*geo.nvecs(p,:).' - geo.planeeqs(p,4) <= geo.geomacc
            continue
        end
        if any(edchainfactor(P,row,planesamples{p},geo) ~= 0)
            yes = true;
            return
        end
    end
end

end

%------------------------------------------------------------------------
function samples = edplanesamples(geo)
% Sample points inside each plane polygon: the corners and the midpoints
% of the polygon sides, moved 1e-3 of the way towards the centroid, the
% centroid, and a 6x6 grid over the bounding rectangle. Only points
% inside the polygon are kept.

nplanes = size(geo.planeeqs,1);
samples = cell(nplanes,1);
ngrid = 6;
for p = 1:nplanes
    nc = geo.ncornersperplanevec(p);
    C = geo.corners(geo.planecorners(p,1:nc),:);
    c0 = mean(C,1);
    M = (C + C([2:nc 1],:))/2;
    X = [c0; C + 1e-3*(c0-C); M + 1e-3*(c0-M)];
    nv = geo.nvecs(p,:);
    e1 = C(2,:) - C(1,:);
    e1 = e1/norm(e1);
    e2 = cross(nv,e1);
    uv = (C - c0)*[e1.' e2.'];
    umin = min(uv(:,1)); umax = max(uv(:,1));
    vmin = min(uv(:,2)); vmax = max(uv(:,2));
    [ug,vg] = meshgrid(umin + ((1:ngrid)-0.5)/ngrid*(umax-umin), vmin + ((1:ngrid)-0.5)/ngrid*(vmax-vmin));
    X = [X; c0 + ug(:)*e1 + vg(:)*e2];
    hitvec = EDpoinpla(X,p*ones(size(X,1),1),geo.minvals,geo.maxvals,...
        geo.planecorners,geo.corners,geo.ncornersperplanevec,geo.nvecs);
    samples{p} = X(hitvec ~= 0,:);
end

end

%------------------------------------------------------------------------
function txt = edwedgetext(geo)
% A text about the concave edges: an edge with the interior (open) angle
% theta < pi allows up to ceil(pi/theta) reflections in a row.

txt = '';
openangle = 2*pi - geo.closwedangvec(:);
iv = find(openangle < pi - 1e-6 & openangle > 1e-6);
if isempty(iv)
    return
end
[minangle,~] = min(openangle(iv));
nmax = ceil(pi/minangle - 1e-9);
txt = [' (The most acute concave edge has the open angle ',num2str(minangle*180/pi,'%.1f'),...
    ' degrees, which allows up to ',int2str(nmax),' reflections between its two planes.)'];

end
