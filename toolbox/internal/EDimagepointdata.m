function [Sdata_ie,Rdata_ie] = EDimagepointdata(edgedata,Sdata,Rdata,firstorderpathdata)
% EDimagepointdata - Adds image sources and image receivers to Sdata and Rdata, for the integral equation.
% For a non-convex object (geoinputdata.allownonconvex = 1), the edges can
% be reached via specular reflections: from the source to the first edge
% (image sources) and from the last edge to the receiver (image
% receivers). EDfindGApaths finds these image points and the parts of each
% edge that they can see. This function turns them into additional
% "sources" and "receivers" for the integral equation:
%   - every image source becomes an extra source, whose contribution to the
%     source term is added to that of its parent source,
%   - every image receiver becomes an extra receiver, whose propagated
%     sound pressure is added to that of its parent receiver.
% Each column ("point") sees each edge over a single interval, given by
% the fields .vispartedgesfroms_start/_end (and ..r_start/_end), as
% fractions 0-1 of the edge length. A point that sees an edge over several
% separate intervals is represented by several columns.
% The sources and receivers themselves come first, in their original
% order, and their visible parts of the edges are also replaced by those
% from EDfindGApaths.
%
% Input parameters:
%   edgedata, Sdata, Rdata, firstorderpathdata   see help EDdatastructs.
%       firstorderpathdata must come from EDfindGApaths.
%
% Output parameters:
%   Sdata_ie, Rdata_ie  The structs Sdata and Rdata with these fields
%                   changed or added (for Sdata: ...s/...S; for Rdata
%                   ...r/...R):
%       .coordinates            [npoints,3] the points (sources/receivers
%                               first, then the extra columns)
%       .parentsource           [npoints,1] (Sdata_ie) or .parentreceiver
%                               (Rdata_ie): the source/receiver that each
%                               point belongs to
%       .vispartedgesfroms      [nedges,npoints] fully visible value
%                               (2^nedgesubs-1) where the point sees the
%                               edge (partly), otherwise 0
%       .vispartedgesfroms_start, .vispartedgesfroms_end  [nedges,npoints]
%                               the visible part, fractions 0-1 (double)
%       .reftoshortlistS, .rSsho, .thetaSsho, .zSsho  edge-related
%                               coordinates of all the points
%       .sourceamplitudes       [npoints,nfrequencies] (Sdata_ie) the
%                               parent's amplitude times the product of the
%                               reflection factors
%       .EDinputdatahash        a new hash, which depends on the original
%                               hash and that of firstorderpathdata
%
% See also EDfindGApaths, EDintegralequation_convex_tf
%
% Peter Svensson (peter.svensson@ntnu.no) 20 Nov. 2026
%
% [Sdata_ie,Rdata_ie] = EDimagepointdata(edgedata,Sdata,Rdata,firstorderpathdata);

% 20 Nov. 2026 First version (by Claude, prompted by PS).

if ~isfield(firstorderpathdata,'Spoints')
    error('ERROR: EDimagepointdata needs a firstorderpathdata struct from EDfindGApaths')
end

nedges = size(edgedata.edgecorners,1);
maxvisibilityval = 2^Sdata.nedgesubs - 1;

Sdata_ie = edaddpoints(Sdata,firstorderpathdata.Spoints,edgedata,nedges,maxvisibilityval,'s');
Rdata_ie = edaddpoints(Rdata,firstorderpathdata.Rpoints,edgedata,nedges,maxvisibilityval,'r');

nfreqamp = size(Sdata.sourceamplitudes,2);
Sdata_ie.sourceamplitudes = Sdata.sourceamplitudes(Sdata_ie.parentsource,:).*...
    repmat(Sdata_ie.pointamplitude,1,nfreqamp);
Sdata_ie = rmfield(Sdata_ie,'pointamplitude');
% The image receivers' reflection factors are applied when their results
% are added to their parent receivers
Rdata_ie.receiveramplitude = Rdata_ie.pointamplitude;
Rdata_ie = rmfield(Rdata_ie,'pointamplitude');

allversions = EDstageversions;
thisversion = allversions.versions.EDimagepointdata;
Sdata_ie.EDinputdatahash = DataHash(struct('parent',EDparenthash(Sdata),...
    'paths',firstorderpathdata.EDinputdatahash,'type','Sdata_ie','version',thisversion));
Rdata_ie.EDinputdatahash = DataHash(struct('parent',EDparenthash(Rdata),...
    'paths',firstorderpathdata.EDinputdatahash,'type','Rdata_ie','version',thisversion));

end

%------------------------------------------------------------------------
function D = edaddpoints(D,P,edgedata,nedges,maxvisibilityval,typeofcoords)

nbase = size(D.coordinates,1);
intervals = P.edgeintervals;

% Columns: first one per base point, then the extra ones. For each point
% (row of P), the k-th column takes the k-th interval on each edge.
colpoint = [];
colnumber = [];
nrows = size(P.coordinates,1);
for row = 1:nrows
    nint = 0;
    if ~isempty(intervals)
        ivrow = intervals(:,1)==row;
        if any(ivrow)
            nint = max(accumarray(intervals(ivrow,2),1));
        end
    end
    if row <= nbase
        ncols = max(nint,1);
    else
        ncols = nint;
    end
    colpoint = [colpoint;row*ones(ncols,1)];
    colnumber = [colnumber;(1:ncols).'];
end
% The base points first, in their order
[~,isort] = sortrows([colpoint>nbase colnumber>1 colpoint colnumber]);
colpoint = colpoint(isort);
colnumber = colnumber(isort);
npoints = numel(colpoint);

vispart = zeros(nedges,npoints);
visstart = zeros(nedges,npoints);
visend = zeros(nedges,npoints);
for ii = 1:size(intervals,1)
    row = intervals(ii,1);
    edgenumber = intervals(ii,2);
    % Which interval number is this, for this point and edge?
    k = sum(intervals(1:ii,1)==row & intervals(1:ii,2)==edgenumber);
    col = find(colpoint==row & colnumber==k);
    vispart(edgenumber,col) = maxvisibilityval;
    visstart(edgenumber,col) = intervals(ii,3);
    visend(edgenumber,col) = intervals(ii,4);
end

coordinates = P.coordinates(colpoint,:);

% Edge-related coordinates, as in EDSorRgeo
rcomplete = zeros(nedges,npoints);
thetacomplete = zeros(nedges,npoints);
zcomplete = zeros(nedges,npoints);
for ii = 1:nedges
    iv = find(vispart(ii,:) > 0);
    if ~isempty(iv)
        [r,theta,z] = EDcoordtrans1(coordinates(iv,:),...
            [edgedata.edgestartcoords(ii,:);edgedata.edgeendcoords(ii,:)],...
            edgedata.edgenvecs(ii,:),reshape(edgedata.edgerelatedcoordsysmatrices(ii,:),3,3));
        rcomplete(ii,iv) = r(:).';
        thetacomplete(ii,iv) = theta(:).';
        zcomplete(ii,iv) = z(:).';
    end
end
[reftoshortlist,rsho,thetasho,zsho,~] = EDcompress3(rcomplete,thetacomplete,zcomplete);

D.coordinates = coordinates;
D.pointamplitude = P.amplitude(colpoint);
% The plane visibility fields keep their meaning for the original points;
% the extra columns get zeros (they are not used by the integral equation)
if typeofcoords == 's' && isfield(D,'visplanesfroms')
    D.visplanesfroms = [D.visplanesfroms zeros(size(D.visplanesfroms,1),npoints-nbase)];
elseif typeofcoords == 'r' && isfield(D,'visplanesfromr')
    D.visplanesfromr = [D.visplanesfromr zeros(size(D.visplanesfromr,1),npoints-nbase)];
end
if typeofcoords == 's'
    D.parentsource = P.parent(colpoint);
    D.vispartedgesfroms = vispart;
    D.vispartedgesfroms_start = visstart;
    D.vispartedgesfroms_end = visend;
    D.reftoshortlistS = reftoshortlist;
    D.rSsho = rsho;
    D.thetaSsho = thetasho;
    D.zSsho = zsho;
else
    D.parentreceiver = P.parent(colpoint);
    D.vispartedgesfromr = vispart;
    D.vispartedgesfromr_start = visstart;
    D.vispartedgesfromr_end = visend;
    D.reftoshortlistR = reftoshortlist;
    D.rRsho = rsho;
    D.thetaRsho = thetasho;
    D.zRsho = zsho;
end

end
