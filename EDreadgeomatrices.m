function [planedata,elapsedtimegeomatrices] = EDreadgeomatrices(corners,planecorners,planerefltypes)
% EDreadgeomatrices - Reads the geometrydata (corners and planecorners) from
% two input matrices, and builds the planedata from these. Rigid surfaces
% are assumed.
%
% Input parameters:
%   corners         A matrix of size [ncorners,3] with the coordinates of
%                   the ncorners corners.
%   planecorners    A matrix of size [nplanes,nmaxnumberofcorners] with the
%                   corner numbers of the nplanes plane definitions.
%   planerefltypes  A vector, [nplanes,1], with values +1,0,-1 specifying
%                   if the planes should be 'RIGID ','TOTABS','SOFT  '
%
% Output parameters:
%	planedata       Structwith fields:
%       corners         Matrix [ncorners,3] with the corner coordinates
%       planecorners    Matrix [planes,nmaxcornersperplane] with the corner numbers that make up each plane.
%                   Since planes can have different numbers of corners, the number of columns is the 
%                   maximum number of corners that any plane has. NB! The corner numbers are in the
%                   renumbered system, not in the CAD-file system.
%       planeabstypes   Matrix [nplanes,nmaxcharacters2] (sparse), with the absorber names 'RIGID ','TOTABS',
%                       'SOFT  ', according to input data planerefltypes.
%       planeeqs        Matrix [nplanes,4] of the plane equations as derived  from the plane definitions. 
%                   Each row has the values [A B C D] of the plane equation on the form Ax + By + Cz = D
%       ncornersperplanevec     Vector [nplanes,1] which gives the number of corners for each plane.
%       minvals         Matrix [nplanes,3] which for each plane gives the smallest coordinate in the x-, y- and z-direction.
%       maxvals         Matrix [nplanes,3] which for each plane gives the largest coordinate in the x-, y- and z-direction.
%       planehasindents Vector [nplanes,1] which for each plane gives the
%                   number of indeting corners
%       indentingcorners  Matrix [nplanes,max(ncornersperplanevec)] which for each plane gives the number to the first corner
%                   in a corner triplet which identifies an indenting corner. The number of the corner is the
%                   order given in planecorners for that plane.
%       cornerinfrontofplane Matrix 
%       modeltype
%  elapsedtimegeomatrices
%
% Uses the functions EDinfrontofplane
% 
% Peter Svensson (peter.svensson@ntnu.no) 10 Oct 2022
%
% [planedata,elapsedtimegeomatrices] = EDreadgeomatrices(corners,planecorners,planerefltypes);

% 29 Nov. 2017 First version
% 12 Jan. 2018 Increased the bounding boxes a bit - doesn't hurt to make
% them a bit bigger.
% 22 Jan 2018 Removed the input parameter planecornerstype since it is not
% used anywhere.
% 12 Apr 2018 Fixed a problem: with a mix of 3-corner planes, and 4-corner
% planes, the last zero of the 3-corner planes need to be replaced by a
% repetition.
% 28 May 2018 Introduced the input parameter planerefltypes
% 13 Aug 2021 Fixed a bug: error messages used "planenumbers" but it was
% corrected to "planecorners" in three locations.
% 10 Oct 2022 Added the timing inside this function and returned as a
% second output parameter.

t00 = clock;

if nargin < 2
    error('ERROR: the corners and planecorners matrices must be specified')
end
if nargin < 3
    nplanes = size(planecorners,1);
    planerefltypes = ones(nplanes,1);
end

% geomacc is only used to make the bounding boxes a bit bigger than the
% corner coordinates. Was 1e-10 earlier.

geomacc = 1e-4;

%---------------------------------------------------------------
% corners is given as an input matrix

ncorners = size(corners,1);

%---------------------------------------------------------------
% planecorners is given as an input matrix

[nplanes,maxcornersperplane] = size(planecorners);
ncornersperplanevec = sum(sign(planecorners),2);

planeabstypes = 'RIGID ';
planeabstypes = planeabstypes(ones(nplanes,1),:);

iv = find(planerefltypes == 0);
if ~isempty(iv)
   for ii = 1:length(iv)
       planeabstypes(iv(ii),:) = 'TOTABS';
   end
end

iv = find(planerefltypes == -1);
if ~isempty(iv)
   for ii = 1:length(iv)
       planeabstypes(iv(ii),:) = 'SOFT  ';
   end
end

if max(max(planecorners)) > ncorners
    error('ERROR: One plane definition in the input matrix planecorners used a higher corner number than was defined in the corners input matrix')
end

% The section below was removed on 22 Jan 2018
% ... and a part was reactivated on 12 Apr 2018. And modified.
%
% % % %---------------------------------------------------------------
% % % % Go through all planes. If there is a plane definition including
% % % % zeros, and planecornerstype == 'circ', expand it repeating the
% % % % same corner order again.
% % % 
% % % if isempty(planecornerstype)
% % % 	planecornerstype = 'circ';
% % % else
% % % 	planecornerstype = char(lower(planecornerstype(1)));
% % % 	if planecornerstype(1) == 'z'
% % % 		planecornerstype = 'zero';
% % %     else
% % % 		planecornerstype = 'circ';
% % % 	end
% % % end
% % % 
% % % if strcmp(planecornerstype,'circ') == 1
	for ii = 1:nplanes
		if ncornersperplanevec(ii) ~= maxcornersperplane
            planecorners(ii,ncornersperplanevec(ii)+1) = planecorners(ii,1);
%             iv = [1:ncornersperplanevec(ii)];		
%             pattern = planecorners(ii,iv);
% 			nrepeatings = ceil(ncornersperplanevec(ii)/ncornersatplane);
% 			for jj = 1:nrepeatings-1
% 				pattern = [pattern planecorners(ii,iv)];
% 			end
% 			planecorners(ii,:) = pattern(1:ncornersperplane);
		end
	end
% % % end

%---------------------------------------------------------------
% Find the normal vectors to the planes using the cross products
%
% The normal vector is basically the cross product between two vectors
% connecting three consecutive points. If there are indents though
% they will cause reversed normal vectors, so one must go through all
% combinations and choose the majority normal vector.
%
% 26mar09  Use the fact described above for detecting indention corners

planenvecs = zeros(nplanes,3);
planehasindents = zeros(nplanes,1);
indentingcorners = sparse(zeros(nplanes,max(ncornersperplanevec)));

for ii = 1:nplanes
    
	iv = planecorners(ii,:)~=0;
	cornerlist = planecorners(ii,iv);
	iv = find(cornerlist == cornerlist(1));
	if length(iv) > 1
		cornerlist = cornerlist(1:iv(2)-1);
	end
	ncorners = length( cornerlist );
	cornerlist = [cornerlist cornerlist(1) cornerlist(2)];

	nvectorlist = zeros(ncorners,3);
	nveclen = zeros(ncorners,1);	
	
	for jj = 1:ncorners
		co1numb = cornerlist(jj);
		co2numb = cornerlist(jj+1);
		co3numb = cornerlist(jj+2);
		vec1 = (corners(co1numb,:) - corners(co2numb,:)).';
		vec2 = (corners(co3numb,:) - corners(co2numb,:)).';
% 		nvec = EDcross(vec1.',vec2.').';
        nvec = [vec1(2,:).*vec2(3,:)-vec1(3,:).*vec2(2,:)
            vec1(3,:).*vec2(1,:)-vec1(1,:).*vec2(3,:)
            vec1(1,:).*vec2(2,:)-vec1(2,:).*vec2(1,:)];
        nvec = reshape(nvec,size(vec1));
		nveclen(jj) = norm(nvec);
		if nveclen(jj) > 0
			nvectorlist(jj,:) = nvec./nveclen(jj);
		end
	end
	
	iv = nveclen < max(nveclen)*0.001;
	nvectorlist(iv,:) = [];
	nvecref = nvectorlist(1,:);
	
    [n1,~] = size(nvectorlist);
    
	nvecsigns = round(sum(     (nvectorlist.').*(nvecref(ones(n1,1),:).')      ));
    
    if sum(nvecsigns) == 0
        disp(' ')        
       error(['ERROR: Plane ',int2str(planecorners(ii)),' (plane numbering as in the CAD file) seems to be twisted.'])        
    end
    
    if abs(sum(nvecsigns)) ~= n1
       nindents = (n1 - abs(sum(nvecsigns)))/2;
       disp(['Plane ',int2str(ii),' has ',int2str(nindents),' indents!']) 
       planehasindents(ii) = nindents;
       if sum(nvecsigns) > 0
           ivindent = find(nvecsigns == -1);
       else
            ivindent = find(nvecsigns == 1);           
       end
       if length(ivindent) == nindents
            indentingcorners(ii,ivindent) = 1;
       else
           if length(ivindent) == ncorners - nindents
               ivindent = nvecsigns == 1;
                indentingcorners(ii,ivindent) = 1;           
           else
              error(['ERROR: An unexpected problem. Please report to the developer'])
           end
       end
    end
    
    nvecdiff = [nvectorlist(2:n1,1).*nvecsigns(2:n1).' nvectorlist(2:n1,2).*nvecsigns(2:n1).' nvectorlist(2:n1,3).*nvecsigns(2:n1).'] - nvecref(ones(n1-1,1),:);
        
    if n1 > 2
        nvecdiff = sum(nvecdiff.'.^2).';    
    else
        nvecdiff = norm(nvecdiff);
    end
    
    if any(nvecdiff>1e-4)
        nvecdiff
        error(['ERROR: Normal vectors for plane ',int2str(planecorners(ii)),' (in the CAD file, = ',int2str(ii),' in the ESIE2 file), get different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    elseif any(nvecdiff>1e-8)
        nvecdiff
        disp(['WARNING: Normal vectors for plane ',int2str(planecorners(ii)),' (in the CAD file, = ',int2str(ii),' in the ESIE2 file), get somewhat different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    end
    
	if ncorners > 5 && abs(sum(nvecsigns)) <= 1
		disp(['WARNING for plane number ',int2str(planecorners(ii)),' in the CAD-file'])
		disp(['   with the name ',strtrim(char(full(planenames(ii,:))))])
		disp('   The normal vector can not be determined for this plane because there are')
		disp('   the same number of inwards and outwards corners')
		disp('   This is a list of the possible normal vectors:')
		[nv1,~] = size(nvectorlist);
		for kk = 1:nv1
			vecstr = ['   ',int2str(kk),'. ',num2str(-nvectorlist(kk,1)),' ',num2str(-nvectorlist(kk,2)),' ',num2str(-nvectorlist(kk,3))];
			disp(vecstr)
		end
		disp(' ')
	
      preferredsign = input('   Give the number of a correct normal vector for this plane please ');
      switchsign = nvecsigns.'./nvecsigns(preferredsign);
		nvectorlist = nvectorlist.*switchsign(:,ones(1,3));	
    else
		mostcommonsign = sign(sum(nvecsigns));

		switchsign = nvecsigns.'./mostcommonsign;
		nvectorlist = nvectorlist.*switchsign(:,ones(1,3));
    end

	planenvecs(ii,:) = mean(nvectorlist);
end

planenvecs = -planenvecs;

%---------------------------------------------------------------
% Plane equations, as Ax + By + Cz = D for each plane

planeeqs = zeros(nplanes,4);
planeeqs(:,1:3) = planenvecs;
planeeqs(:,4) = sum( (planenvecs.').*(corners(planecorners(:,1),:).')  ).';

%---------------------------------------------------------------
% Useful data: planesatcorners, minvals and maxvals

[ncorners,~] = size(corners);
planesatcornerhits = zeros(ncorners,nplanes);

for ii = 1:nplanes
	cornerlist = planecorners(ii,1:ncornersperplanevec(ii));
	planesatcornerhits(cornerlist,ii) = planesatcornerhits(cornerlist,ii) + 1;
end

maxplanespercorner = 0;
for ii = 1:ncorners
	nplanes = length(find(planesatcornerhits(ii,:) ~= 0));
	if nplanes > maxplanespercorner
		maxplanespercorner = nplanes;
	end	
end

planesatcorners = zeros(ncorners,maxplanespercorner);
nplanespercorners = zeros(ncorners,1);
for ii = 1:ncorners
	iv = find(planesatcornerhits(ii,:)~=0);
	planesatcorners(ii,1:length(iv)) = iv;
	nplanespercorners(ii) = length(iv);
end

% find cubic boxes inside which the planes are placed

[nplanes,~] = size(planeeqs);

minvals = zeros(nplanes,3);
maxvals = zeros(nplanes,3);

for ii = 1:nplanes
	cornerlist = planecorners(ii,:);
	cornerlist = cornerlist( cornerlist~= 0 );
	cornercoord = corners(cornerlist,:);
	minvals(ii,:) = min(cornercoord);
	maxvals(ii,:) = max(cornercoord);
end

minvals = minvals - geomacc;
maxvals = maxvals + geomacc;

%---------------------------------------------------------------
%
%			cornerinfrontofplane
%
%---------------------------------------------------------------
%
% Generate a matrix which can be used to find out if the model is convex
%  + be useful later:
% cornerinfrontofplane, size [nplanes,ncorners]
% Values are:
%   1 means that a point is in front of the plane
%   0 means that a point is aligned with a plane
%   -1 means that a point is behind a plane
% All corners that belong to a plane will have the value 0.

% Corner number is given by the col no.
if ncorners < 256
    cornernumb = uint8([1:ncorners]);
elseif ncorners < 65536
    cornernumb = uint16([1:ncorners]);
else
    cornernumb = uint32([1:ncorners]);        
end

% Plane numbers is given by the row no.
if nplanes < 256
    planenumb = uint8([1:nplanes].');
elseif nplanes < 65536
    planenumb = uint16([1:nplanes].');
else
    planenumb = uint32([1:nplanes].');
end

cornerinfrontofplane = EDinfrontofplane(corners,planeeqs(:,1:3),planecorners,[],cornernumb,planenumb);
cornerinfrontofplane = reshape(cornerinfrontofplane,nplanes,ncorners);

if any(any(cornerinfrontofplane==1)) == 0
   modeltype = 'convex_ext'; 
elseif any(any(cornerinfrontofplane==-1)) == 0
   modeltype = 'convex_int'; 
else
    modeltype = 'other';
end
if any(any(abs(cornerinfrontofplane))) == 0
   if ncorners == max(ncornersperplanevec)
      modeltype = 'singleplate'; 
   else
       modeltype = 'thinplates';       
   end
end

%---------------------------------------------------------------
% Store the relevant variables in the output struct

% NB! Sparse matrices can not get a non-double format

if ncorners < 256
    planecorners = uint8(planecorners);
elseif ncorners < 65536
    planecorners = uint16(planecorners);    
end   

if max(ncornersperplanevec) <= 255
    ncornersperplanevec = uint8(ncornersperplanevec);
    planehasindents = uint8(planehasindents);
else
    ncornersperplanevec = uint16(ncornersperplanevec);
    planehasindents = uint16(planehasindents);
end

planeabstypes = sparse(planeabstypes+1-1);

planedata = struct('corners',corners,'planecorners',planecorners,...
    'planeabstypes',planeabstypes,...
    'planeeqs',planeeqs,'ncornersperplanevec',ncornersperplanevec,...
    'minvals',minvals,'maxvals',maxvals,...
    'planehasindents',planehasindents,'indentingcorners',indentingcorners,...
    'cornerinfrontofplane',cornerinfrontofplane,'modeltype',modeltype);

elapsedtimegeomatrices = etime(clock,t00);
