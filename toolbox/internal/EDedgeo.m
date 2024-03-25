function [edgedata,planedata,elapsedtimeedgeo,existingfilename] = EDedgeo(planedata,...
    geoinputdata,EDversionnumber,filehandlingparameters)
% EDedgeo - Calculates some plane- and edge-related geometrical
% parameters, based on corners and planes in the planedata struct.
%
% Input parameters:
%	planedata               A struct with the corners and plane data.
%   geoinputdata            A struct where these fields are optionally used
%	.firstcornertoskip 	    All edges including at least one corner with this number or
%                   		higher will be excluded from the calculations. Default: 1e6
%   .listofcornerstoskip    All edges including at least one corner
%                           in this list will be excluded from the
%                           calculations. Default: [].
%   .listofedgestoskip      This is a vertical matrix with one or two
%                           columns. If it has one column, then edges with
%                           these numbers will be turned off. If the matrix
%                           has two columns, then edges with the corner
%                           pairs will be turned off. .listofedgestoskip
%                           takes precedence over .firstcornertoskip and
%                           .listofcornerstoskip. Default: [].
%   .planeseesplanestrategy If this parameter is given the value 1, a plane-to-plane
%                           visibility check is done by checking the plane-midpoint to plane-midpoint
%                           for obstructions, which might be enough for
%                           some geometries. Default: 0.
%  EDversionnumber 
%  filehandlingparameters   a struct which contains the field showtext.
% 
% Output parameters:
%	edgedata				Struct with fields as below
%
%       edgecorners             Matrix [nedges,2] with the corner numbers that
%	                        define each edge, in the reference order.
%       edgestartcoords         Matrix [nedges,3] with the coordinates of the
%                           start corner for each edge.
%       edgeendcoords           Matrix [nedges,3] with the coordinates of the
%                           end corner for each edge.
%       edgelengthvec           List [nedges,1] with the lengths of each edge.
%       planesatedge            Matrix [nedges,2] with the plane numbers of the
%                           two planes that connect to each edge. The first
%                           plane is considered as the reference plane for
%                           each edge. Rule: if the RH thumb is placed
%                           along the edge, from start to end, and the hand
%                           is rotated, the fingers should "come up
%                           through" the reference plane.
%       closwedangvec           List [nedges,1] with the angles in radians of
%                           each edge - for the closed part of each edge.
%       edgenvecs               Matrix [nedges,3] with the normal vectors of
%                           the reference plane for each edge.
%       offedges                List [nlength,1] where nlength = 0-nedges,
%                           containing the edge numbers that should not
%                           used.
%       reflfactors             List [nplanes,1] with the reflection factors
%                           of each plane. The only supported values are
%                           +1 for rigid planes, 0 for perfectlt absorbing
%                           planes and -1 for pressure-release planes. The
%                           pressure-release cases might not be implemented
%                           fully.
%       planeisthin             List [nplanes,1] with the value 1 or 0 stating
%                           whether the plane is thin or not. 
%       rearsideplane           List [nplanes,1] with the plane number that is at
%                           the rear side. If a plane is non-thin, the
%                           value in this list is zero.
%       planehasindents         List [nplanes,1] with the value 1 or 0 stating
%                           whether a plane has indents or not.
%       indentingedgepairs      Matrix [nindentingedgepairs,2] where each row
%                           gives two connected edges that form an indent
%                           in the shared plane. Consequently, those two
%                           edges can never see each other.
%       canplaneobstruct        List [nplanes,1] with the value 1 or 0 stating
%                           whether a plane has the potential to obstruct
%                           or not.
%       planeseesplane          A matrix, [nplanes,nplanes], with the value 1,
%                           0,-1,-2 stating:
%                           1 A plane is front of another plane and could
%                             then potentially see it.
%                           0 A plane is completely behind another plane
%                             but they are not aligned
%                           -1 A plane is aligned with another plane
%                           -2 A plane is back-to-back with another plane.
%       edgesatplane            A matrix, [nplanes,nmaxedgesperplane] which for
%                           row N contains the edge numbers that connect to
%                           to plane N.
%       edgeseesplane           A matrix, [nplanes,nedges], which contains the
%                           values 0,1,-2 or -2, meaning:
%                           0   The edge is completely behind a plane
%                           -1  The edge can never see the plane and the
%                               edge is aligned with the plane
%                           -2  The edge can never see the plane and the
%                               edge belongs to the plane
%                           1   The edge has at least one point in front of the plane.
%       edgerelatedcoordsysmatrices    A matrix, [nedges,9], containing, for each edge, the 3*3
%                           matrices needed for finding the edge-related
%                           coordinate system. Each row has to be reshaped
%                           by
%                           Bmatrix =
%                           reshape(edgerelatedcoordsysmatrices(edgenumber,:),3,3);    
%   planedata               Struct with some added fields.
%       planeisthin 
%       planeseesplane
%       rearsideplane
%       canplaneobstruct
%       reflfactors
%   elapsedtimeedgeo        This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimeedgeo has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        For v2 of this function, if an existing file was
%                           found that was reused, then the reused file
%                           name is given here. If no existing file could be 
%                           reused then this variable is empty. For v1 of
%                           this function, this variable is also empty.
% 
% Uses the functions EDcoordtrans1, EDinfrontofplane , EDrecycleresultfiles
% from EDtoolbox
% Uses the function Datahash from Matlab Central
% 
% Peter Svensson (peter.svensson@ntnu.no) 25 March 2024
%
% [edgedata,planedata,elapsedtimeedgeo,existingfilename] = EDedgeo...
% (planedata,geoinputdata,EDversionnumber,filehandlingparameters);

% 9 July 2009   Stable version
% 28 Sep. 2014  Added the optional input parameter listofcornerstoskip
% 1 Dec. 2014   Added the new output parameter edgerelatedcoordsysmatrices
% 3 Dec. 2014   Changed the input from cageofile to the struct planedata.
%               Also changed the output to a struct edgedata + an extension
%               of the struct planedata
% 31 March 2015 Added '' in the file saving, to handle file and directory
%               names with blanks
% 29 Dec. 2016  Added two fields to the struct edgedata: ismodelconvex and
%               ismodelthinplate
% 12 Jan. 2017  Corrected one missing structform of planedata.corners.
%               Also moved def. of planenvecs to outside an if.
% 14 Nov. 2017  A small mistake was fixed with the firstcornertoskip
% 26 Nov. 2017  Trimmed down version from ESIE2edgeox: removed the input
%               parameters open_or_closed_model, int_or_ext_model,
%               specorder, nedgesubs. Removed the file saving. 
% 27 Nov. 2017  Removed the difforder parameter, because quite little computation is
%               saved. Changed name for one parameter to firstcornertoskip
%               (from firstskipcorner).
% 28 Nov. 2017  Introduced the non-global showtext parameter. Made some
%               code improvements.
% 29 Nov. 2017  Small correction with findstr( and modeltype. Changed call,
%               from ESIE2checkobstrpaths to ED
% 22 Jan 2018 The output parameters ismodelconvex and ismodelthinplate were
%               removed (weren't used other places).
% 23 Jan 2018 Error from yesterday's change: a long line with ... at the
%             end can not be followed by a commented out line!
% 8 Feb 2018 Introduced the EDinputdatahash
% 15 Mar 2018 Fixed a but that occurred when some plane was TOTABS (in a
%             cad-file) or SOFT
% 15 Mar 2018 Fixed a bug which occured if the planecorners matrix had some
%             zeros; a definition of zeros used type int8 instead of 'int8'
% 28 May 2018 Fixed a small bug: planeabstypes = TOTABS or SOFT were not
%             recognized.
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 28 Oct. 2023 Changed the set of input parameters somewhat (removed the
% info about "version 1".
% 30 Oct. 2023 Moved the EDversionnumber position among input parameters.
% Fine-tuned the EDinputdatahash
% 25 March 2024 Implemented the new input parameter listofedgestoskip

t00 = clock;
geomacc = 1e-10;

showtext = filehandlingparameters.showtext;
if isempty(geoinputdata.firstcornertoskip)
    firstcornertoskip = 1e6;
else
    firstcornertoskip = geoinputdata.firstcornertoskip;
end
if isempty(geoinputdata.listofcornerstoskip)
    listofcornerstoskip = [];
else
    listofcornerstoskip = geoinputdata.listofcornerstoskip;
end
if isempty(geoinputdata.listofedgestoskip)
    listofedgestoskip = [];
else
    listofedgestoskip = geoinputdata.listofedgestoskip;
end
if isempty(geoinputdata.planeseesplanestrategy)
    planeseesplanestrategy = 0;
else
    planeseesplanestrategy = geoinputdata.planeseesplanestrategy;
end

EDinputdatastruct = struct('corners',planedata.corners,'planecorners',planedata.planecorners,...
    'planeabstypes',planedata.planeabstypes,...
'firstcornertoskip',firstcornertoskip,'listofcornerstoskip',listofcornerstoskip,...
'listofedgestoskip',listofedgestoskip,...
'planeseesplanestrategy',planeseesplanestrategy,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. Should the data be saved in a file? 

if filehandlingparameters.suppressresultrecycling == 1
	foundmatch = 0;
	existingfilename = '';
else
	[foundmatch,existingfilename] = ... 
		EDrecycleresultfiles(filehandlingparameters.outputdirectory,...
		'_eddata',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_eddata.mat'];

if foundmatch == 1
%		eval(['load ''',existingfilename,''''])
	eval(['load(''',existingfilename,''')'])
	if ~strcmp(existingfilename,desiredname)
		copyfile(existingfilename,desiredname);
	end
	elapsedtimeedgeo_new = etime(clock,t00);
	elapsedtimeedgeo = [elapsedtimeedgeo_new elapsedtimeedgeo];
	return
end

%----------------------------------------------------------------
% No existing file can be used

ncorners = size(planedata.corners,1);

[nplanes,maxncornersperplane] = size(planedata.planecorners);
minncornersperplane = min(planedata.ncornersperplanevec);
nedgesperplanevec = double(planedata.ncornersperplanevec);
ncornersperplanevec = double(planedata.ncornersperplanevec);

% 14.11 2017: it was: if .... firstcornertoskip < ncorners.
% Should be <= ncorners

if (firstcornertoskip < 1e6) && (firstcornertoskip <= ncorners)
%     switchoffcorners = find( cornernumbers >= firstcornertoskip );
    switchoffcorners = firstcornertoskip:ncorners;
else
   switchoffcorners =  listofcornerstoskip;
end

% onesvec2 = ones(nplanes,1);
zerosvec1 = zeros(nplanes,1);    
    
zerosvec2 = zeros(nplanes,nplanes,'int8');

%----------------------------------------------------------------------------
%
%		GEOMETRICAL ACOUSTICS PARAMETERS
%
%----------------------------------------------------------------------------
% 
% Go through all planeabstypes. If it is:
%	'SOFT' or 'soft'  then reflfactors should be -1.
%	'TOTA' or 'tota' then reflfactors should be 0 (short for TOTABS)
%   Otherwise it is assumed to be rigid, so reflfactors = 1.
% 
% 28 May 2018 Fixed a bug: planeasbtypes = lower(... didn't work without
% the lower(setstr(...

nchars = size(planedata.planeabstypes,2);
reflfactors = ones(nplanes,1);
if nchars >= 4
    planeabstypes = lower(setstr(full(planedata.planeabstypes(:,1:min([4,nchars])))));
    
    comptxt = 'soft';
    ivpotential = find(planeabstypes(:,1)==comptxt(1));
    if ~isempty(ivpotential)
        comptxt = 'oft';
        compmat = comptxt(ones(length(ivpotential),1),:);
        % Fixed bug 15 Mar 2018: wrong use of logical indexing. Need the
        % "find"
%         ivsoft   = ivpotential(prod( (planeabstypes(ivpotential,2:4).*compmat).' ).');
        ivsoft   = ivpotential(find(prod( (planeabstypes(ivpotential,2:4).*compmat).' ).'));
        reflfactors(ivsoft)   = -1*ones(size(ivsoft));
    end
    
    comptxt = 'tota';
    ivpotential = find(planeabstypes(:,1)==comptxt(1));
    if ~isempty(ivpotential)
        comptxt = 'ota';
        compmat = comptxt(ones(length(ivpotential),1),:);
        % Fixed bug 15 Mar 2018: wrong use of logical indexing. Need the
        % "find"
%         ivtotabs =  ivpotential(prod( double(planeabstypes(ivpotential,2:4)==compmat).' ).');
        ivtotabs =  ivpotential(find(prod( double(planeabstypes(ivpotential,2:4)==compmat).' ).'));
        reflfactors(ivtotabs) = zeros(size(ivtotabs));
    end
end

%----------------------------------------------------------------------------
%
%		EDGE PARAMETERS
%
%----------------------------------------------------------------------------
%
% Check that the planecorners matrix contains no zeros. In such case,
% add the circular repetition of coordinates.

planecorners = planedata.planecorners;

if sum(sum( planecorners == 0 )) ~= 0
	for ii = 1:nplanes
		iv = find( planecorners(ii,:) ~= 0);
		ncornersatplane = length(iv);
		if ncornersatplane ~= maxncornersperplane
			pattern = planecorners(ii,iv);
			nrepeatings = ceil(maxncornersperplane/ncornersatplane);
			for jj = 1:nrepeatings-1
				pattern = [pattern planecorners(ii,iv)];
			end
			planecorners(ii,:) = pattern(1:maxncornersperplane);
		end
	end
end
    
%--------------------------------------------------------------------------------
%
% Find all edges
%
% Go through all planes. All consecutive pairs of corners form an edge.
% The list planesatedge gives the one or two planes that connect
% to each edge. If the second plane number, for a certain edge,
% is zero, then that specific edge is a plate edge.

% edgesatplane = zeros( size(planecorners) );

nedgesguess = sum(ncornersperplanevec);
edgecorners = zeros(nedgesguess,2);
tempplanesatedge = zeros(nedgesguess,1);
% nedges = 0;
% thirdpoint = [];
edgecounter = 0;
multfac = round(10^(floor(log10(nedgesguess))+2));

% First go through all planes, and construct edges from the list of corners
% that define every plane. This way, many edges will occur twice, unless
% a thin plane has no rear side (then the edge occurs a single time)
% or if two thin planes are connected (then the edge occurs four times)

for ii = 1:nplanes

   for jj = 1:nedgesperplanevec(ii)
      edgecounter = edgecounter + 1;   
      corner1 = planecorners(ii,jj);
      if jj + 1 <= ncornersperplanevec(ii)
         corner2 = planecorners(ii,jj+1);
      else
         corner2 = planecorners(ii,1);
      end
	  edgecorners(edgecounter,:) = [corner1 corner2];
	  tempplanesatedge(edgecounter) = ii;
   end

end

if edgecounter < nedgesguess
	edgecorners(edgecounter+1:nedgesguess,:) = [];
	tempplanesatedge(edgecounter+1:nedgesguess,:) = [];
end

% To find all duplicates, the edge list is sorted in numerical order
% If an edge is defined twice or four times, then this is an OK edge.
% If an edge is defined a single time, it is connected to a plane which is not connected, and then the
% edge should be switched off.
% If an edge is defined three times, there must be something like a reflector array where one of
% the rear planes is missing.

[test,flipvec] = sort(edgecorners.');
test = test.';
flipvec = flipvec(1,:).';

test = test(:,1)*multfac + test(:,2);
[test,sortvec] = sort(test);
tempplanesatedge = tempplanesatedge(sortvec);
flipvec = flipvec(sortvec);
planesatedge = [tempplanesatedge.*(flipvec==1) tempplanesatedge.*(flipvec==2)];

% Check if there are loners; remove and give a warning.
 
dtest = diff([0;test;test(length(test))+100]);
ntest = length(dtest);
loners = find(dtest(1:ntest-1)~=0 & dtest(2:ntest)~=0);

if ~isempty(loners)
	disp('WARNING! Some edges had only one plane connected and they will be excluded.')
	disp('         Check if they should be corrected. ')
	disp('         Remember that reflectors must have a plane at the rear side too')
	
	test(loners) = [];
	planesatedge(loners,:) = [];
	
end

ntest = length(test);
if ntest >= 2
    if test(ntest) ~= test(ntest-1)
	    test(ntest) = [];
	    planesatedge(ntest,:) = [];
    end
end

% Check if there are triplets

triptest = test(1:2:length(test))-test(2:2:length(test));
triplets = find(triptest~=0);

if ~isempty(triplets)
	disp('ERROR: Triplets: some edges are defined three times which means that two thin planes have a correct')
	disp('       definition but some error on the rear side. You must correct these. ')
	disp('       Only the first one can be detected; it is:')
	disp('       (ED corner numbers for the two edge points are given)')
	[floor(test(triplets(1)*2-1)/multfac) test(triplets(1)*2-1)-floor(test(triplets(1)*2-1)/multfac)*multfac]
	pause
end

% The planesatedge matrix is twice as long as it should so it has to be
% compiled. If the number of rows with zeros in the first column is
% different from the number of rows w zeros in the second column, then the
% cad model is wrong: most probably some plane has the wrong order of
% corners.

edgecounter = length(test);

iv1 = find(planesatedge(:,1)~=0);
iv2 = 1:edgecounter;
iv2(iv1) = [];

if length(iv2) ~= length(iv1)
   disp('ERROR: The cad-file seems to have some mistake. Please check')
   disp('       that all the planes have the same order of the corners')
   error
end

edgecorners = [floor(test(iv1)/multfac) test(iv1)-floor(test(iv1)/multfac)*multfac];

planesatedge = planesatedge(iv1,:) + planesatedge(iv2,:);

nedges = size(edgecorners,1);

zerosvec3 = zeros(nedges,1);
zerosvec5 = zeros(nedges,3);

thirdpoint = zerosvec5;
for ii = 1:nedges
% 	refplane = planesatedge(ii,1);
	secondplane = planesatedge(ii,2);
	if secondplane ~= 0
		edgeco1 = planedata.corners(edgecorners(ii,1),:);
		edgeco2 = planedata.corners(edgecorners(ii,2),:);
		edgemidpoint = edgeco1 + (edgeco2-edgeco1)/2;
%		intoplanevec = cross( edgeco2-edgeco1,planenvecs(secondplane,:));
% 		intoplanevec = ESIE2cross( (edgeco2-edgeco1).',(planedata.planeeqs(secondplane,1:3)).').';
        vec1 = (edgeco2-edgeco1).';
        vec2 = (planedata.planeeqs(secondplane,1:3)).';
        intoplanevec = [vec1(2,:).*vec2(3,:)-vec1(3,:).*vec2(2,:)
            vec1(3,:).*vec2(1,:)-vec1(1,:).*vec2(3,:)
            vec1(1,:).*vec2(2,:)-vec1(2,:).*vec2(1,:)];
        intoplanevec = reshape(intoplanevec,size(vec1)).';
        
		inplanepoint = edgemidpoint + intoplanevec*0.1;
		if sum(abs(inplanepoint)) == 0
			inplanepoint = edgemidpoint + intoplanevec*0.2;		
		end
		thirdpoint(ii,:) = inplanepoint;
	end
end

%---------------------------------------------------------------
% Calculate the closwedang values for all edges
%
% 1 Dec 2014 Added storing the B-matrices for edge-related coordinate
%            transformations

closwedangvec = zerosvec3;

ivec = find( planesatedge(:,1).*planesatedge(:,2) ~= 0);
xwedgevec = [planedata.corners(edgecorners(ivec,1),:) planedata.corners(edgecorners(ivec,2),:)];
nvec1vec  = planedata.planeeqs( planesatedge(ivec,1),1:3);
xsouvec   = thirdpoint(ivec,:);

for jj = 1:length(ivec)
  
   ii = ivec(jj);
   [~,thetas,~,~] = EDcoordtrans1(xsouvec(jj,:),[xwedgevec(jj,1:3);xwedgevec(jj,4:6)],nvec1vec(jj,:));
   if thetas == 0
      closwedangvec(ii) = 0;
   else
      closwedangvec(ii) = 2*pi - thetas;
   end
   
end

%-------------------------------------------------------------------
% Now we check for duplicates of the edge definitions
% Some edge definitions might need to be swapped: if there are duplicate edges, 
% and one of them has closwedangvec = 0,
% then there is a mistake in the edge definition, so a swap will be made.

planenvecs = planedata.planeeqs(:,1:3);

edgecosinglelist = edgecorners(:,1)*multfac + edgecorners(:,2);
nsteps = diff([0;edgecosinglelist]);
ivduplicates = find(nsteps == 0);
nduplicates = length(ivduplicates);
if nduplicates > 0
	for ii = 1:nduplicates
		comb1 = ivduplicates(ii)-1;
		comb2 = comb1+1;
		if closwedangvec(comb1) == 0 || closwedangvec(comb2) == 0		% edge definitions should be swapped
			plane1 = planesatedge(comb1,2);
			plane2 = planesatedge(comb2,2);
			planesatedge(comb1,2) = plane2;
			planesatedge(comb2,2) = plane1;		
			temp1 = thirdpoint(comb1,:);
			temp2 = thirdpoint(comb2,:);
			thirdpoint(comb1,:) = temp2;
			thirdpoint(comb2,:) = temp1;
            
   			xwedge = [planedata.corners(edgecorners(comb1,1),:);planedata.corners(edgecorners(comb1,2),:)];
   			nvec1   = planenvecs( planesatedge(comb1,1),: );
%    			nvec1   = planedata.planeeqs( planesatedge(comb1,1),1:3 );

            xsou   = thirdpoint(comb1,:);
   			[~,thetas,~] = EDcoordtrans1(xsou,xwedge,nvec1);
   			if thetas == 0
    	  		closwedangvec(comb1) = 0;
   			else
    	  		closwedangvec(comb1) = 2*pi - thetas;
   			end

   			xwedge = [planedata.corners(edgecorners(comb2,1),:);planedata.corners(edgecorners(comb2,2),:)];
   			nvec1   = planenvecs( planesatedge(comb2,1),: );
            
   			xsou   = thirdpoint(comb2,:);
   			[~,thetas,~] = EDcoordtrans1(xsou,xwedge,nvec1);
   			if thetas == 0
    	  		closwedangvec(comb2) = 0;
   			else
    	  		closwedangvec(comb2) = 2*pi - thetas;
   			end

		end   % if closwedangvec(comb1) == 0 | 

	end		% for ii = 1:nduplicates
end

if nplanes < 256
    planesatedge = uint8(planesatedge);
elseif nplanes < 65535
    planesatedge = uint16(planesatedge);
else
    planesatedge = uint32(planesatedge);
end

%-------------------------------------------------------------------
% Find which edges each plane is connected to

% if showtext >= 4,
% 	disp('         find edgesatplane')
% end

edgesatplane = zeros(nplanes,double(max(ncornersperplanevec)));

for ii= 1:nplanes
	iv = find(planesatedge==ii);
	iv = iv - floor((iv-1)/nedges)*nedges;
	if ~isempty(iv)
		edgesatplane(ii,1:length(iv)) = iv.';
	end
end

% Check how many edges are defined for each plane. There must be at
% least three edges per plane.

tempnedgesperplanevec = sum(edgesatplane.'~=0).';

iv = find(tempnedgesperplanevec<3);

if ~isempty(iv)
    
    disp(' ')
    if showtext >= 4
        for ii = 1:length(iv)
            disp(['         Plane ',int2str(iv(ii)),' has only ',int2str(tempnedgesperplanevec(iv(ii))),' edges connected to it'])        
        end
        error('ERROR: The planes above have less than three edges.')
    else
        error('ERROR: Some planes have less than three edges. Set showtext >= 4 to see a list of those planes.')
    end
        
end

%-------------------------------------------------------------------
% Make a special list of thin planes

planeisthin = zerosvec1;
rearsideplane = zerosvec1;

iv = find(closwedangvec==0);
if ~isempty(iv)
	planenumberlist = planesatedge(iv,1);
	planeisthin(planenumberlist) = planeisthin(planenumberlist) + 1;
	rearsideplane(planenumberlist) = planesatedge(iv,2);
end

iv = find(closwedangvec == 0 & planesatedge(:,2)~=0);
if ~isempty(iv)
	planenumberlist = planesatedge(iv,2);
	planeisthin(planenumberlist)   = planeisthin(planenumberlist) + 1;		
	rearsideplane(planenumberlist) = planesatedge(iv,1); 
end

planeisthin = sign(planeisthin);
listofthinplanes = find(planeisthin);

%---------------------------------------------------------------
% Closwedangvec should be calculated into nyvec

nyvec = pi./(2*pi - closwedangvec);
integerny = ( abs(nyvec - round(nyvec)) < 1e-10 );

%-----------------------------------------------------------
% Construct some other useful edge variables

% edgestartcoords      = zerosvec5;
% edgestartcoordsnudge = zerosvec5;
% edgeendcoords        = zerosvec5;
% edgeendcoordsnudge   = zerosvec5;
% edgemidcoords        = zerosvec5;
% edgenvecs            = zerosvec5;

edgestartcoords = [planedata.corners(edgecorners(:,1),:)];
edgeendcoords   = [planedata.corners(edgecorners(:,2),:)];
% edgemidcoords   = (edgestartcoords+edgeendcoords)/2;
edgenvecs       = planedata.planeeqs(planesatedge(:,1),1:3);

% Changed normvec direction!!! 050504
% % %     edgenormvecs = edgestartcoords - edgeendcoords;
edgenormvecs = edgeendcoords - edgestartcoords;

edgestartcoordsnudge = edgestartcoords + 1e-10*edgenormvecs;
edgeendcoordsnudge   = edgeendcoords   - 1e-10*edgenormvecs;

edgelengthvec = sqrt(sum( ((edgenormvecs).^2).' )).';
edgenormvecs = edgenormvecs./edgelengthvec(:,ones(1,3));

%---------------------------------------------------------------
%
%		BACK TO SOME PURE PLANE RELATED PARAMETERS
%
%---------------------------------------------------------------
%
% First, make a list of which planes that are absorbing/non-reflective

planeisabsorbing = zerosvec1;

% listofactiveplanes = find(reflfactors~=0);
listofabsorbingplanes = find(reflfactors == 0);
planeisabsorbing(listofabsorbingplanes) = ones(size(listofabsorbingplanes));

%---------------------------------------------------------------
% Help variables: planealignedwplane and planeconnectstoplane
%
% Preparatory for the check of which planes see which plane: check
% if any planes are aligned with each other. If two planes are aligned,
% then planealignedwplane = 1 for that combination.
% However, if two planes are back-to back (thin planes), then
% planealignedwplane = 2 for that combination.
%
% Also make a matrix of which planes connect to which planes and via
% which edges.

planealignedwplane = zerosvec2;

% First, check which planes are aligned with each other

for ii = 1:nplanes-1
	if planeisabsorbing(ii) ~= 1
		oneplaneeq = planedata.planeeqs(ii,:);
		diffvec1 = oneplaneeq(ones(nplanes-ii,1),:) - planedata.planeeqs(ii+1:nplanes,:);
		diffvec2 = oneplaneeq(ones(nplanes-ii,1),:) + planedata.planeeqs(ii+1:nplanes,:);
		diffvec = min( [sum( diffvec1.'.^2 ).'   sum( diffvec2.'.^2 ).'].' ).';
		alignedplanes = find(diffvec < geomacc) + ii;
		if ~isempty(alignedplanes)
			planealignedwplane(alignedplanes,ii) = int8(double(planealignedwplane(alignedplanes,ii)) + double((planeisabsorbing(alignedplanes)~=1))); 
		end	
	end
end

planealignedwplane = (sparse(sign(double(planealignedwplane) + double(planealignedwplane).')));

% Second, check which planes are back to back

% % % ivec = find(planeisthin);   listofthinplanes
if ~isempty(listofthinplanes)
	for ii = 1:length(listofthinplanes)
		plane = listofthinplanes(ii);
		rearplane = rearsideplane(plane);
		if planeisabsorbing(plane) ~= 1 && planeisabsorbing(rearplane) ~= 1
			planealignedwplane(plane,rearplane) = 2;
			planealignedwplane(rearplane,plane) = 2;
		end
	end
end

planeconnectstoplane = zerosvec2;
clear zerosvec2

for ii = 1:nplanes
	if planeisabsorbing(ii) ~= 1
%      disp(['Plane number ',int2str(planenumbers(ii)),' (CAD) ',int2str(ii),' (ESIE2): '])
% 		edgelist = edgesatplane(ii,:);
 		edgelist = edgesatplane(ii,1:double(ncornersperplanevec(ii)));
		ivec = planesatedge(edgelist,:);
		ivec = reshape(ivec.',length(edgelist)*2,1);
		ivec = ivec(ivec~=ii);
		okplanes = find(planeisabsorbing(ivec)~=1);
		ivec = ivec(okplanes);
		if ~isempty(ivec)
			planeconnectstoplane(ii,ivec) = edgelist(okplanes);
		end
	end
end
    
%---------------------------------------------------------------
%
%			planeseesplane
%
%---------------------------------------------------------------
% Check which planes see which planes:
% 1. If one of the planes has reflfactors = 0, then the two can't see each other.
% 2. Reflective planes that are aligned with each other can not reach each other
% 3. Reflective planes that are not aligned but have the same normal vectors
%    can not reach each other
% 4. Planes that have all its corners behind another plane can not be seen by that plane.
%
% planeseesplane = 0 means that either, one of the planes is non-reflective or, that
%					 two planes never can see each other, but they are not aligned with each other.
% planeseesplane = 1 means that two reflective planes might see each other, but there could
%                    be obstructions
% planeseesplane = -1 means that two reflective planes are aligned (but they are not back-to-back)
%					  and thus they can never see each other.
% planeseesplane = -2 means that two reflective planes are back-to-back with each other and thus
%					  they can never see each other.

if showtext >= 4
	disp('         Check which planes see which planes')
end

planeseesplane = int8(ones(nplanes,nplanes));

% 1. If one of the planes has reflfactors = 0, then the two can't see each other.

if ~isempty(listofabsorbingplanes)
	zerosvec = zeros(length(listofabsorbingplanes),nplanes);
	planeseesplane(listofabsorbingplanes,:) = zerosvec;
	planeseesplane(:,listofabsorbingplanes) = zerosvec.';
end

% 2. Reflective planes that are aligned with each other can not reach each other

ivec = find( planealignedwplane == 1);
if ~isempty(ivec)
	planeseesplane(ivec) = - 1*ones(size(ivec));
end

ivec = find( planealignedwplane == 2);
if ~isempty(ivec)
	planeseesplane(ivec) = - 2*ones(size(ivec));
end

% 3. Reflective planes that have the same normal vectors can not reach each other

numvec = (1:nplanes);
for ii = 1:nplanes
	if planeisabsorbing(ii) ~= 1
		nvec1 = planedata.planeeqs(ii,1:3);
		ivec = find(planealignedwplane(ii,:)==0 & planeisabsorbing.'==0 & numvec>ii);
		if ~isempty(ivec)
			similarvec = abs( nvec1(ones(length(ivec),1),:) - planedata.planeeqs(ivec,1:3)) < geomacc;	
			similarvec = prod(  double(similarvec.') ).';
			ivec2 = ivec(similarvec==1);
			if ~isempty(ivec2)
				zerosvec = zeros(size(ivec2));
				planeseesplane(ii,ivec2) = zerosvec;	
				planeseesplane(ivec2,ii) = zerosvec.';
			end		
		end
	end
end

% 3.5 A plane does not see itself

%planeseesplane = int8(double(planeseesplane).*(1-eye(nplanes)));
iv = (0:nplanes-1)*nplanes + (1:nplanes);
planeseesplane(iv) = zeros(size(iv));

% 4. Planes that have all its corners behind another plane can not be seen by that plane.
%    First, we construct a list referring to the complete matrix of size [nplanes,nplanes]
%    The index vector is called iv. Then we find which combinations that
%    could be cleared. After having cleared a number of combinations (in
%    iv, fromplanenumb and toplanenumb) we pick out out only the non-zero
%    indices in iv.

if nplanes*nplanes < 128
    iv = int8([1:nplanes*nplanes].');                
elseif nplanes*nplanes < 32768
    iv = int16([1:nplanes*nplanes].');                
else
    iv = int32([1:nplanes*nplanes].');                
end
iv = reshape(iv,nplanes,nplanes);

% If there are any absorbing planes, we zero all plane-to-plane
% combinations that involve such a plane

if ~isempty(listofabsorbingplanes)
    nabsorbingplanes = length(listofabsorbingplanes);
    iv(listofabsorbingplanes,:) = zeros(nabsorbingplanes,nplanes);
    iv(:,listofabsorbingplanes) = zeros(nplanes,nabsorbingplanes);
end

% Check connecting planes first: if two planes are connected and their
% shared edge has a closwedangvec > pi, then the two planes can see each
% other. So, we can zero the plane-pair combinations where closwedangvec
% <= pi but larger than zero.

edgecombstozero = find(closwedangvec<=pi & closwedangvec > 0);
if ~isempty(edgecombstozero)
    ivzerocombs = (double(planesatedge(edgecombstozero,1))-1)*nplanes + double(planesatedge(edgecombstozero,2));
    iv(ivzerocombs) = zeros(size(ivzerocombs));
    ivzerocombs = (double(planesatedge(edgecombstozero,2))-1)*nplanes + double(planesatedge(edgecombstozero,1));
    iv(ivzerocombs) = zeros(size(ivzerocombs));
end

% Now we should keep only the index numbers in iv that are still non-zero
% We also take the chance to zero planeseesplane for the combinations that
% could be cleared.
ivclear = uint32(find(iv==0));
planeseesplane(ivclear) = zeros(size(ivclear));
iv(ivclear) = [];

% Of the remaining plane-to-plane combinations, we pick out the ones 
% for which we need to check if the "to-plane" is in front of the "from-plane"
iv = iv(planeseesplane(iv)==1 & planeconnectstoplane(iv)==0);

% We create full matrices, fromplanenumb and toplanenumb, and then keep
% only the combinations that remain in iv,
% to keep fromplanenumb and toplanenumb aligned with iv.

% To-plane numbers is given by the row no.
if nplanes < 256
    toplanenumb = uint8([1:nplanes].');
elseif nplanes < 65536
    toplanenumb = uint16([1:nplanes].');
else
    toplanenumb = uint32([1:nplanes].');
end
toplanenumb = reshape(toplanenumb(:,ones(1,nplanes)),nplanes*nplanes,1);
toplanenumb = toplanenumb(iv);

% From-plane numbers is given by the col no.
if nplanes < 256
    fromplanenumb = uint8(1:nplanes);
elseif nplanes < 65536
    fromplanenumb = uint16(1:nplanes);
else
    fromplanenumb = uint32(1:nplanes);
end
fromplanenumb = reshape(fromplanenumb(ones(nplanes,1),:),nplanes*nplanes,1);
fromplanenumb = fromplanenumb(iv);

% Now, if a plane has *all* its corners behind another plane, those two
% planes can not see each other.
% We check the corners of the "to-plane" (columns) and check if they are in front of
% the "from-plane" (rows).
% The ivreftocoplamatrix is the index number to the
% cornerinfrontofplanematrix, which has the size [nplanes,ncorners].
% NB! In order to save some space, we don't construct the ivreftocoplamatrix specifically.

% First we check the first 3 corners since all planes have at least 3
% corners.
% Plane corner 1
% ivreftocoplamatrix = double(fromplanenumb) + ( double(planecorners( toplanenumb,1 ))-1 )*nplanes;
% cornersbehind = cornerinfrontofplane(ivreftocoplamatrix)<1;
cornersbehind = planedata.cornerinfrontofplane(double(fromplanenumb) + ( double(planecorners( toplanenumb,1 ))-1 )*nplanes)<1;
% Plane corner 2
% ivreftocoplamatrix = double(fromplanenumb) + ( double(planecorners( toplanenumb,2 ))-1 )*nplanes;
% cornersbehind = cornersbehind &(cornerinfrontofplane(ivreftocoplamatrix) < 1);
cornersbehind = cornersbehind &(planedata.cornerinfrontofplane(double(fromplanenumb) + ( double(planecorners( toplanenumb,2 ))-1 )*nplanes) < 1);
% Plane corner 3
% ivreftocoplamatrix = double(fromplanenumb) + ( double(planecorners( toplanenumb,3 ))-1 )*nplanes;
% cornersbehind = cornersbehind &(cornerinfrontofplane(ivreftocoplamatrix) < 1);
cornersbehind = cornersbehind &(planedata.cornerinfrontofplane(double(fromplanenumb) + ( double(planecorners( toplanenumb,3 ))-1 )*nplanes) < 1);

if minncornersperplane == 4 && maxncornersperplane == 4
%   Plane corner 4
%     ivreftocoplamatrix = double(fromplanenumb) + ( double(planecorners( toplanenumb,4 ))-1 )*nplanes;
%     cornersbehind = cornersbehind &(cornerinfrontofplane(ivreftocoplamatrix) < 1);
    cornersbehind = cornersbehind &(planedata.cornerinfrontofplane(double(fromplanenumb) + ( double(planecorners( toplanenumb,4 ))-1 )*nplanes) < 1);
elseif not(minncornersperplane == 3 & maxncornersperplane == 3)
    ivall = uint32([1:length(fromplanenumb)].');
    if minncornersperplane == 3
        iv3cornerplanes = ncornersperplanevec(toplanenumb)==3;
        ivall(iv3cornerplanes) = [];
    end
    if maxncornersperplane == 4
% Plane corner 4
%        ivreftocoplamatrix = double(fromplanenumb(ivall)) + ( double(planecorners( toplanenumb(ivall),4 ))-1 )*nplanes;
%        cornersbehind(ivall) = cornersbehind(ivall) &(cornerinfrontofplane(ivreftocoplamatrix) < 1);
        cornersbehind(ivall) = cornersbehind(ivall) &(planedata.cornerinfrontofplane(double(fromplanenumb(ivall)) + ( double(planecorners( toplanenumb(ivall),4 ))-1 )*nplanes) < 1);
    else
        ivmorethan4 = find(ncornersperplanevec(toplanenumb(ivall))>4);
        temp = ivall(ivmorethan4);
        ivall(ivmorethan4) = [];
        ivmorethan4 = temp;
        cornersbehind(ivall) = cornersbehind(ivall) &(planedata.cornerinfrontofplane(double(fromplanenumb(ivall)) + ( double(planecorners( toplanenumb(ivall),4 ))-1 )*nplanes) < 1);
        for ii = 5:maxncornersperplane
            ivselection = find(ncornersperplanevec(toplanenumb(ivmorethan4))>=ii);
            cornersbehind(ivmorethan4(ivselection)) = cornersbehind(ivmorethan4(ivselection)) &(planedata.cornerinfrontofplane(double(fromplanenumb(ivmorethan4(ivselection))) + ( double(planecorners( toplanenumb(ivmorethan4(ivselection)),ii ))-1 )*nplanes) < 1);
        end
        clear ivselection
    end
end
clear toplanenumb fromplanenumb

ivclear = iv(cornersbehind==1);
clear iv

planeseesplane(ivclear) = zeros(size(ivclear));
clear ivclear

temp1 = (planeseesplane~=0);
temp2 = (planeseesplane.'~=0);
planeseesplane = int8(double(planeseesplane).*(temp1.*temp2));

% New addition:
% If the user has asked for it, check plane-mid-point to plane-mid-point
% for obstruction. If there are obstructions, set planeseesplane to 0 for
% that combination.

if planeseesplanestrategy == 1
	ivorig = uint32(find(planeseesplane==1));
    iv = ivorig;
	fromplane = ceil(double(iv)/nplanes);
	toplane = uint16(double(iv) - (fromplane-1)*nplanes);
    fromplane = uint16(fromplane);
	listofplanesthatcanobscure = find( sum(planeseesplane==1) );

    iv4planes = find(ncornersperplanevec == 4);
    planemidpoints = zeros(nplanes,3);
    
    if any(ncornersperplanevec~=4)
        disp('WARNING: planeseesplanestrategy 1 not implemented for models with non-4-corner planes')
        planeseesplanestrategy = 0;
    else
        xvalues = planedata.corners(planecorners.',1);
        xvalues = reshape(xvalues,4,length(xvalues)/4);
        yvalues = planedata.corners(planecorners.',2);
        yvalues = reshape(yvalues,4,length(yvalues)/4);
        zvalues = planedata.corners(planecorners.',3);
        zvalues = reshape(zvalues,4,length(zvalues)/4);
        planemidpoints = [mean(xvalues).' mean(yvalues).' mean(zvalues).'];
    end
end

if planeseesplanestrategy == 1
	for ii = listofplanesthatcanobscure
        planeobsc = ii;
        ivcheckvis1 = uint32( planeobsc + (double(fromplane)-1)*nplanes);
        ivcheckvis2 = uint32( planeobsc + (double(toplane)-1)*nplanes);        
        iv3 = find(toplane~=planeobsc & fromplane~=planeobsc & ((planeseesplane(ivcheckvis1) == 1) | (planeseesplane(ivcheckvis2) == 1)) & (planeseesplane(ivcheckvis1) >= 0) & (planeseesplane(ivcheckvis2) >= 0) );
        tempcanplaneobstruct = zerosvec1.';
        tempcanplaneobstruct(planeobsc) = 1;
        [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstrpaths(planemidpoints(fromplane(iv3),:),planemidpoints(toplane(iv3),:),[],[],tempcanplaneobstruct,planeseesplane,...
            planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);
        if length(iv3) > length(nonobstructedpaths)
            iv3(nonobstructedpaths) = [];
            iv(iv3) = [];
	        toplane(iv3) = [];
	        fromplane(iv3) = [];
        end
%        disp(['Obsc.plane ',int2str(planeobsc),': ',int2str(length(iv)),' p2p combos left'])        
	end
    planeseesplane(ivorig) = 0;    
    planeseesplane(iv) = 1;    
    clear iv ivorig toplane fromplane ivcheckvis ivcheckvis2
end

%---------------------------------------------------------------
%
%		EDGE RELATED PARAMETERS AGAIN
%
%---------------------------------------------------------------
%
% Which edges should be switched off?
% (1) If the list listofedgestoskip contains values, then the variable
%     firstcornerskip and the list listofcornerstoskip are ignored. If
%     listofedgestoskip is empty, then step (2) will be carried out.
% (2) Go through the list edgecorners. If any edge contains one of the
%     corners, whose numbers are in switchoffcorners, then that edge should be
%     stored in the list offedges.
% (3) If an edge has an integer ny-value, then it should be switched off.
% (4) If one of the two connecting planes has reflfactors = 0, then it should
%     be switched off.

offedges = zerosvec3;

% (1) If the list listofedgestoskip contains values, then the variable
%     firstcornerskip and the list listofcornerstoskip are ignored. If
%     listofedgestoskip is empty, then step (2) will be carried out.
% (2) Go through the list edgecorners. If any edge contains one of the
%     corners, whose numbers are in switchoffcorners, then that edge should be
%     stored in the list offedges.

if ~isempty(listofedgestoskip)   
    ncols = size(listofedgestoskip,2);
    if ncols == 1  % listofedgestoskip simply contains edge numbers
        if max(listofedgestoskip) > nedges
            error('ERROR: listofedgestoskip contains an edge number which is too large.')
        end
        if min(listofedgestoskip) < 1
            error('ERROR: listofedgestoskip contains an edge number which is < 1')
        end
        offedges(listofedgestoskip) = 1;
    else  % listofedgestoskip contains corner pairs
        if max(max(listofedgestoskip)) > ncorners
            error('ERROR: listofedgestoskip contains a corner number which is too large.')
        end
        if min(min(listofedgestoskip)) < 1
            error('ERROR: listofedgestoskip contains a corner number which is < 1')
        end
        for ii = 1:size(listofedgestoskip,1)
            cornerpair = listofedgestoskip(ii,:);
            [lia,locb] = ismember(cornerpair,edgecorners,'rows');
            if lia == 1
                offedges(locb) = 1;
            else
                [lia,locb] = ismember(cornerpair([2 1]),edgecorners,'rows');
                if lia == 1
                    offedges(locb) = 1;               
                end
            end        
        end

    end
else
    if ~isempty(switchoffcorners)   
  	    for ii = 1:nedges
	       corner1 = edgecorners(ii,1);
	       corner2 = edgecorners(ii,2);
	       remove = sum( corner1 == switchoffcorners ) + sum( corner2 == switchoffcorners );
	       offedges(ii) = sign(remove);
	    end
    end
end

% (2) If an edge has an integer ny-value, then it should be switched off.

offedges = offedges + integerny;

% (3) If one of the two connecting planes has reflfactors = 0, then it should
%     be switched off.

reflectingedges = reshape(reflfactors(planesatedge),nedges,2);
reflectingedges = reflectingedges(:,1).*reflectingedges(:,2);
offedges = offedges + (1-reflectingedges);

offedges = find(offedges ~= 0);

%--------------------------------------------------------------------------------
%
%			edgeseesplane
%
%--------------------------------------------------------------------------------
%
% Now check which planes each edge can see. 
%    1. Set all combinations with reflfactors = 0 to -3
%    2. Switch off all combinations with offedges
% 	 3. For all edges that are aligned with a plane, set edgeseesplane to:
%			-1 if the edge doesn't belong to the plane
%			-2 if the edge does belong to the plane
%    4. If at least one edge end point is in front of a plane, then the
%       plane is potentially visible (edgeseesplane = 1)
%    5. If at least one corner point is in front of:
%           *one* of the two edge planes, when closwedang < pi
%           *both* edge planes, when closwedangvec > pi

if showtext >= 4
	disp('         Check which edges see which planes')
end

edgeseesplane = int8(ones(nplanes,nedges));

%    1. Set all edges belonging to planes with reflfactors = 0 to
%    edgeseesplane = -3.

% % if ~isempty(listofabsorbingplanes),
% % 	edgeseesplane(listofabsorbingplanes,:) = -3*ones(length(listofabsorbingplanes),nedges);
% % end

%    2. Switch off all combinations with offedges

if ~isempty(offedges)
	edgeseesplane(:,offedges) = zeros(nplanes,length(offedges));
end

% 	3.  For all edges that belong to a plane, set edgeseesplane to -2
%    	Also, for all other edges that are aligned with a plane, set
%    	edgeseesplane to -1

edgelist = [1:nedges].';
edgelist(offedges) = [];

plane1list = planesatedge(edgelist,1); 
plane2list = planesatedge(edgelist,2); 
indexvec = uint32( double(plane1list) + (edgelist-1)*nplanes);
edgeseesplane(indexvec) = -2*(reflfactors(plane1list)~=0);
indexvec = uint32( double(plane2list) + (edgelist-1)*nplanes);
edgeseesplane(indexvec) = -2*(reflfactors(plane2list)~=0);

for ii = 1:length(edgelist)

	aligners1 = find(planealignedwplane(:,plane1list(ii))==1);
	if ~isempty(aligners1)
		edgeseesplane(aligners1,edgelist(ii)) = -1*ones(size(aligners1));
	end

	aligners2 = find(planealignedwplane(:,plane2list(ii))==1);
	if ~isempty(aligners2)
		edgeseesplane(aligners2,edgelist(ii)) = -1*ones(size(aligners2));
	end

end

%    4. If at least one edge end point is in front of a plane, then the
%       plane is potentially visible.
%       Do this for all plane-edge combinations for which edgeseesplane = 1
%       up til now.
%    5. If at least one corner point is in front of:
%           *one* of the two edge planes, when closwedang < pi
%           *both* edge planes, when closwedangvec > pi

ntot = nplanes*nedges;
ivclear = find(edgeseesplane~=1);
% Edge number is given by the col no.
if nedges < 256
    edgenumb = uint8(1:nedges);
elseif nedges < 65536
    edgenumb = uint16(1:nedges);
else
    edgenumb = uint32(1:nedges);        
end
edgenumb = reshape(edgenumb(ones(nplanes,1),:),nplanes*nedges,1);
edgenumb(ivclear) = [];
% Plane numbers is given by the row no.
if nplanes < 256
    planenumb = uint8([1:nplanes].');
elseif nplanes < 65536
    planenumb = uint16([1:nplanes].');
else
    planenumb = uint32([1:nplanes].');
end
planenumb = reshape(planenumb(:,ones(1,nedges)),nplanes*nedges,1);
planenumb(ivclear) = [];

if nplanes*nedges < 128
    iv = int8([1:nplanes*nedges].');                
elseif nplanes*nedges < 32768
    iv = int16([1:nplanes*nedges].');                
else
    iv = int32([1:nplanes*nedges].');                
end
iv(ivclear) = [];

if ~isempty(edgenumb)
	% For the "edgecorners in front of plane"-test,
	% we can use the cornerinfrontofplane matrix, but we need to construct
	% index vectors, based on the planenumb and edgenumb values.
	
	% ivreftocoplamatrix = double(planenumb) + ( double( edgecorners(edgenumb,1) )-1 )*nplanes;
	% edgeinfrontofplane = cornerinfrontofplane( ivreftocoplamatrix );
	% ivreftocoplamatrix = double(planenumb) + ( double( edgecorners(edgenumb,2) )-1 )*nplanes;
	% edgeinfrontofplane = edgeinfrontofplane | cornerinfrontofplane( ivreftocoplamatrix );
	% edgeseesplane(iv) = edgeseesplane(iv) & edgeinfrontofplane; 
	
	edgeinfrontofplane =                      planedata.cornerinfrontofplane( double(planenumb) + ( double( edgecorners(edgenumb,1) )-1 )*nplanes );
	edgeinfrontofplane = (edgeinfrontofplane==1) | (planedata.cornerinfrontofplane( double(planenumb) + ( double( edgecorners(edgenumb,2) )-1 )*nplanes )==1);
	edgeseesplane(iv) = int8(double(edgeseesplane(iv)).*double(edgeinfrontofplane)); 
	
	%edgeseesplane(iv) = int8(   double(edgeseesplane(iv)).*(  (EDinfrontofplane(corners(edgecorners(edgenumb,1),:),planenvecs(planenumb,:),corners(planecorners(planenumb,1),:),corners(planecorners(planenumb,2),:)) >0) | ( EDinfrontofplane(corners(edgecorners(edgenumb,2),:),planenvecs(planenumb,:),corners(planecorners(planenumb,1),:),corners(planecorners(planenumb,2),:)) >0 ) ));
	clear edgeinfrontofplane
end

if ~isempty(edgenumb)
	% For the "planecorners in front of edge-plane"-test,
	% we can use the cornerinfrontofplane matrix, but we need to construct
	% index vectors, based on the planenumb and edgenumb values.
	%
	% First we split up the iv into two halves: the one with edges < pi and the
	% one with edges > pi;
	
	ivsmaller = uint32(find(closwedangvec(edgenumb)<pi));
	ivlarger = uint32([1:length(edgenumb)].');
	ivlarger(ivsmaller) = [];
	% ivsmaller = iv(ivsmaller);
	% ivlarger = iv(ivlarger);

    if ~isempty(ivsmaller)
		% First the edges smaller than pi (closwedang < pi): at least one plane corner needs
		% to be in front of one or the other plane.
		% Plane corner 1
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),1)  ) + ( double(planecorners( planenumb(ivsmaller),1 ))-1 )*nplanes;
		planecornerinfrontofedge =                            planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),2)  ) + ( double(planecorners( planenumb(ivsmaller),1 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		% Plane corner 2
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),1)  ) + ( double(planecorners( planenumb(ivsmaller),2 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),2)  ) + ( double(planecorners( planenumb(ivsmaller),2 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		% Plane corner 3
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),1)  ) + ( double(planecorners( planenumb(ivsmaller),3 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),2)  ) + ( double(planecorners( planenumb(ivsmaller),3 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		
		if minncornersperplane == 4 && maxncornersperplane == 4
            % Plane corner 4
			ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),1)  ) + ( double(planecorners( planenumb(ivsmaller),4 ))-1 )*nplanes;
			planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
			ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller),2)  ) + ( double(planecorners( planenumb(ivsmaller),4 ))-1 )*nplanes;
			planecornerinfrontofedge = planecornerinfrontofedge | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		elseif not(minncornersperplane == 3 & maxncornersperplane == 3)
            ivall = uint32(1:length(ivsmaller));
            if minncornersperplane == 3
                iv3cornerplanes = ncornersperplanevec(planenumb(ivsmaller))==3;
                ivall(iv3cornerplanes) = [];
            end
            if maxncornersperplane == 4
                % Plane corner 4
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivall)),1)  ) + ( double(planecorners( planenumb(ivsmaller(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivall)),2)  ) + ( double(planecorners( planenumb(ivsmaller(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
            else
                ivmorethan4 = find(ncornersperplanevec(planenumb(ivsmaller(ivall)))>4);
                temp = ivall(ivmorethan4);
                ivall(ivmorethan4) = [];
                ivmorethan4 = temp;
                % Plane corner 4
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivall)),1)  ) + ( double(planecorners( planenumb(ivsmaller(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivall)),2)  ) + ( double(planecorners( planenumb(ivsmaller(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
                for ii = 5:maxncornersperplane
                    % Plane corner 5,6,7,etc
%                     ivselection = find(ncornersperplanevec(planenumb(ivsmaller(ivmorethan4)))>=ii);
					ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivmorethan4)),1)  ) + ( double(planecorners( planenumb(ivsmaller(ivmorethan4)),4 ))-1 )*nplanes;
					planecornerinfrontofedge(ivmorethan4) = planecornerinfrontofedge(ivmorethan4) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
					ivreftocoplamatrix = double(  planesatedge(edgenumb(ivsmaller(ivmorethan4)),2)  ) + ( double(planecorners( planenumb(ivsmaller(ivmorethan4)),4 ))-1 )*nplanes;
					planecornerinfrontofedge(ivmorethan4) = planecornerinfrontofedge(ivmorethan4) | (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
                end
                clear ivselection
            end
		end
		clear ivreftocoplamatrix
		
		ivclear = ivsmaller(planecornerinfrontofedge==0);
		edgeseesplane(iv(ivclear)) = 0;
    end
    
    if ~isempty(ivlarger)
		% Second the edges larger than pi (closwedang > pi): at least one plane corner needs
		% to be in front of *both* edge planes.
		% Plane corner 1
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),1)  ) + ( double(planecorners( planenumb(ivlarger),1 ))-1 )*nplanes;
		planecornerinfrontofedge =                            planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),2)  ) + ( double(planecorners( planenumb(ivlarger),1 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1);
		% Plane corner 2
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),1)  ) + ( double(planecorners( planenumb(ivlarger),2 ))-1 )*nplanes;
		temp =                                                planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),2)  ) + ( double(planecorners( planenumb(ivlarger),2 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
		% Plane corner 3
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),1)  ) + ( double(planecorners( planenumb(ivlarger),3 ))-1 )*nplanes;
		temp =                                                planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
		ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),2)  ) + ( double(planecorners( planenumb(ivlarger),3 ))-1 )*nplanes;
		planecornerinfrontofedge = planecornerinfrontofedge | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
		
		if minncornersperplane == 4 && maxncornersperplane == 4
            % Plane corner 4
			ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),1)  ) + ( double(planecorners( planenumb(ivlarger),4 ))-1 )*nplanes;
			temp =                                                planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
			ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger),2)  ) + ( double(planecorners( planenumb(ivlarger),4 ))-1 )*nplanes;
			planecornerinfrontofedge = planecornerinfrontofedge | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
		elseif not(minncornersperplane == 3 & maxncornersperplane == 3)
            ivall = uint32(1:length(ivlarger));
            if minncornersperplane == 3
                iv3cornerplanes = ncornersperplanevec(planenumb(ivlarger))==3;
                ivall(iv3cornerplanes) = [];
            end
            if maxncornersperplane == 4
                % Plane corner 4
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivall)),1)  ) + ( double(planecorners( planenumb(ivlarger(ivall)),4 ))-1 )*nplanes;
				temp =                                                              planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivall)),2)  ) + ( double(planecorners( planenumb(ivlarger(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
            else
                ivmorethan4 = find(ncornersperplanevec(planenumb(ivlarger(ivall)))>4);
                temp = ivall(ivmorethan4);
                ivall(ivmorethan4) = [];
                ivmorethan4 = temp;
              % Plane corner 4
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivall)),1)  ) + ( double(planecorners( planenumb(ivlarger(ivall)),4 ))-1 )*nplanes;
				temp =                                                              planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
				ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivall)),2)  ) + ( double(planecorners( planenumb(ivlarger(ivall)),4 ))-1 )*nplanes;
				planecornerinfrontofedge(ivall) = planecornerinfrontofedge(ivall) | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
                for ii = 5:maxncornersperplane
                    % Plane corner 5,6,7,etc
%                     ivselection = find(ncornersperplanevec(planenumb(ivlarger(ivmorethan4)))>=ii);
					ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivmorethan4)),1)  ) + ( double(planecorners( planenumb(ivlarger(ivmorethan4)),4 ))-1 )*nplanes;
					temp =                                                          planedata.cornerinfrontofplane(ivreftocoplamatrix)==1;
					ivreftocoplamatrix = double(  planesatedge(edgenumb(ivlarger(ivmorethan4)),2)  ) + ( double(planecorners( planenumb(ivlarger(ivmorethan4)),4 ))-1 )*nplanes;
					planecornerinfrontofedge(ivmorethan4) = planecornerinfrontofedge(ivmorethan4) | (temp & (planedata.cornerinfrontofplane(ivreftocoplamatrix)==1));
                end
                clear ivselection
            end
		end
		
		ivclear = ivlarger(planecornerinfrontofedge==0);
		edgeseesplane(iv(ivclear)) = 0;
    end
end

clear ivlarger ivsmaller ivreftocoplamatrix edgenumb planenumb iv

%----------------------------------------------------------------------------
%
%           indentingedgepairs
%
%--------------------------------------------------------------------------------
%
% If there are any planes with indents, make a list of all consecutive edge pairs that
% form an indent, because such edge pairs could never see each other.

indentingedgepairs = [];

if any(planedata.planehasindents)
    iv = find(planedata.indentingcorners);
    indentingedgepairs = zeros(length(iv),2);
    indentingplanenumbers = mod(iv,nplanes);
    ivzeros = find(indentingplanenumbers==0);
    if ~isempty(ivzeros)
       indentingplanenumbers(ivzeros) = nplanes; 
    end
    conumbers = (iv - indentingplanenumbers)/nplanes+1;
    
    % We expand planecorners temporarily, cyclically
    planecorners = [planecorners planecorners(:,1:2)];
    
    for ii = 1:length(iv)
        edge1cornerpair = sort(planecorners(indentingplanenumbers(ii),conumbers(ii):conumbers(ii)+1));
        edge2cornerpair = sort(planecorners(indentingplanenumbers(ii),conumbers(ii)+1:conumbers(ii)+2));
        
        [~,edge1] = ismember(edge1cornerpair,edgecorners,'rows');
        [~,edge2] = ismember(edge2cornerpair,edgecorners,'rows');
        
        edgepair = sort([edge1 edge2]);
        indentingedgepairs(ii,:) = edgepair;
    
    end
    
    % Remove the temporary expansion
    
    planecorners = planecorners(:,1:size(planecorners,2)-2);
end

%----------------------------------------------------------------------------
%
%           canplaneobstruct
%
%--------------------------------------------------------------------------------
%
% Check which planes that are potentially obstructing:
%
% For an internal, convex problem, no planes can obstruct
%
% For exterior problems, all active planes can potentially obstruct.

canplaneobstruct = [];

if ~isempty(strfind(planedata.modeltype,'convex_ext')) || ~isempty(strfind(planedata.modeltype,'plate'))
%     canplaneobstruct = ones(1,nplanes,int8);
    canplaneobstruct = int8(reflfactors.'~=0);
elseif strfind(planedata.modeltype,'convex_int')
    % Bug fixed 15 Mar 2018
    % The line below used "nplanes,int8);" instead of "nplanes,'int8');"     
    canplaneobstruct = zeros(1,nplanes,'int8');   
end

% For an interior problem we know for sure that if a plane 
% has all other points in front of itself, it can not obstruct.
% NB! We need to check only points that do not belong to planes
% that are aligned with the plane, or planes that are non-active.

% if int_or_ext_model == 'e',
%     canplaneobstruct = (reflfactors.'~=0);
% else
    if isempty(canplaneobstruct)
    
    canplaneobstruct = (reflfactors.'~=0);
    maskvec1 = canplaneobstruct.';
    listofactiveplanes = find(canplaneobstruct);
    zerosvec = zeros(nplanes,1);
    
    for ii = 1:length(listofactiveplanes)
        iplane = listofactiveplanes(ii);
        maskvec2 = zerosvec;
        maskvec2(iplane) = 1;
        otherplanestocheck = find((double(planeseesplane(:,iplane))+maskvec2)~=1 & planeseesplane(:,iplane)~=-1 & maskvec1~=0);
        if ~isempty(otherplanestocheck)      
            cornerstocheck = unique(planecorners(otherplanestocheck,:));
            cornerstocheck = setdiff(unique(cornerstocheck),planecorners(iplane,1:ncornersperplanevec(iplane)));
%             pointinfront = EDinfrontofplane(corners(cornerstocheck,:),planenvecs(iplane,:),corners(planecorners(iplane,1),:),corners(planecorners(iplane,2),:));
            pointinfront = EDinfrontofplane(planedata.corners(cornerstocheck,:),planenvecs(iplane,:),planedata.corners(planecorners(iplane,1),:),planedata.corners(planecorners(iplane,2),:));
            if isempty(find(pointinfront==-1, 1))
                canplaneobstruct(iplane) = 0;    
            end
        else
            canplaneobstruct(iplane) = 0;    
        end
        
    end
end

%---------------------------------------------------------------
% For each thin plane, make sure that all edges have the same normal
% vector. If not, switch the direction of that edge, and all the other
% relevant parameters.

ivthin = find(planeisthin);

for ii = 1:length(ivthin)
    plane = ivthin(ii);
    if rearsideplane(plane) > plane
        edgelist = unique(edgesatplane(plane,1:nedgesperplanevec(plane)));
        iv = find(closwedangvec(edgelist,:)==0);
        if ~isempty(iv)
            edgelist = edgelist(iv);
            nedgestemp = length(edgelist);
            edgenveclist = edgenvecs(edgelist,:);
            refnvec = edgenveclist(1,:);
            nvecdiff = edgenveclist(2:nedgestemp,:) - refnvec(ones(nedgestemp-1,1),:);            
            nvecdiff = sum(abs(nvecdiff.')).';
            ivswitch = find(nvecdiff)+1;
            if ~isempty(ivswitch)
                edgenumber = edgelist(ivswitch);
                edgecorners(edgenumber,:)  = [edgecorners(edgenumber,2) edgecorners(edgenumber,1)];
                planesatedge(edgenumber,:) = [planesatedge(edgenumber,2) planesatedge(edgenumber,1)];
                edgenvecs(edgenumber,:) = -edgenvecs(edgenumber,:);
                edgestartcoords(edgenumber,:) = planedata.corners(edgecorners(edgenumber,1),:);
                edgeendcoords(edgenumber,:)   = planedata.corners(edgecorners(edgenumber,2),:);
                tempvec = edgeendcoordsnudge(edgenumber,:);    
                edgeendcoordsnudge(edgenumber,:) = edgestartcoordsnudge(edgenumber,:);
                edgestartcoordsnudge(edgenumber,:) = tempvec;
            end
        end        
    end

end

%-------------------------------------------------------------------
% Now that edges have been settled, find their 
% edgerelatedcoordsysmatrices

edgerelatedcoordsysmatrices = zeros(length(ivec),9);

for ii = 1:length(closwedangvec)
    xneworigo = edgestartcoords(ii,:);
    xknown1 = edgeendcoords(ii,:) - xneworigo;
    xknown1 = xknown1 / norm(xknown1);
    nvec1 = edgenvecs(ii,:);    
    A = [nvec1(2)*xknown1(3)-nvec1(3)*xknown1(2) ; nvec1(3)*xknown1(1)-nvec1(1)*xknown1(3) ; nvec1(1)*xknown1(2)-nvec1(2)*xknown1(1)];
    A = inv([xknown1.' A nvec1.']);  
    Bmatrix = A([2 3 1],:);
    edgerelatedcoordsysmatrices(ii,:) = reshape(Bmatrix,1,9);   
end

% The section below was removed 22 Jan 2018
% 
% % % %-------------------------------------------------------------------
% % % % Add some fields
% % % 
% % % ismodelconvex = 1;
% % % if any(any(edgeseesplane+2)) 
% % %     ismodelconvex = 0;
% % % end
% % % if any(planedata.planehasindents)
% % %     ismodelconvex = 0;
% % % end
% % % 
% % % if prod(planeisthin) == 1
% % %     ismodelthinplate = 1;
% % %    if any(any(sign(edgeseesplane)+1))
% % %        ismodelthinplate = 0;        
% % %    end
% % % else
% % %     ismodelthinplate = 0;
% % % end

%----------------------------------------------------------------------------
%
%		STORE THE VARIABLES IN THE OUTPUT STRUCTS
%
%----------------------------------------------------------------------------
%
% Also the variables from cadgeo??
% Yes: planeeqs planenvecs ncornersperplanevec planecorners corners
%      minvals maxvals

edgeseesplane = int8(edgeseesplane);
planeseesplane = int8(planeseesplane);

% if ncorners < 256
%     planecorners = uint8(planecorners);
% elseif ncorners < 65536
%     planecorners = uint16(planecorners);    
% end   
planeisthin = uint8(planeisthin);

% if max(ncornersperplanevec) <= 255
%     ncornersperplanevec = uint8(ncornersperplanevec);
% else
%     ncornersperplanevec = uint16(ncornersperplanevec);
% end

edgedata = struct('edgecorners',edgecorners,'planesatedge',planesatedge,...
    'closwedangvec',closwedangvec,'indentingedgepairs',indentingedgepairs,...
    'edgeseesplane',edgeseesplane,'edgestartcoords',edgestartcoords,...
    'edgeendcoords',edgeendcoords,'edgenvecs',edgenvecs,...
    'edgesatplane',edgesatplane,'edgelengthvec',edgelengthvec,...
    'offedges',offedges,'edgerelatedcoordsysmatrices',edgerelatedcoordsysmatrices,...
    'edgenormvecs',edgenormvecs,...
    'edgestartcoordsnudge',edgestartcoordsnudge,'edgeendcoordsnudge',edgeendcoordsnudge);

planedata.planeisthin = planeisthin;
planedata.planeseesplane = planeseesplane;
planedata.rearsideplane = rearsideplane;
planedata.canplaneobstruct = canplaneobstruct;
planedata.reflfactors = reflfactors;

elapsedtimeedgeo = etime(clock,t00);

if filehandlingparameters.saveeddatafile == 1
	eval(['save(''',desiredname,''',''planedata'',''edgedata'',''EDinputdatahash'',''elapsedtimeedgeo'');'])
end
