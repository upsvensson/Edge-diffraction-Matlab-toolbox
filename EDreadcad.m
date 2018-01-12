function [planedata,extraCATTdata] = EDreadcad(CADfile,planecornerstype,checkgeom)
% EDreadcad - Reads a file of type .CAD (made by e.g. CATT-Acoustic) and returns all the geometry data 
%  + extra geometrical data based on planes, in two structs.
% CAD-files v6,v7,v8 are read.
%
% Input parameters:
%  CADfile     	(optional) The input file, with or without the .CAD extension.
%					If this file is not specified, a file opening window will appear.
%  planecornerstype (optional) Could have the value 'zero' or 'circ'.Default: 'circ'.
%               	Affects the matrix planecorners, see below.
%  checkgeom		(optional) If this parameter is given the value 'check', then a few checks
%                   of the geometry consistency will be done: a check for duplicate corners
%                   for redundant corners and for corners that are connected to only one plane.
%					Only warnings are given. As default, no check is done.
%
% Output parameters:
%	planedata       Structwith fields:
%       corners         Matrix [ncorners,3] with the corner coordinates
%       planecorners    Matrix [planes,nmaxcornersperplane] with the corner numbers that make up each plane.
%                   Since planes can have different numbers of corners, the number of columns is the 
%                   maximum number of corners that any plane has. NB! The corner numbers are in the
%                   renumbered system, not in the CAD-file system.
%       planeabstypes   Matrix [nplanes,nmaxcharacters2] (sparse), with the absorber names in the CAD file. 
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
%   extraCATTdata   Struct, which presently has only one field: Sdirectivitynames, but could
%   easily get these fields as well (because they are extracted from the
%   CAD-file:
%           planenumbers    Vector [nplanes,1] with the plane numbers used in the CAD file.
%           cornercrossref  Vector [nmaxnumberinCADfile,1] (sparse), which gives the matlab-file corner
%                   numbers for each CAD-file corner number. For instance, if a CAD-file corner
%                   number 1200 is translated to corner number 72 in the matlab-file, then
%                   cornercrossref(1000) = 72.
%           Snumbers        Vector [nsources,1] of the source numbers in the CAD-file. NB! Only CAD-files
%                   v6 use source numbers; later versions use source names instead of numbers. 
%                   For CAD-file v7 and v8, Snumbers is empty.
%           Snames          Matrix [nsources,2] of the source names in the CAD-file. NB! Only CAD-files v7
%                   or later use source names (such as 'A0' etc). Older versions use source numbers, 
%                   in which case Snames is empty.
%           Sdirectivitynames   Matrix [nsources,nmaxcharacters3] (sparse) with the source directivity names 
%                   in the CAD file. 
%           Sdelays         Vector [nsources,1] of the source delays in the CAD-file, in milliseconds.
%           Scoords         Matrix [nsources,3] of the source coordinates.
%           Sdirections     Matrix [nsources,3] of the aim point coordinates for the sources.
%           Srotations      Vector [nsources,1] of the source rotations, in degrees.
%           Slevels         Matrix [nsources,6/8] of the source levels for six (CAD-files v6) or eight
%                   octave bands (CAD-files v7 or later).
%           Rnumbers        Vector [nreceivers,1] of the receiver numbers in the CAD-file.
%           Rcoords         Matrix [nreceivers,3] of the receiver coordinates.
%           CATTversionnumber   6,7 or 8
%           cornernumbers   Vector [ncorners,1] with the corner numbers used in the CAD-file. That is
%                   in the outputfile, the corners will effectively be renumbered from 1 to ncorners
%                   in the order they appeared in the CAD-file. 
%           planenames      Matrix [nplanes,nmaxcharacters1] (sparse), with the planenames in the CAD file.
%
% Uses the functions EDextrnums EDinfrontofplane
% 
% Peter Svensson (peter.svensson@ntnu.no) 12 Jan. 2018
%
% [planedata,extraCATTdata] = EDreadcad(CADfile,planecornerstype,checkgeom);

% 18 July 2009 Last previous fix
% 29 Oct. 2014 Fixed a bug for planes with indents, which gave the message
%              'Contact the developer'
% 3 Dec. 2014  Changed to the use of a struct "geometrydata" for the output data
% 31 March 2015 Added '' in the file saving, to handle file and directory
%               names with blanks
% 24 Nov. 2017  Copied, with very small changes (removed the SHOWTEXT), from ESIE2readcad
%               Added the fields cornerinfrontofplane and modeltype.
%               Removed the file saving.
% 28 Nov. 2017 Did some code improving as suggested by Matlab
% 3 Dec. 2017 Fixed a bug: if a plane definition had all zeros, then
% an infinite loop happened, around line 460.
% 12 Jan. 2018 Increased the bounding boxes a bit - doesn't hurt to make
% them a bit bigger.

if nargin == 0
	CADfile = '';
	planecornerstype = '';
	checkgeom = '';
elseif nargin == 1
	planecornerstype = '';
	checkgeom = '';
elseif nargin == 2
	checkgeom = '';
end

% geomacc is only used to make the bounding boxes a bit bigger than the
% corner coordinates. Was 1e-10 earlier.

geomacc = 1e-4;

%---------------------------------------------------------------
% If no CAD-file was specified, present a file opening window

if isempty(CADfile)
	[CADfile,filepath] = uigetfile('*.cad','Please select the CADfile');
    [filepath,~,~,~] = fileparts(filepath);
	if ~isstr(CADfile)
		return
	end
    [~,filestem,~,~] = fileparts(CADfile);

    CADfile = [[filepath,filesep],filestem,'.cad'];
else
    [filepath,filestem,fileext] = fileparts(CADfile);
    CADfile = [[filepath,filesep],filestem,fileext];
end


if exist(CADfile) ~= 2
	error(['ERROR: CAD-file: ',CADfile,' can not be found'])
end

% % %---------------------------------------------------------------
% % % If no output file was specified, construct an automatic file name
% % 
% % if isempty(outputfile)
% % 	outputfile = [[filepath,filesep],filestem,'_cadgeo.mat'];
% % end

%---------------------------------------------------------------
% Read in the entire file into the string B
% Find the line starts and ends

fid = fopen([CADfile],'r');
if fid == 0
	error(['ERROR: The CADfile ',CADfile,' could not be opened.'])
end
B = fread(fid,inf,'char').';
fclose(fid);

% Cut out sequential spaces to save some space.
% iv = find(B(1:length(B)-1)==32 & B(2:length(B))==32);
% B(iv) = [];
iv = B(1:length(B)-1)==32 & B(2:length(B))==32;
B(iv) = [];

% If the text file was generated on a Mac, each line has ASCII 13 at the end.
% On a PC, each line ends with ASCII 13, and the next line starts with ASCII 10.
% In Unix, each line ends with ASCII 10.
%
% Change this to a single ASCII 13

% New 17 Jul 09: use the builtin function regexprep instead:

% 1. Replace all ASCII 10 by ASCII 13

B = regexprep(char(B),'\n','\r');

% 2. Cut out sequential CRs (ASCII13) to save some space.
iv = B(1:length(B)-1)==13 & B(2:length(B))==13;
B(iv) = [];

% 3. Cut out sequences of CR, space, CR because they "survived" the
%    search for consecutive blanks.
iv = B(1:length(B)-1)==32 & B(2:length(B))==13;
B(iv) = [];
iv = B(1:length(B)-1)==13 & B(2:length(B))==13;
B(iv) = [];

% 4. Convert all text to lower case
B = lower(char(B));

% Some special characters might give negative ASCII values
% (Comment by PS 17Jul09: Will this ever happen??)

iv = find(B<0);
if ~isempty(iv)
	B(iv) = B(iv)+256;   
end

%---------------------------------------------------------------
% Look for the text 'CATT-Acoustic v' anywhere, and read the numerical
% value right after this text.

iv = regexp(B,'catt-acoustic v');

if isempty(iv)
	CATTversionnumber = [];
else

 	CATTversionnumber = str2double(char(B(iv(1)+15)));
end

%---------------------------------------------------------------
% Find the four sections of the file:
%    %CORNERS      %PLANES      %SOURCES   %RECEIVERS
% (or small letters)
% We assume that they come in this order

stringposCORNERS = regexp(B,'%corners');
if isempty(stringposCORNERS)
 	   error('ERROR: The .cad file must contain the word %CORNERS');    
end

stringposPLANES = regexp(B,'%planes');
if isempty(stringposPLANES)
 	   error('ERROR: The .cad file must contain the word %PLANES');    
end

stringposSOURCES = regexp(B,'%sources');

stringposRECEIVERS = regexp(B,'%receivers');

%---------------------------------------------------------------
% Look for the corners definitions between the lines containing
% %SOURCES and another keyword  (%PLANES or %SOURCES or %RECEIVERS)
% The lines should look like:
%	1   -0.2000   0.34   2.3


if stringposPLANES < stringposCORNERS
    error('ERROR: Sorry, but in the .cad file, the %CORNERS section must come before the %PLANES section')
end

C = textscan(B(stringposCORNERS:stringposPLANES-1),'%d%f%f%f','Headerlines',1);
cornernumbers = (C{1});
corners = [C{2} C{3} C{4}];
ncorners = length(cornernumbers);

%---------------------------------------------------------------
% Look for the planes definitions between the lines containing
% %PLANES and %SOURCES, or the end
% The lines should look like:
%  1 / planename /RIGID
%  1  4  3  2
%  2 / /RIGID

%  1 / planename /RIGID
%  1  4  3  2
%  2 / /RIGID

if ~isempty(stringposSOURCES)
    if stringposSOURCES < stringposPLANES
        error('ERROR: Sorry, but in the .cad file, the %PLANES section must come before the %SOURCES section')
    end    
    Str = B(stringposPLANES:stringposSOURCES-1);
else
    Str = B(stringposPLANES:end);    
end

% To make it easier with textscan, we put the two types of line onto
% a single line, with a new '/' inserted. Then all lines have the same
% format. 

% Find all CR and replace every second with a '/'
ivendlines = find(Str==13);
Str(ivendlines(2:2:end)) = '/';

C = textscan(Str,'%d%s%s%s','Headerlines',1,'Delimiter','/');

planenumbers = C{1};
Str2 = C{2};
Str3 = C{3};
Str = C{4};

nplanes = size(Str,1);
ncornersperplanevec = zeros(nplanes,1);
planecorners = zeros(nplanes,4);
planenames = zeros(nplanes,30);
planeabstypes = zeros(nplanes,6);

for ii = 1:nplanes
    listofplanes = EDextrnums(Str{ii});
    ncornersperplanevec(ii) = length(listofplanes);
    if ncornersperplanevec(ii) > size(planecorners,2)
       planecorners = [planecorners zeros(nplanes, ncornersperplanevec(ii)-size(planecorners,2))];
    end
    planecorners(ii,1:ncornersperplanevec(ii)) = listofplanes;

    txtstr2 = Str2{ii};
    if length(txtstr2) > size(planenames,2)
       planenames = [planenames zeros(nplanes,length(txtstr2)-size(planenames,2))]; 
    end
    planenames(ii,1:length(txtstr2)) = txtstr2;

    txtstr3 = Str3{ii};
    if length(txtstr3) > size(planeabstypes,2)
       planeabstypes = [planeabstypes zeros(nplanes,length(txtstr3)-size(planeabstypes,2))]; 
    end
    planeabstypes(ii,1:length(txtstr3)) = txtstr3;
end

if max(max(planecorners)) > max(cornernumbers)
    error('ERROR: One plane definition in the CAD-file used a higher corner number than was defined in the CORNERS section')
end

%---------------------------------------------------------------
% Look for the sources definitions between the lines containing
% %SOURCES and %RECEIVERS, or the end
% The lines should look like: (v6)
% 0 OMNI:SD0
%   0.0000000   0.0000000   1.0000000
%   0.0000000   0.0000000   0.0000000
%  85.0  88.0  91.0  94.0  97.0  100.0
%
% or: (v7)
% A0 OMNI:SD0
%   0.0000000   0.0000000   1.0000000
%   0.0000000   0.0000000   0.0000000  15.00000
%  85.0  88.0  91.0  94.0  97.0  100.0 : 103.0 106.0

% The lines should look like: (v6)
% 0 OMNI:SD0
%   0.0000000   0.0000000   1.0000000
%   0.0000000   0.0000000   0.0000000
%  85.0  88.0  91.0  94.0  97.0  100.0
%
% or: (v7)
% A0 OMNI:SD0   0
%   0.0000000   0.0000000   1.0000000
%   0.0000000   0.0000000   0.0000000  15.00000
%  85.0  88.0  91.0  94.0  97.0  100.0 : 103.0 106.0

if ~isempty(stringposSOURCES)

    if ~isempty(stringposRECEIVERS)
        if stringposRECEIVERS < stringposSOURCES
            error('ERROR: Sorry, but in the .cad file, the %SOURCES section must come before the %RECEIVERS section')
        end    
        Str = B(stringposSOURCES:stringposRECEIVERS-1);
    else
        Str = B(stringposSOURCES:end);    
    end

    % To make it easier with textscan, we put the four types of line onto
    % a single line, with new '/' inserted. Then all lines have the same
    % format. 
    % Find all CR and replace number 2,3,4, 6,7,8, etc with a '/'
    
    ivendlines = find(Str==13);
    ivkeep = (1:4:length(ivendlines));
    ivremove = (1:length(ivendlines));
    ivremove(ivkeep) = [];
    Str(ivendlines(ivremove)) = '/';

    C = textscan(Str,'%s%s%s%s','Headerlines',1,'Delimiter','/');

    nsources = size( C{1},1 );
    Scoords = zeros(nsources,3);
    Sdirections = zeros(nsources,3);
    Srotations = zeros(nsources,1);
    Slevels = zeros(nsources,8);
    Snumbers = zeros(nsources,1);
    Snames = zeros(nsources,4);
    Sdirectivitynames = ( zeros(nsources,30) );
    Sdelays = zeros(nsources,1);

    Str1 = C{1};
    Str2 = C{2};
    Str3 = C{3};
    Str4 = C{4};
    for ii = 1:nsources    
        Scoords(ii,:) = EDextrnums(Str2{ii});
        listofvalues3 = EDextrnums(Str3{ii});
        listofvalues4 = EDextrnums(Str4{ii});
        if ii==1
           if length(listofvalues3)==3
               CATTversionnumber = 6;
           elseif length(listofvalues3)==4
               CATTversionnumber = 7;
           else
               error('ERROR: The source data should have three or four values for the direction')
           end
        end

        Sdirections(ii,1:3) = listofvalues3(1:3);

        % The first field contains three items which are a mix of text and numerical values.
        % First we remove possible consecutive blanks.
        textstr = Str1{ii};
        iv = textstr(1:length(textstr)-1)==32 & textstr(2:length(textstr))==32;
        textstr(iv) = [];
        iv = find(textstr==32);


        if CATTversionnumber == 7
           Srotations(ii) = listofvalues3(4); 
           Snames(ii,1:length(textstr(1:iv(1)-1))) = textstr(1:iv(1)-1);
            Sdirectivitynames(ii,1:iv(2)-iv(1)-1) = textstr(iv(1)+1:iv(2)-1);
            Sdelays(ii) = str2double( textstr(iv(2)+1:end) );
        else
            Snumbers(ii) = str2double( textstr(1:iv(1)-1) );
            Sdirectivitynames(ii,1:length(textstr(iv(1)+1:end))) = textstr(iv(1)+1:end);
        end
        Slevels(ii,1:length(listofvalues4)) = listofvalues4;      
    end
    Snames = char(Snames);
    Sdirectivitynames = char(Sdirectivitynames);

else
    
    Scoords = [];
    Sdirections = [];
    Srotations = [];
    Slevels = [];
    Snumbers = [];
    Snames = [];
    Sdirectivitynames = [];
    Sdelays = [];    
    
end

%---------------------------------------------------------------
% Look for the receivers definitions after the line containing
% %RECEIVERS
% The lines should look like:
%	1   -0.2000   0.34   2.3



if ~isempty(stringposRECEIVERS)
    C = textscan(B(stringposRECEIVERS:end),'%d%f%f%f','Headerlines',1);
    Rnumbers = (C{1});
    Rcoords= [C{2} C{3} C{4}];
else
    Rnumbers = [];
    Rcoords = [];
end

%---------------------------------------------------------------
% Change the numbering in the CAD file to a contiguous one

% Corner numbers
% First we make an inelegant crossreference vector giving the
% resulting corner number in the position given by the CAD-file
% corner number.

maxnumb = max(cornernumbers);
cornercrossref = sparse(zeros(1,maxnumb));
cornercrossref(cornernumbers) = (1:ncorners);

[nplanes,ncornersperplane] = size(planecorners);
iv = find(planecorners~=0);
planecorners(iv) = full(cornercrossref(planecorners(iv)));

%---------------------------------------------------------------
% Go through all planes. If there is a plane definition including
% zeros, and planecornerstype == 'circ', expand it repeating the
% same corner order again.

if isempty(planecornerstype)
	planecornerstype = 'circ';
else
	planecornerstype = char(lower(planecornerstype(1)));
	if planecornerstype(1) == 'z'
		planecornerstype = 'zero';
    else
		planecornerstype = 'circ';
	end
end

% In the next loop, before the && ~isempty(iv) was inserted, an infinite
% loop was generated if a plane definition had all zeros.

if strcmp(planecornerstype,'circ') == 1
	for ii = 1:nplanes
		iv = find( planecorners(ii,:) ~= 0);
		ncornersatplane = length(iv);
		if (ncornersatplane ~= ncornersperplane) && ~isempty(iv)
			pattern = planecorners(ii,iv);
			nrepeatings = ceil(ncornersperplane/ncornersatplane);
			for jj = 1:nrepeatings-1
				pattern = [pattern planecorners(ii,iv)];
			end
			planecorners(ii,:) = pattern(1:ncornersperplane);
		end
	end
end

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
       error(['ERROR: Plane ',int2str(planenumbers(ii)),' (plane numbering as in the CAD file) seems to be twisted.'])        
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
        error(['ERROR: Normal vectors for plane ',int2str(planenumbers(ii)),' (in the CAD file, = ',int2str(ii),' in the ESIE2 file), get different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    elseif any(nvecdiff>1e-8)
        nvecdiff
        disp(['WARNING: Normal vectors for plane ',int2str(planenumbers(ii)),' (in the CAD file, = ',int2str(ii),' in the ESIE2 file), get somewhat different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    end
    
	if ncorners > 5 && abs(sum(nvecsigns)) <= 1
		disp(['WARNING for plane number ',int2str(planenumbers(ii)),' in the CAD-file'])
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
% Check the geometry a bit

if isempty(checkgeom)
	checkgeom = 0;
else
	if char(checkgeom(1)) == 'c'
		checkgeom = 1;
    else
		checkgeom = 0;
	end
end

if checkgeom

	badproblem = 0;
	disp(' ')
	disp('   Checking if there are corners with identical coordinates')
	for ii = 1:ncorners-1
		othercorners = [ii+1:ncorners];
		iv = find((corners(othercorners,1)==corners(ii,1)) & (corners(othercorners,2)==corners(ii,2)) & (corners(othercorners,3)==corners(ii,3)));
		if ~isempty(iv)
			disp(['      WARNING: corner ',int2str(ii),' and ',int2str(othercorners(iv(1))),' have the same coordinates'])
			disp(['      This should be fixed in the CAD-file or in the CATT-model'])
			badproblem = 1;
		end
	end

	disp('   Checking if there are unused corners')
	iv = find(planesatcorners(:,1)==0);
	if ~isempty(iv)
		disp('      WARNING: The following corners in the CAD-file are not used:')
		printvec = int2str(cornernumbers(iv(1)));
		for iiprint = 2:length(iv)
			printvec = [printvec,' ',int2str(cornernumbers(iv(iiprint)))];
		end
		disp(['            ',printvec])
		disp(['      This is not a big problem'])
		planesatcorners(iv,:) = [];
	end

	disp('   Checking if any corners seem redundant')
	iv = find(planesatcorners(:,2)==0);
	if ~isempty(iv)
		disp('      WARNING: The following corners in the CAD-file belong to only one plane:')
		printvec = int2str(cornernumbers(iv(1)));
		for iiprint = 2:length(iv)
			printvec = [printvec,' ',int2str(cornernumbers(iv(iiprint)))];
		end
		disp(['            ',printvec])
		disp(['      This should be fixed in the CAD-file or in the CATT-model'])
		badproblem = 1;
	end

	if badproblem == 1
		disp(' ');
		error(['ERROR: Problems with the geometry defined in the CAD-file: ',CADfile])
	end

end

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
% Save the relevant variables in the output file

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

% Varlist = [' planedata cornernumbers planenames '];
% Varlist = [Varlist,' planenumbers cornercrossref Snumbers Snames Sdirectivitynames Sdelays'];
% Varlist = [Varlist,' Scoords Sdirections Srotations Slevels Rnumbers Rcoords CATTversionnumber'];
% Varlist = [Varlist,' planenvecs   planesatcorners nplanespercorners cornerinfrontofplane modeltype'];

planenames = sparse(planenames+1-1);
planeabstypes = sparse(planeabstypes+1-1);
Sdirectivitynames = sparse(Sdirectivitynames+1-1);

extraCATTdata = struct('Sdirectivitynames',Sdirectivitynames);
planedata = struct('corners',corners,'planecorners',planecorners,...
    'planeabstypes',planeabstypes,...
    'planeeqs',planeeqs,'ncornersperplanevec',ncornersperplanevec,...
    'minvals',minvals,'maxvals',maxvals,...
    'planehasindents',planehasindents,'indentingcorners',indentingcorners,...
    'cornerinfrontofplane',cornerinfrontofplane,'modeltype',modeltype);

