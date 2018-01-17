function [corners,planecorners,ncornersperplanevec,radius] = EDmakegeo_cylinder(radius,width,numberofcorners,intextgeom,fullcirclefract,angleoffset,ztranslation)
% EDmakegeo_cylinder generates the geometry for a polyhedral cylinder.
% The end caps of the cylinder will be parallel to the z=0 plane,
% symmetrically placed. Optionally, the whole cylinder can be translated in
% the z-direction.
% 
% The polygonal shape will be scaled so that regardless of numbers of
% edges, the volume of the polygonal cylinder will be the same as a truly
% circular cylinder with the radius given as input parameter.
%
% Input parameters:
%	radius			The radius of the equivalent circular cylinder
%	width			The length of the cylinder
%	numberofcorners	The number of corners approximating the circular cross-section
%	intextgeom		'int' or 'ext' defining if an interior or exterior geometry
%					should be constructed
%	fullcirclefract	1 or 0.5, indicating whether a full or half cylinder should be created.
%   angleoffset		(optional) 1 or 0, indicating whether or not the first corner
%   				should be placed in y = 0 (angleoffset = 0) or shifted
%					half an angle step (angleoffset = 1).
%   ztranslation    (optional) A value which will be added to all the
%                   z-coordinates.
%
% Output parameters:
%   corners         Matrix, [2*numberofcorners,3], with the corner
%                   coordinates
%   planecorners    Matrix, [numberofcorners+2,numberofcorners], with the
%                   plane corners. The first numberofcorners rows have the
%                   planes of the cylindrical surface, and the two last
%                   rows have the flat circular endsurfaces.
%   ncornersperplanevec     A vector which gives the number of corners per
%                   plane
%	radius			The actual radius of the corners; this will be different
%					from the (desired) input parameter 'radius', since the
%					polygonal approximation of the circle is scaled so that
%					it gets the same cross-section area as the real circle.
%
% Peter Svensson (peter.svensson@ntnu.no)   17 Jan 2018
%
% [corners,planecorners,ncornersperplanevec,radius] =
% EDmakegeo_cylinder(radius,width,numberofcorners,intextgeom,fullcirclefract,angleoffset,ztranslation);

% 9 September 2006  Stable version
% 16 October 2012   Added ztranslation optional input parameter
% 17 Jan 2018 Stripped off the CAD-file writing, and saved to the
% EDtoolbox, renaming from GEOmakewheelcad

if nargin < 7
    ztranslation = 0;
end
if nargin < 6
   angleoffset = 0; 
end
if nargin < 5
   fullcirclefract = 1;
end
if nargin < 4
    intextgeom = 'ext';
end

if angleoffset ~= 0 && angleoffset ~= 1
   error('ERROR: only values 0 and 1 are allowed for angleoffset') 
end

intextgeom = lower(intextgeom(1));
if intextgeom ~= 'e' && intextgeom ~= 'i'
	error(['ERROR: The parameter intextgeom must have the value ''int'' or ''ext''.'])
end
if fullcirclefract ~= 1 && fullcirclefract ~= 0.5
	error(['ERROR: The parameter fullcirclefract must have the value 1 or 0.5'])
end
xtranslation = [0 0 ztranslation];

constantarea = 1;

if constantarea == 1
	radius = radius*sqrt(2*pi/numberofcorners/sin(2*pi/numberofcorners));
end

%------------------------------------------------------------------------------
% Calculate the corner coordinates

corners = zeros(numberofcorners*2,3);

fivec = [0:numberofcorners-1].'*2*pi/numberofcorners;
if angleoffset == 1
   fivec = fivec + 0.5*2*pi/(numberofcorners); 
end
onesvec = ones(size(fivec));
corners(1:numberofcorners,1) = radius*cos(fivec) + xtranslation(1)*onesvec;
corners(1:numberofcorners,2) = radius*sin(fivec) + xtranslation(2)*onesvec;
corners(numberofcorners+1:2*numberofcorners,1) = radius*cos(fivec) + xtranslation(1)*onesvec;
corners(numberofcorners+1:2*numberofcorners,2) = radius*sin(fivec) + xtranslation(2)*onesvec;
corners(1:numberofcorners,3) = width/2*ones(numberofcorners,1) + xtranslation(3)*onesvec;
corners(numberofcorners+1:2*numberofcorners,3) = -width/2*ones(numberofcorners,1) + xtranslation(3)*onesvec;;

%------------------------------------------------------------------------------
% Define the planes

planecorners = zeros(numberofcorners+2,numberofcorners);
for ii = 1:numberofcorners-1
	planecorners(ii,1:4) = [ii+1 ii numberofcorners+ii numberofcorners+ii+1];
end
planecorners(numberofcorners,1:4) = [1 numberofcorners 2*numberofcorners numberofcorners+1];

% sideplanecorners = [1:numberofcorners;2*numberofcorners:-1:numberofcorners+1];
planecorners(numberofcorners+1,:) =  [1:numberofcorners];
planecorners(numberofcorners+2,:) =  [2*numberofcorners:-1:numberofcorners+1];

ncornersperplanevec = 4*ones(numberofcorners+2,1);
ncornersperplanevec(end-1:end) = numberofcorners;

% surfplanecorners = zeros(numberofcorners,4);
% for ii = 1:numberofcorners-1
% 	surfplanecorners(ii,:) = [ii+1 ii numberofcorners+ii numberofcorners+ii+1];
% end
% surfplanecorners(numberofcorners,:) = [1 numberofcorners 2*numberofcorners numberofcorners+1];


% % % % %------------------------------------------------------------------------------
% % % % % Create a CAD file
% % % % 
% % % % comptype = computer;
% % % % comptype = comptype(1:2);
% % % % if comptype == 'MA',
% % % % 	Lineending = 13;
% % % % elseif comptype(1) == 'S',
% % % % 	Lineending = 10;
% % % % else,
% % % % 	if comptype == 'PC',
% % % % 		Lineending = [13,10];
% % % % 	else,
% % % % 		error('ERROR: Not implemented for this computer type yet')	
% % % % 	end
% % % % end
% % % % 
% % % % fid = fopen([CADfilepath,CADfile],'w');
% % % % if fid == -1,
% % % % 	disp('This file is not possible to open - check that it isn''t opened by any program!')
% % % % 	return
% % % % end
% % % % 
% % % % %------------------------------------------------------------------------------
% % % % 
% % % % fwrite(fid,['% ',CADfile,' created on ',date,' by GEOmakewheelcad',Lineending],'char');
% % % % fwrite(fid,[' ',Lineending],'char');
% % % % 
% % % % %------------------------------------------------------------------------------
% % % % % Write the corners
% % % % 
% % % % fwrite(fid,['%CORNERS',Lineending],'char');
% % % % fwrite(fid,[' ',Lineending],'char');
% % % % if fullcirclefract == 1,
% % % % 	for ii = 1: 2*numberofcorners + 2*numberofcorners*reflplane,
% % % % 		fwrite(fid,	[int2str(ii),' ',num2str(corners(ii,1),12),' ',num2str(corners(ii,2),12),' ',num2str(corners(ii,3),12),Lineending],'char');
% % % % 	end
% % % % elseif fullcirclefract == 0.5,
% % % % 	for ii = [[1:numberofcorners/2+1] [numberofcorners+1:numberofcorners+numberofcorners/2+1]],
% % % % 		fwrite(fid,	[int2str(ii),' ',num2str(corners(ii,1),12),' ',num2str(corners(ii,2),12),' ',num2str(corners(ii,3),12),Lineending],'char');
% % % % 	end
% % % % end
% % % % %------------------------------------------------------------------------------
% % % % % Write the planes
% % % % 
% % % % fwrite(fid,[' ',Lineending],'char');
% % % % fwrite(fid,['%PLANES',Lineending],'char');
% % % % fwrite(fid,[' ',Lineending],'char');
% % % % covec1 = int2str(sideplanecorners(1,1));
% % % % if fullcirclefract == 1,
% % % % 	covec2 = int2str(sideplanecorners(2,1));
% % % % elseif fullcirclefract == 0.5,
% % % % 	covec2 = int2str(sideplanecorners(2,numberofcorners/2));
% % % % end
% % % % if intextgeom == 'e',
% % % % 	if fullcirclefract == 1,
% % % % 		for ii = 2:numberofcorners,
% % % % 			covec1 = [covec1,' ',int2str(sideplanecorners(1,ii))];
% % % % 			covec2 = [covec2,' ',int2str(sideplanecorners(2,ii))];
% % % % 		end
% % % % 	elseif fullcirclefract == 0.5,
% % % % 		for ii = 2:numberofcorners/2+1,
% % % % 			covec1 = [covec1,' ',int2str(sideplanecorners(1,ii))];
% % % % 			covec2 = [covec2,' ',int2str(sideplanecorners(2,ii+numberofcorners/2-1))];
% % % % 		end	
% % % % 	end
% % % % else,
% % % % 	if fullcirclefract == 1,
% % % % 		for ii = numberofcorners:-1:2,
% % % % 			covec1 = [covec1,' ',int2str(sideplanecorners(1,ii))];
% % % % 			covec2 = [covec2,' ',int2str(sideplanecorners(2,ii))];
% % % % 		end
% % % % 	elseif fullcirclefract == 0.5,
% % % % 		for ii = numberofcorners/2+1:-1:2,
% % % % 			covec1 = [covec1,' ',int2str(sideplanecorners(1,ii))];
% % % % 			covec2 = [covec2,' ',int2str(sideplanecorners(2,ii+numberofcorners/2-1))];
% % % % 		end	
% % % % 	end
% % % % end
% % % % if endcapsabsorbing == 0,
% % % %     fwrite(fid,	['1 / /RIGID',Lineending],'char');
% % % %     fwrite(fid,	[covec1,Lineending],'char');
% % % %     fwrite(fid,[' ',Lineending],'char');
% % % %     fwrite(fid,	['2 / /RIGID',Lineending],'char');
% % % %     fwrite(fid,	[covec2,Lineending],'char');
% % % %     fwrite(fid,[' ',Lineending],'char');
% % % % else,
% % % %     fwrite(fid,	['1 / /TOTABS',Lineending],'char');
% % % %     fwrite(fid,	[covec1,Lineending],'char');
% % % %     fwrite(fid,[' ',Lineending],'char');
% % % %     fwrite(fid,	['2 / /TOTABS',Lineending],'char');
% % % %     fwrite(fid,	[covec2,Lineending],'char');
% % % %     fwrite(fid,[' ',Lineending],'char');    
% % % % end
% % % % if fullcirclefract == 1,
% % % % 	lastcorner = numberofcorners;
% % % % elseif fullcirclefract == 0.5,
% % % % 	lastcorner = numberofcorners/2;
% % % % end
% % % % for ii = 1:lastcorner,
% % % % 	if intextgeom == 'e',
% % % % 		covec = [int2str(surfplanecorners(ii,1)),' ',int2str(surfplanecorners(ii,2)),' ',int2str(surfplanecorners(ii,3)),' ',int2str(surfplanecorners(ii,4))];
% % % % 	else,
% % % % 		covec = [int2str(surfplanecorners(ii,4)),' ',int2str(surfplanecorners(ii,3)),' ',int2str(surfplanecorners(ii,2)),' ',int2str(surfplanecorners(ii,1))];	
% % % % 	end
% % % % 	fwrite(fid,	[int2str(ii+2),' / /RIGID',Lineending],'char');
% % % % 	fwrite(fid,[covec,Lineending],'char');
% % % % 	fwrite(fid,[' ',Lineending],'char');
% % % % end
% % % % 
% % % % if reflplane == 1,
% % % % 	covec1 = int2str(sideplanecorners(3,1));
% % % % 	covec2 = int2str(sideplanecorners(4,1));
% % % % 	for ii = 2:numberofcorners,
% % % % 		covec1 = [covec1,' ',int2str(sideplanecorners(3,ii))];
% % % % 		covec2 = [covec2,' ',int2str(sideplanecorners(4,ii))];
% % % % 	end
% % % % 	II1 = int2str(3+numberofcorners);
% % % % 	II2 = int2str(4+numberofcorners);
% % % % 	fwrite(fid,	[II1,' / /RIGID',Lineending],'char');
% % % % 	fwrite(fid,	[covec1,Lineending],'char');
% % % % 	fwrite(fid,[' ',Lineending],'char');
% % % % 	fwrite(fid,	[II2,' / /RIGID',Lineending],'char');
% % % % 	fwrite(fid,	[covec2,Lineending],'char');
% % % % 	fwrite(fid,[' ',Lineending],'char');	
% % % % 	for ii = 1:numberofcorners,
% % % % 		covec = [int2str(surfplanecorners(ii+nset2,1)),' ',int2str(surfplanecorners(ii+nset2,2)),' ',int2str(surfplanecorners(ii+nset2,3)),' ',int2str(surfplanecorners(ii+nset2,4))];
% % % % 		fwrite(fid,	[int2str(ii+4+numberofcorners),' / /RIGID',Lineending],'char');
% % % % 		fwrite(fid,[covec,Lineending],'char');
% % % % 		fwrite(fid,[' ',Lineending],'char');
% % % % 	end
% % % % end
% % % % 
% % % % fwrite(fid,['%EOF',Lineending],'char');
% % % % 
% % % % fclose(fid);
