function [lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom] = ...
    EDfindISEStree(planedata,edgedata,edgetoedgedata,S,specorder,difforder,visplanesfromS,vispartedgesfromS,nedgesubs,showtext)
% EDfindISEStree - Constructs a list ("an ISES tree") of all possible specular-diffraction combinations.
% EDfindISEStree constructs a list ("an ISES tree") of all possible specular-diffraction combinations that
% are possible for the given source. This list is based on the visibility
% matrices from source to planes/edges and between planes/edges. Visibility
% checks are employed wherever possibble. NB! Visibility checks that depend
% on the receiver position are not used at this stage.
% Planes are numbered 1-nplanes, and edges have numbers
% [nplanes+1:nplanes+nedges].
%
% Input parameters:
%   planedata, edgedata,edgetoedgedata   structs
% 	S
% 	specorder                   The max. desired refl. order
% 	difforder                   The max. desired diffr. order
% 	visplanesfromS, vispartedgesfromS       From the Sdata struct
% 	nedgesubs                   The number of subdivisions of each edge for
%                               a simplified edge visibility test
%   showtext
% GLOBAL parameters
%   SUPPRESSFILES               If this global parameter has the value 1,
%                               then the results from this function will be returned in a struct
%                               rather than in a file.
%
% Output parameters:
% 	lengthNspecmatrix           A list, [1,specorder], with the number of
% 	                            entries in each column of IVNSPECMATRIX.
% 	lengthNdiffmatrix           A list, [1,specorder], with the number of
% 	                            entries in each column of IVNDIFFMATRIX.
% 	singlediffcol               A list, [nfirstorderdiff,1], containing the column number
%                               that the first-order diffracted component can be found in (in
% 	                            POTENTIALISES).
% 	startindicessinglediff      A list, [specorder,1], where the value in position N contains
%                               which row in IVNDIFFMATRIX(:,1) that has the first combination
%                               of order N and with one diffraction component.
% 	endindicessinglediff        A list, [specorder,1], where the value in position N contains
%                               which row in IVNDIFFMATRIX(:,1) that has the last combination
%                               of order N and with one diffraction 
% 	ndecimaldivider             See below.
% 	PointertoIRcombs            A sparse list which for position N contains a row number
%                               (in POTENTIALISES) that has a combination
%                               that ends with a specular reflection in
%                               plane N, after a diffraction. A row with one of the specular
%                               reflections M,N after a diffraction can be
%                               found in position M*ndecimaldivider+N etc
%                               for higher orders.
% 	IRoriginsfrom               A list, [npossiblecombs,1], where the value in position N
%                               states which other row number (in POTENTIALISES) that the
%                               image receiver in row N origins from.
% Global output data:
% 	POTENTIALISES               A matrix, [npossiblecombs,specorder], of all the 
% 	                            possible specular-diffraction combos that are possible for the
% 	                            given source. The combos are denoted by
% 	                            plane numbers and edge numbers, but edge
% 	                            numbers have a constant number (=nplanes)
% 	                            added to them, i.e., if there are 20 planes
% 	                            in the model, edge number 1 will be denoted
% 	                            by number 21.
% 	ORIGINSFROM                 A list, [npossiblecombs,1], where the value
% 	                            in row N states which other row in ISCOORDS that
% 	                            the image source in row N originates from.
% 	ISCOORDS                    A matrix, [npossiblecombs,3], containing
% 	                            the image source coordinates, where this is applicable,
%                               i.e., where the combo starts with a specular reflection. 
% 	ISESVISIBILITY              A list, [npossiblecombs,1], containing the visibility of
%                               the first edge in a sequence, seen from the source. 
% 	IVNSPECMATRIX               A matrix, [nspeccombs,specorder], where each column contains
%                               the row numbers in POTENTIALISES that contain combos with
%                               1 or 2 or... or "specorder" purely specular reflections.
% 	REFLORDER                   A list, [npossiblecombs,1], containing the
% 	                            order of reflection for each row of POTENTIALISES
% 	IVNDIFFMATRIX               A matrix, [ndiffcombs,specorder], where each column contains
%                               the row numbers in POTENTIALISES that contain combos with
%                               1 or 2 or... or "specorder" diffractions.
%
% Uses functions  EDfindis EDgetedgepoints EDchkISvisible EDcheckobstrpaths
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 16 March 2021
%
% [POTENTIALISES,ORIGINSFROM.ISCOORDS ,ISESVISIBILITY,IVNSPECMATRIX,lengthNspecmatrix,IVNDIFFMATRIX,lengthNdiffmatrix, ...
%    singlediffcol,REFLORDER,startindicessinglediff,endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom] = ...
%    EDfindISEStree(planedata,edgedata,edgetoedgedata,S,specorder,difforder,visplanesfroms,vispartedgesfroms,nedgesubs);

% 10.2.2010 Functioning version
% 29.11.2015 Found a bug related to "exist" for struct variables
% 29 Apr 2016 Made a small change of the PointertoIRcombs handling: expand
%             its size in steps.
% 6 Apr 2017 When difforder = 1, the edgeseesedge should not be included in
% the PLANEseesPLANE matrix
% 19 Oct. 2017 Had forgotten one struct
% 27 Nov. 2017 Copied from ESEI2toolbox. Removed the doconetrace-part.
%              Removed the input parameter isou.
% 28 Nov. 2017 Cleaned up code a bit. Removed global parameters from the
% output list.
% 16 Mar 2021 Adapted to change in EDchkISvisible.

global POTENTIALISES ISCOORDS IVNDIFFMATRIX
global IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER

if nargin < 10
    showtext = 0; 
    if nargin < 9
        nedgesubs = 2; 
    end
end

GEOMACC = 1e-10;

% planedata.corners = size(planedata.corners,1);
nplanes = length(planedata.planeisthin);
nedges = length(edgedata.closwedangvec);
nhighval = nplanes + nedges;

onesvec1 = ones(1,nedgesubs);
obstructtestneeded = (sum(planedata.canplaneobstruct)~=0);
% onesvec3 = ones(1,3);

planenvecs = planedata.planeeqs(:,1:3);
        
%--------------------------------------------------------------------------
% We modify edgeseesplane because in ESIE2edgeo, edgeseesplane was set to 1
% for totabs planes (this was necessary for ESIE2srgeo).

% Set all edges belonging to planes with reflfactors = 0 to edgeseesplane = 0.

% % listofabsorbingplanes = find(planedata.reflfactors == 0);
% % 
% % if ~isempty(listofabsorbingplanes)
% % 	edgeseesplane(listofabsorbingplanes,:) = zeros(length(listofabsorbingplanes),nedges);
% % end

%--------------------------------------------------------------------------
% Combine planeseesplane with edgeseesplane
% and visplanesfroms with vispartedgesfroms
% so that edges appear as planes

if difforder > 0
    visPLANESfroms = [visplanesfromS;2*sign(double(vispartedgesfromS))];
    % source sees planes that have visplanesfroms:
    %   2 => in front of reflective planes
    %   4 => inside plane which is reflective
    %   old version:  sousees = (visPLANESfroms==1 | visPLANESfroms==-2);
    sousees = (visPLANESfroms==2 | visPLANESfroms==4);
    
    % Bug found 29.11.2015: previous code used "exist" on the struct field
    % but exist can not be used like that 
    
    if isfield(edgetoedgedata,'edgeseespartialedge') ~= 1
        edgeseesedge = zeros(nedges,nedges,'int8');    
    else
                
        if ~isempty(edgetoedgedata.edgeseespartialedge)

            edgeseesedge = int8(full(abs(double(edgetoedgedata.edgeseespartialedge))>0));
            
        else
            edgeseesedge = zeros(nedges,nedges,'int8');    
        end
    end

    if difforder > 1
        PLANEseesPLANE = [[planedata.planeseesplane edgedata.edgeseesplane];[edgedata.edgeseesplane.' edgeseesedge]];
    else
        PLANEseesPLANE = [[planedata.planeseesplane edgedata.edgeseesplane];[edgedata.edgeseesplane.' edgeseesedge*0]];        
    end
    clear edgeseesplane edgeseesedge
else
    visPLANESfroms = visplanesfromS;
    % 19 Oct. 2017 har forgotten the planedata struct in the line below
    PLANEseesPLANE = planedata.planeseesplane;
%   old version: sousees = (visPLANESfroms==1 | visPLANESfroms==-2);
    sousees = (visPLANESfroms==2 | visPLANESfroms==4);

    clear visplanesfromS
end

% PLANEseesPLANE = int8(full(PLANEseesPLANE));
startindices = zeros(specorder,1);
endindices = zeros(specorder,1);

%--------------------------------------------------------------------------
% Construct a list of which planes a plane sees so that a search isn't
% needed later.

nPLANES = size(PLANEseesPLANE,1);
listofvisPLANES = zeros(nPLANES,nPLANES-1);    
nvisPLANES = zeros(nPLANES,1);
for ii = 1:nPLANES
    visPLANES = find(PLANEseesPLANE(ii,:)==1);
    nvisPLANES(ii) = length(visPLANES);
    listofvisPLANES(ii,1:nvisPLANES(ii)) = visPLANES;
end

%##################################################################
%##################################################################
%##################################################################
%
%         First order
%
%------------------------------------------------------------------

startindices(1) = 1;

possiblefirsts  =  find( sousees );
npossiblefirsts = length(possiblefirsts);

if nhighval < 255
    POTENTIALISES = uint8( [[possiblefirsts ] zeros(npossiblefirsts,specorder-1)] );
else
    POTENTIALISES = uint16( [[possiblefirsts ] zeros(npossiblefirsts,specorder-1)] );
end

ORIGINSFROM = zeros(npossiblefirsts,1,'uint32');
if nedgesubs < 8 
    ISESVISIBILITY = zeros(npossiblefirsts,1,'uint8');
elseif nedgesubs < 16
    ISESVISIBILITY = zeros(npossiblefirsts,1,'uint16');
else    
    ISESVISIBILITY = zeros(npossiblefirsts,1,'uint32');
end
iv = find(possiblefirsts>nplanes);
if ~isempty(iv)
    ISESVISIBILITY(iv) = vispartedgesfromS(possiblefirsts(iv)-nplanes);    
end

endindices(1) = npossiblefirsts;

% Compute the IS coords

ivdiff = (startindices(1):endindices(1));
ivspec = find(POTENTIALISES(ivdiff,1)<=nplanes);
ivdiff(ivspec) = [];

ISCOORDS = zeros(npossiblefirsts,3);
%ISCOORDS(ivspec,:) = ESIE2findis(S,POTENTIALISES(ivspec,1),planedata.planeeqs);
ISCOORDS(ivspec,:) = EDfindis(S,POTENTIALISES(ivspec,1),planedata.planeeqs);

%##################################################################
%##################################################################
%##################################################################
%
% Second order
%
%---------------------------------------------------------------------------------------------

startindices(2) = startindices(1) + length(possiblefirsts);

if specorder > 1
    if showtext >= 3
        disp(['     Order number 2'])
    end
    
    if nhighval < 255
        startplanes = uint8(possiblefirsts);
    else
        startplanes = uint16(possiblefirsts);
    end
        
	addtocomb = startplanes;
	nadds = size(addtocomb,1);
    if nadds > 0
        
       	maxnumberofvispla = max(sum(PLANEseesPLANE(:,startplanes)==1));
	
        addtocomb = (reshape(addtocomb(:,ones(1,maxnumberofvispla)).',length(addtocomb)*maxnumberofvispla,1));
		addtocomb = [addtocomb zeros(size(addtocomb))];
        
        addtoISESVISIBILITY = reshape(ISESVISIBILITY(:,ones(1,maxnumberofvispla)).',nadds*maxnumberofvispla,1);
        
        startpos = [[1:maxnumberofvispla:length(addtocomb)]  length(addtocomb)+1  ].';
		
		for ii = 1:length(startpos)-1
%			possibleplanes = find( (PLANEseesPLANE(:,addtocomb(startpos(ii)))==1) );
            if nvisPLANES(addtocomb(startpos(ii))) > 0
    			possibleplanes = listofvisPLANES(addtocomb(startpos(ii)),1:nvisPLANES(addtocomb(startpos(ii)))).';
			    nposs = length(possibleplanes);
			    addtocomb( startpos(ii):startpos(ii)+nposs-1,2) = possibleplanes;
            end
            
		end
	
        addtooriginsfrom = uint32([1:startindices(2)-1].');
		addtooriginsfrom = addtooriginsfrom(:,ones(1,maxnumberofvispla));
		addtooriginsfrom = reshape(addtooriginsfrom.',maxnumberofvispla*(startindices(2)-1),1);
		
		cutout = find(addtocomb(:,2)==0);
		addtocomb(cutout,:) = [];
		addtooriginsfrom(cutout) = [];
        addtoISESVISIBILITY(cutout) = [];
% %         beamangle(cutout) = [];
% %         beamdirection1(cutout) = [];
        clear cutout
        nadds = size(addtocomb,1);
		    
        if nadds > 0

            %--------------------------------------------------------------------------
			% Try to get rid of some combinations that we know can not be propagated
                        
			% Find those combinations of specular-specular where the IS is behind the second reflecting plane
			% because they can then never ever be propagated, or valid
			
            if difforder > 0
                ivss = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)<=nplanes));
            else
                ivss = uint32([1:nadds].');
            end

			imsoucheck = dot((ISCOORDS(addtooriginsfrom(ivss),1:3) - planedata.corners(planedata.planecorners(addtocomb(ivss,2),1),:)).',planenvecs(addtocomb(ivss,2),:).').';
			imsounotOK = find(imsoucheck < GEOMACC);
			addtocomb(ivss(imsounotOK),:) = [];
			addtooriginsfrom(ivss(imsounotOK),:) = [];
            addtoISESVISIBILITY(ivss(imsounotOK)) = [];
% % %             beamangle(ivss(imsounotOK)) = [];
% % %             beamdirection(ivss(imsounotOK)) = [];
            % %             ivss(imsounotOK) = [];

            nadds = size(addtocomb,1);            

			% Combinations of diffraction-specular don't need to be checked because
			% 'edgeseesplane' should have taken care of this.
			%%%%ivds = find(addtocomb(:,1)>nplanes  & addtocomb(:,2)<=nplanes);
			
			%--------------------------------------------------------------------------
			% Combinations of specular-diffraction: check if the IS is behind both of
			% the two planes that the reflecting edge connects to.
            %
            % Bug found 080711: The visibility test must be split into two, depending on the wedge angle!! 
            % If the open wedge angle > pi, then it is correct to check if the IS behind both of the planes. 
            % If the open wedge angle < pi, then it should be checked if the IS behind one of the planes!!             
            
            if difforder > 0 && nadds > 0

                % First we check the sd combinations where the open wedge angle > pi
                % For those cases, we can exclude combinations where the IS
                % is behind both planes
                
  				ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                ivsd1 = ivsd( edgedata.closwedangvec( double(addtocomb(ivsd,2))-nplanes ) < pi );

                if ~isempty(ivsd1)
                    imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)    ,:).').';
                    imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)    ,:).').';
                    GEOMACC = 1e-10;
                    imsounotOK = find(imsoucheck1 < GEOMACC & imsoucheck2 < GEOMACC);
                    clear imsoucheck1 imsoucheck2

                    addtocomb(ivsd1(imsounotOK),:) = [];
                    addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                    addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                end
                
                % Second we check the sd combinations where the open wedge angle < pi
                % For those cases, we can exclude combinations where the IS
                % is behind one of the two planes
                
  				ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                if ~isempty(ivsd)
                    ivsd1 = ivsd( edgedata.closwedangvec( double(addtocomb(ivsd,2))-nplanes ) > pi );

                    if ~isempty(ivsd1)
                        imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)    ,:).').';
                        imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)    ,:).').';
                        GEOMACC = 1e-10;
                        imsounotOK = find(imsoucheck1 < GEOMACC | imsoucheck2 < GEOMACC);
                        clear imsoucheck1 imsoucheck2

                        addtocomb(ivsd1(imsounotOK),:) = [];
                        addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                        addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                    end
                end
                
				nadds = size(addtocomb,1);
                
                ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                
                if ~isempty(ivsd)
		
                    %----------------------------------------------------------------------
                    % Combinations of specular-diffraction that are not visible at
                    % all should be removed and then they are not propagated either.
				
                    if nedges < 255
                        possibleedges = uint8(double(addtocomb(ivsd,2)) - nplanes);
                    else
                        possibleedges = uint16(double(addtocomb(ivsd,2)) - nplanes);
                    end
                     
                    possiblecombs = addtocomb(ivsd,1);
                    reftoISCOORDS = addtooriginsfrom(ivsd);
				
                    % Expand to take the edge subdivisions into account
                    
                    nposs = length(ivsd);
                    nposs = nposs*nedgesubs;        % We treat all the little edge subdivisions as separate edges
				
                    expandedposscombs = possiblecombs(:,onesvec1);
                    clear possiblecombs
                    expandedposscombs = reshape(expandedposscombs.',nposs,1);			
                    
                    expandedreftoISCOORDS = reftoISCOORDS(:,onesvec1);
                    expandedreftoISCOORDS = reshape(expandedreftoISCOORDS.',nposs,1);
                    expandedpossedges = possibleedges(:,onesvec1);
                    expandedpossedges = reshape(expandedpossedges.',nposs,1);
                    expandedivsd = ivsd(:,onesvec1);
                    expandedivsd = reshape(expandedivsd.',nposs,1);
                    
                    if showtext >= 3
						disp(['         ',int2str(nposs),' IS+edge segments found initially '])
                    end
                    
                    if nposs > 0
                        fromcoords = ISCOORDS(expandedreftoISCOORDS,:);
                        [tocoords,edgeweightlist,~] = EDgetedgepoints(edgedata.edgestartcoords(possibleedges,:),edgedata.edgeendcoords(possibleedges,:),edgedata.edgelengthvec(possibleedges,:),nedgesubs);
                        clear possibleedges
                        
                         [hitplanes,reflpoints,edgehits,edgehitpoints,edgehitnumbers,cornerhits,cornerhitpoints,cornerhitnumbers] = ...
                             EDchkISvisible(fromcoords,tocoords,planedata.planeeqs(expandedposscombs(:,1),4),planenvecs(expandedposscombs(:,1),:),planedata.minvals(expandedposscombs(:,1),:),...
 						    planedata.maxvals(expandedposscombs(:,1),:),planedata.planecorners(expandedposscombs(:,1),:),planedata.corners,planedata.ncornersperplanevec(expandedposscombs(:,1)));
                        if ~isempty(edgehits) || ~isempty(cornerhits)
                            disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                            disp('         handled correctly yet.')
                        end

%                         eval(['reflpoints',JJ(1,1),' = reflpoints;'])
%                         reflpoints1 = reflpoints;
				
                        expandedivsd          = expandedivsd(hitplanes);
                        expandedposscombs     = expandedposscombs(hitplanes,:);
                        expandedreftoISCOORDS = expandedreftoISCOORDS(hitplanes);
                        expandedpossedges     = expandedpossedges(hitplanes);
                        edgeweightlist        = edgeweightlist(hitplanes);
                        toedgecoords          = tocoords(hitplanes,:);
				
                        nposs = length(expandedivsd);
				
                    end
                    
                    if showtext >= 3
						disp(['         ',int2str(nposs),' IS+edge segments visible '])
                    end
			
                    if obstructtestneeded && nposs > 0
                        
                        for jj = 1:2
                            if nposs > 0
                                if jj==1
                                    fromcoords = S;    
                                    startplanes = [];    
                                else
                                    startplanes = expandedposscombs(:,jj-1);
                                    eval(['fromcoords = reflpoints',int2str(jj-1),';'])
                                end
                                if jj == 2
                                    tocoords = toedgecoords;
                                    endplanes = [edgedata.planesatedge(expandedpossedges,1) edgedata.planesatedge(expandedpossedges,2)];
                                else
                                    eval(['tocoords = reflpoints',int2str(jj),';'])    
                                    endplanes = expandedposscombs(:,jj);    
                                end

                                [nonobstructedpaths,nobstructions] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,planedata.canplaneobstruct,planedata.planeseesplane,...
                                    planedata.planeeqs,planenvecs,planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.rearsideplane);
                                
                                if nobstructions > 0
                                    expandedivsd          = expandedivsd(nonobstructedpaths);
                                    expandedposscombs     = expandedposscombs(nonobstructedpaths,:);
                                    expandedreftoISCOORDS = expandedreftoISCOORDS(nonobstructedpaths);
                                    expandedpossedges     = expandedpossedges(nonobstructedpaths);
                                    edgeweightlist        = edgeweightlist(nonobstructedpaths);
                                    toedgecoords          = tocoords(nonobstructedpaths,:);
                                    nposs = length(expandedivsd);
                                    
                                    for kk = 1:1 %norder,
                                        eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(nonobstructedpaths,:);'])    
                                    end
                                    
                                end
                            end
                 
                        end
                    end        
                    
                    if showtext >= 3
					    disp(['         ',int2str(nposs),' IS+edge segments survived the obstruction test'])
                    end
                    
                    % There are repetitions of each combination since each edge segment has
                    % been treated separately. Pack them together now
			
                    test = [expandedposscombs expandedpossedges];
                    ncombs = length(expandedpossedges);
                    dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                    ivremove = find(dtest==1);
                    
                    while ~isempty(ivremove) && ~isempty(edgeweightlist)
                        edgeweightlist(ivremove+1) = double(edgeweightlist(ivremove+1)) + double(edgeweightlist(ivremove));
                        edgeweightlist(ivremove) = [];
                        expandedpossedges(ivremove) = [];
                        expandedposscombs(ivremove,:) = [];
                        expandedivsd(ivremove) = [];
                        nposs = length(expandedivsd);   
                        
                        test = [expandedposscombs expandedpossedges];
                        ncombs = length(expandedpossedges);
                        dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                        ivremove = find(dtest==1);
			
                    end
				       
                    if showtext >= 3
						disp(['         ',int2str(nposs),' IS+edge combinations after edge segments are packed together'])
                    end
                    combstoremove = setdiff(ivsd,expandedivsd);
                    addtocomb(combstoremove,:) = [];
                    addtooriginsfrom(combstoremove) = [];  
                    addtoISESVISIBILITY(combstoremove) = [];
                    nadds = length(addtooriginsfrom);
                end
                
            end
            
			% Combinations of diffraction-diffraction don't need to be checked because
			% 'edgeseesedge' should have taken care of this.
			
            %----------------------------------------------------------------------
            % Add the new combinations to the master list
		
            if nadds > 0
                POTENTIALISES = [POTENTIALISES;[addtocomb zeros(nadds,specorder-2)]];
				ORIGINSFROM = [ORIGINSFROM;addtooriginsfrom];
			                
                if difforder > 0
                    ivtemp = find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes);
                    if ~isempty(ivtemp)       
                        addtoISESVISIBILITY(ivtemp) = edgeweightlist;
                    end
%                     ivtemp = find(addtocomb(:,1)>nplanes);
                end
                ISESVISIBILITY = [ISESVISIBILITY;addtoISESVISIBILITY];                

				endindices(2) = length(ORIGINSFROM);
			
				% Compute the IS coords of the combinations with only specular reflections
				
				ISCOORDStoadd = zeros(nadds,3);
				if difforder > 0
                    ivss = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)<=nplanes));
				else
                    ivss = uint32([1:nadds].');    
				end
				soutoreflect = ISCOORDS(ORIGINSFROM(double(ivss)+endindices(1)),1:3);
%				ISCOORDStoadd(ivss,:) = ESIE2findis(soutoreflect,POTENTIALISES(double(ivss)+endindices(1),2),planedata.planeeqs);
				ISCOORDStoadd(ivss,:) = EDfindis(soutoreflect,POTENTIALISES(double(ivss)+endindices(1),2),planedata.planeeqs);
                clear soutoreflect
				ISCOORDS = [ISCOORDS;ISCOORDStoadd];
		
            end
        end
    end
end

%##################################################################
%##################################################################
%##################################################################
%
% Third order and higher
%
%---------------------------------------------------------------------------------------------


for ordernum = 3:specorder
    if showtext >= 3
        disp(['     Order number ',int2str(ordernum)])
    end
        
    startindices(ordernum) = startindices(ordernum-1) + length(addtocomb);

    % The lines below could use
    % planestoprop = unique(addtocomb(:,ordernum-1));
    
    planelist = sort(addtocomb(:,ordernum-1));
    
    if ~isempty(planelist)           
        planestoprop = planelist([1;find(diff(double(planelist))~=0)+1]);
    else
        planestoprop = [];    
    end
    clear planelist

	maxnumberofvispla = max(sum(PLANEseesPLANE(:,planestoprop)==1));
    oldaddtocomb = addtocomb;

	nadds = size(addtocomb,1);
    if nhighval < 255
        addtocomb = zeros(nadds*maxnumberofvispla,ordernum,'uint8');
    else
        addtocomb = zeros(nadds*maxnumberofvispla,ordernum,'uint16');
    end

    onesvec2 = ones(1,maxnumberofvispla);
	for ii = 1:ordernum-1
		temp = oldaddtocomb(:,ii);
		addtocomb(:,ii) = reshape(temp(:,onesvec2).',length(temp)*maxnumberofvispla,1);
	end
    nrows = size(addtocomb,1);
	clear oldaddtocomb temp
    
    if nrows > 0
       	startpos = [[1:maxnumberofvispla:nrows]  nrows+1  ].';
    else
        startpos = 1;    
    end
    
    for ii = 1:length(startpos)-1
%         addtocomb(startpos(ii),ordernum-1);
% 		possibleplanes = find( (PLANEseesPLANE(:,addtocomb(startpos(ii),ordernum-1))==1) );
% 		nposs = length(possibleplanes);
% 		addtocomb( startpos(ii):startpos(ii)+nposs-1,ordernum) = possibleplanes;
            if nvisPLANES(addtocomb(startpos(ii),ordernum-1)) > 0
    			possibleplanes = listofvisPLANES(addtocomb(startpos(ii),ordernum-1),1:nvisPLANES(addtocomb(startpos(ii),ordernum-1))).';
			    nposs = length(possibleplanes);
			    addtocomb( startpos(ii):startpos(ii)+nposs-1,ordernum) = possibleplanes;
            end
    end

    addtooriginsfrom = uint32([1:nadds].' + startindices(ordernum-1)-1);
    addtooriginsfrom = addtooriginsfrom(:,ones(1,maxnumberofvispla));
    addtooriginsfrom = reshape(addtooriginsfrom.',maxnumberofvispla*nadds,1);

    addtoISESVISIBILITY = reshape(addtoISESVISIBILITY(:,ones(1,maxnumberofvispla)).',nadds*maxnumberofvispla,1);
    
	clear startpos
	cutout = find(addtocomb(:,ordernum)==0);
	addtocomb(cutout,:) = [];
    addtooriginsfrom(cutout) = [];
    addtoISESVISIBILITY(cutout) = [];
    nadds = size(addtocomb,1);
    
	if ordernum == specorder
		clear cutout
	end

    %--------------------------------------------------------------------------
    % Try to get rid of some combinations that we know can not be propagated

    % Find those combinations of specular-specular where the IS is behind the second reflecting plane
    % because they can then never ever be propagated.

    if difforder > 0
        ivss = uint32(find(prod( double(addtocomb<=nplanes).' ).'));
    else
        ivss = uint32([1:nadds].');
    end
    
    imsoucheck = dot((ISCOORDS(addtooriginsfrom(ivss),1:3) - planedata.corners(planedata.planecorners(addtocomb(ivss,ordernum),1),:)).',planenvecs(addtocomb(ivss,ordernum),:).').';
	GEOMACC = 1e-10;
	imsounotOK = uint32(find(imsoucheck < GEOMACC));
    clear imsoucheck
	addtocomb(ivss(imsounotOK),:) = [];
	addtooriginsfrom(ivss(imsounotOK),:) = [];
    addtoISESVISIBILITY(ivss(imsounotOK)) = [];
    clear imsounotOK

    % Addition 10feb10: we should check all combinations with dss in the
    % three last orders of addtocomb: if the edge in the first 'd' belongs to the
    % plane of the last 's' *and* the two 's' planes form a 90 degree
    % corner, then these combos should be removed.
    % 
    % As a first step, we check if there are any interior-90-degree
    % planedata.corners/edges of the model, because only then can these combos occur.
    
    if ordernum == 3
       deg90edges = find( abs(edgedata.closwedangvec-3*pi/2) < GEOMACC );
       if any(deg90edges)
           deg90planepairs = edgedata.planesatedge(deg90edges,:);    
       end
    end    
    
    if difforder > 0 && any(deg90edges)
        ivdss = uint32(find(    double(addtocomb(:,ordernum-2)>nplanes).*double(addtocomb(:,ordernum-1)<=nplanes).*double(addtocomb(:,ordernum)<=nplanes)      ));
    else
        ivdss = [];
    end

    if ~isempty(ivdss)
        
        % First we check if the 'd' (the edge) belongs to the last 's'
        A = addtocomb(ivdss,ordernum-2:ordernum);

        ivsubset = find( (edgedata.planesatedge(A(:,1)-nplanes,1) == A(:,3)) | (edgedata.planesatedge(A(:,1)-nplanes,2) == A(:,3)) );
        ivdss = ivdss(ivsubset);

        A = A(ivsubset,:);
        
        % Bug found 20150919: The A matrix should be reduced in the same way
        % as ivdss! This has been fixed by the line above.

        if ~isempty(ivdss)
            % Second, we check if the two 's' planes form a 90 deg corner
            planesare90degpairs = ismember( A(:,2),deg90planepairs(:,1) ).*ismember( A(:,3),deg90planepairs(:,2) ) + ismember( A(:,2),deg90planepairs(:,2) ).*ismember( A(:,3),deg90planepairs(:,1) );
            ivsubset = planesare90degpairs;

            ivdss = ivdss(ivsubset);

            if ~isempty(ivdss)
                addtocomb(ivdss,:) = [];
                addtooriginsfrom(ivdss,:) = [];
                addtoISESVISIBILITY(ivdss) = [];
                clear ivdss ivsubset planesare90degpairs                
            end
        end
    end    
    
    nadds = size(addtocomb,1);
    
	% Combinations with all-specular plus diffraction as the last one:
    % check if the IS is behind both of the two planes that the reflecting edge connects to.
    %
    % Bug found 080711: The visibility test must be split into two, depending on the wedge angle!! 
    % If the open wedge angle > pi, then it is correct to check if the IS behind both of the planes. 
    % If the open wedge angle < pi, then it should be checked if the IS behind one of the planes!! 

    if difforder > 0
        
        ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));

        if ~isempty(ivsd)

             % First we check the sssssd combinations where the open wedge angle > pi
             % For those cases, we can exclude combinations where the IS
             % is behind both planes
                
             ivsd1 = ivsd(    edgedata.closwedangvec( double(addtocomb(ivsd,ordernum))-nplanes ) < pi );
%            ivsd2 = ivsd( find(   closwedangvec( double(addtocomb(ivsd,ordernum))-nplanes ) > pi ));

             if ~isempty(ivsd1)
                imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)    ,:).').';
                imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)    ,:).').';
                GEOMACC = 1e-10;
                imsounotOK = find(imsoucheck1 < GEOMACC & imsoucheck2 < GEOMACC);
                clear imsoucheck1 imsoucheck2
                addtocomb(ivsd1(imsounotOK),:) = [];
                addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                clear imsounotOK
             end
             
             % Second we check the sssssd combinations where the open wedge angle < pi
             % For those cases, we can exclude combinations where the IS
             % is behind one of the two planes
                
             ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));
             if ~isempty(ivsd)
                 ivsd1 = ivsd(    edgedata.closwedangvec( double(addtocomb(ivsd,ordernum))-nplanes ) > pi );

                 if ~isempty(ivsd1)
                    imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)    ,:).').';
                    imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - planedata.corners(planedata.planecorners(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)     ,1),:)).',planenvecs(    edgedata.planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)    ,:).').';
                    GEOMACC = 1e-10;
                    imsounotOK = find(imsoucheck1 < GEOMACC | imsoucheck2 < GEOMACC);
                    clear imsoucheck1 imsoucheck2
                    addtocomb(ivsd1(imsounotOK),:) = [];
                    addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                    addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                    clear imsounotOK
                 end
             end
             
            % For combinations including dsd, if the two diff. edges are
            % aligned and the intermediate specular reflections are
            % caused by planes that are perpendicular to the edges,
            % remove these combinations
            
            if difforder >= 2
                for kk = 1:ordernum-2
                    ncombs = size(addtocomb,1); 
                    matchvec = ones(ncombs,1);
%                     nprespecs = kk-1;
                    nmidspecs = ordernum-1-kk;
                    matchvec = matchvec.*(addtocomb(:,ordernum)>nplanes).*(addtocomb(:,kk)>nplanes);
                    if nmidspecs == 1
                        matchvec = matchvec.*(addtocomb(:,2)<=nplanes);
                    elseif nmidspecs >= 2
                        matchvec = matchvec.*prod( double(addtocomb(:,kk+1:ordernum-1)<=nplanes).' ).';
                    end
                    ivdsd = uint32(find( matchvec ));
                    if ~isempty(ivdsd)
                        edge1 = double(addtocomb(ivdsd,kk)) - nplanes;
                        edge2 = double(addtocomb(ivdsd,ordernum)) - nplanes;
                        midspecs = double(addtocomb(ivdsd,kk+1:ordernum-1));
                        
                        lookupindmat1 = (edge1-1)*nedges + edge2;
                        matrixcolnumbers = (edge1-1)*nplanes;
                        lookupindmat2 = matrixcolnumbers(:,ones(1,nmidspecs)) + midspecs;                    
                        
                        if nmidspecs == 1
                            ivinvalid = find(edgetoedgedata.edgealignedwithedge(lookupindmat1) & edgetoedgedata.edgeperptoplane(lookupindmat2));
                        else
                            ivinvalid = find(edgetoedgedata.edgealignedwithedge(lookupindmat1) & prod(edgetoedgedata.edgeperptoplane(lookupindmat2).').');
                        end
            			addtocomb(ivdsd(ivinvalid),:) = [];
				        addtooriginsfrom(ivdsd(ivinvalid),:) = [];
                        addtoISESVISIBILITY(ivdsd(ivinvalid)) = [];
                    end
                 end
            end            
            
            % We need to find the valid combinations again after we have removed a
            % number of combinations
	
            ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));
	
            %----------------------------------------------------------------------
            % Combinations of specular-diffraction that are not visible at
            % all should be removed and then they are not propagated either.
		
            if nhighval < 255
                possibleedges = uint8(double(addtocomb(ivsd,ordernum)) - nplanes);
                possiblecombs = uint8(double(addtocomb(ivsd,1:ordernum-1)));
            else
                possibleedges = uint16(double(addtocomb(ivsd,ordernum)) - nplanes);
                possiblecombs = uint16(double(addtocomb(ivsd,1:ordernum-1)));
            end
            reftoISCOORDS = addtooriginsfrom(ivsd);
		
            % Expand to take the edge subdivisions into account
            
            nposs = length(ivsd);
            nposs = nposs*nedgesubs;        % We treat all the little edge subdivisions as separate edges

            expandedposscombs = reshape(repmat(possiblecombs.',[nedgesubs,1]),ordernum-1,nposs).';
            clear possiblecombs                        
            expandedreftoISCOORDS = reftoISCOORDS(:,onesvec1);
            expandedreftoISCOORDS = reshape(expandedreftoISCOORDS.',nposs,1);
            expandedpossedges = possibleedges(:,onesvec1);
            expandedpossedges = reshape(expandedpossedges.',nposs,1);
            expandedivsd = ivsd(:,onesvec1);
            expandedivsd = reshape(expandedivsd.',nposs,1);
            
            if showtext >= 3
				disp(['         ',int2str(nposs),' IS+edge segments found initially '])
            end
	
            % Go through, iteratively, and check if the path from S to edge passes
            % through all reflection planes along the way
            
            for jj = ordernum-1:-1:1
	
                if nposs > 0
		
                    if jj == ordernum-1
                        fromcoords = ISCOORDS(expandedreftoISCOORDS,:);
                        [tocoords,edgeweightlist,edgenumberlist] = EDgetedgepoints(edgedata.edgestartcoords(possibleedges,:),edgedata.edgeendcoords(possibleedges,:),edgedata.edgelengthvec(possibleedges,:),nedgesubs);
                        clear possibleedges
                    else
                        eval(['tocoords = reflpoints',int2str(jj+1),';'])    
                        ivlist = ORIGINSFROM(expandedreftoISCOORDS);
                        for kk = jj:ordernum-3
                            ivlist = ORIGINSFROM(ivlist);                            
                        end
                        fromcoords = ISCOORDS(ivlist,:);
                                           
                    end

                    [hitplanes,reflpoints,edgehits,edgehitpoints,edgehitnumbers,...
                        cornerhits,cornerhitpoints,cornerhitnumbers] = EDchkISvisible(fromcoords,tocoords,planedata.planeeqs(expandedposscombs(:,jj),4),planenvecs(expandedposscombs(:,jj),:),planedata.minvals(expandedposscombs(:,jj),:),...
					    planedata.maxvals(expandedposscombs(:,jj),:),planedata.planecorners(expandedposscombs(:,jj),:),planedata.corners,planedata.ncornersperplanevec(expandedposscombs(:,jj)));
                    if ~isempty(edgehits) || ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                        disp('         handled correctly yet.')
                    end
                    eval(['reflpoints',int2str(jj),' = reflpoints;'])
			
                    expandedivsd          = expandedivsd(hitplanes);
                    expandedposscombs     = expandedposscombs(hitplanes,:);
                    expandedreftoISCOORDS = expandedreftoISCOORDS(hitplanes);
                    expandedpossedges     = expandedpossedges(hitplanes);
                    edgeweightlist        = edgeweightlist(hitplanes);
                    toedgecoords          = tocoords(hitplanes,:);
                    if jj < ordernum-1
                        for kk = jj+1:ordernum-1
                           eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(hitplanes,:);'])
                        end
                    end
		
                    nposs = length(expandedivsd);
                    
                end
                if showtext >= 3
				    disp(['         ',int2str(nposs),' IS+edge segments survived the visibility test in refl plane ',int2str(jj)])
                end
                
            end
	
            if obstructtestneeded && nposs > 0
                for jj = 1:ordernum
                    if nposs > 0
                        
                        if jj==1
                            fromcoords = S;    
                            startplanes = [];    
                        else
                            startplanes = expandedposscombs(:,jj-1);
                            eval(['fromcoords = reflpoints',int2str(jj-1),';'])
                        end
                        if jj == ordernum
                            tocoords = toedgecoords;
                            endplanes = [edgedata.planesatedge(expandedpossedges,1) edgedata.planesatedge(expandedpossedges,2)];
%                             endplanes = [];    
                        else
                            eval(['tocoords = reflpoints',int2str(jj),';'])    
                            endplanes = expandedposscombs(:,jj);    
                        end
             
                        [nonobstructedpaths,nobstructions] = EDcheckobstrpaths(fromcoords,tocoords,startplanes,endplanes,planedata.canplaneobstruct,planedata.planeseesplane,...
                            planedata.planeeqs,planenvecs,planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.rearsideplane);
		
                        if nobstructions > 0
                            expandedivsd          = expandedivsd(nonobstructedpaths);
                            expandedposscombs     = expandedposscombs(nonobstructedpaths,:);
                            expandedreftoISCOORDS = expandedreftoISCOORDS(nonobstructedpaths);
                            expandedpossedges     = expandedpossedges(nonobstructedpaths);
                            edgeweightlist        = edgeweightlist(nonobstructedpaths);
                            toedgecoords          = tocoords(nonobstructedpaths,:);
                            nposs = length(expandedivsd);
                            for kk = 1:ordernum-1
                                eval(['reflpoints',int2str(kk),' = reflpoints',int2str(kk),'(nonobstructedpaths,:);'])    
                            end
                            
                        end
                    end
                end
            end        
        
            if showtext >= 3
		        disp(['         ',int2str(nposs),' IS+edge segments survived the obstruction test'])
            end

            % There are repetitions of each combination since each edge segment has
            % been treated separately. Pack them together now
	
            test = [expandedposscombs expandedpossedges];
            ncombs = length(expandedpossedges);
            dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
            ivremove = find(dtest==1);
            
            while ~isempty(ivremove)
           
                edgeweightlist(ivremove+1) = double(edgeweightlist(ivremove+1)) + double(edgeweightlist(ivremove));
                edgeweightlist(ivremove) = [];
                expandedpossedges(ivremove) = [];
                expandedposscombs(ivremove,:) = [];
                expandedivsd(ivremove) = [];
                nposs = length(expandedivsd);   
                
                test = [expandedposscombs expandedpossedges];
                ncombs = length(expandedpossedges);
                dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                ivremove = find(dtest==1);
	
            end
            	       
            if showtext >= 3
				disp(['         ',int2str(nposs),' IS+edge segments after edge segments are packed together'])
            end
	
            combstoremove = setdiff(ivsd,expandedivsd);
            addtocomb(combstoremove,:) = [];
            addtooriginsfrom(combstoremove) = [];  
            addtoISESVISIBILITY(combstoremove) = [];
            nadds = length(addtooriginsfrom);
        end        
    end   

    POTENTIALISES = [POTENTIALISES;[addtocomb zeros(nadds,specorder-ordernum)]];
    ORIGINSFROM = [ORIGINSFROM;addtooriginsfrom];

    ivtemp = find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes));
    if difforder > 0
        if ~isempty(ivtemp)
            addtoISESVISIBILITY(ivtemp) = edgeweightlist;
        end
    end
    ISESVISIBILITY = [ISESVISIBILITY;addtoISESVISIBILITY];
    
    endindices(ordernum) = length(ORIGINSFROM);

    % Compute the IS coords

    ISCOORDStoadd = zeros(nadds,3);
    if difforder > 0
        ivss = uint32(find(prod( double(addtocomb<=nplanes).' ).'));
    else
        ivss = uint32([1:nadds].');    
    end
    soutoreflect = ISCOORDS(ORIGINSFROM(double(ivss)+endindices(ordernum-1)),1:3);
%    ISCOORDStoadd(ivss,:) = ESIE2findis(soutoreflect,POTENTIALISES(double(ivss)+endindices(ordernum-1),ordernum),planedata.planeeqs);
    ISCOORDStoadd(ivss,:) = EDfindis(soutoreflect,POTENTIALISES(double(ivss)+endindices(ordernum-1),ordernum),planedata.planeeqs);
    clear soutoreflect
    ISCOORDS = [ISCOORDS;ISCOORDStoadd];
    clear ISCOORDStoadd
    
end

% ntot = size(POTENTIALISES,1);

%-------------------------------------------------------
% Add the direct sound to the start of the list

POTENTIALISES = [zeros(1,specorder);POTENTIALISES];
ORIGINSFROM = [0;ORIGINSFROM];
startindices = startindices+1;
endindices = endindices+1;
ORIGINSFROM = uint32(double(ORIGINSFROM)+1);

ISCOORDS = [S;ISCOORDS];
ISESVISIBILITY = [0;ISESVISIBILITY];

%-------------------------------------------------------
% Trim the list if there are zeros-only in the last column(s)

n2 = size(POTENTIALISES,2);

if n2 > 1
    columnstokeep = sum(POTENTIALISES)~=0;
    POTENTIALISES = POTENTIALISES(:,columnstokeep);    
    n2new = size(POTENTIALISES,2);    
    if n2new < n2
%         specorderold = specorder;
        specorder = n2new;
    end
end

%-------------------------------------------------------
% Create index lists

% First, the total reflection order

[n1,n2] = size(POTENTIALISES);

lengthNspecmatrix = [];
lengthNdiffmatrix = [];
singlediffcol = [];
startindicessinglediff = [];
endindicessinglediff = [];

if n1 > 0 && n2 > 0
if n2 > 1
    REFLORDER = uint8(sum(POTENTIALISES.'>0).');    
else
    REFLORDER = uint8(POTENTIALISES>0);
end

if difforder > 0
    if specorder > 1
        ivss = uint32(find(prod( double(POTENTIALISES<=nplanes).' ).'));
    else
        ivss = uint32(find( POTENTIALISES<=nplanes ));    
    end

    ivother = uint32([1:length(ORIGINSFROM)].');
    ivother(ivss) = [];
    ivss(1) = [];
    IVNDIFFMATRIX = zeros(2,specorder,'uint32');
    lengthNdiffmatrix = zeros(1,specorder);
    nrowsdiff = 2;
    IVNSPECMATRIX = zeros(2,specorder,'uint32');
    lengthNspecmatrix = zeros(1,specorder);
    nrowsspec = 2;
    
    for ii = 1:specorder
        if specorder > 1
            ivreftoNdiff = uint32(find(sum( double(POTENTIALISES(ivother,:)> nplanes).' ).'==ii));
        else
            ivreftoNdiff = uint32(find(( double(POTENTIALISES(ivother,:)> nplanes).' ).'==ii));
        end
        ivreftoNspec = uint32(find(REFLORDER(ivss)==ii));

        ivNdiff = ivother(ivreftoNdiff);
        if ~isempty(ivNdiff)
            lengthNdiffmatrix(ii) = length(ivNdiff);
            if lengthNdiffmatrix(ii) > nrowsdiff
               IVNDIFFMATRIX = [IVNDIFFMATRIX;zeros(lengthNdiffmatrix(ii)-nrowsdiff,specorder,'uint32')];
               nrowsdiff = lengthNdiffmatrix(ii);
            end
            IVNDIFFMATRIX(1: lengthNdiffmatrix(ii),ii) =  ivNdiff;     
        end
        ivother(ivreftoNdiff) = [];
        
        ivNspec = ivss(ivreftoNspec);
        if ~isempty(ivNspec)
            lengthNspecmatrix(ii) = length(ivNspec);
            if lengthNspecmatrix(ii) > nrowsspec
                IVNSPECMATRIX = [IVNSPECMATRIX;zeros(lengthNspecmatrix(ii)-nrowsspec,specorder,'uint32')];
                nrowsspec = lengthNspecmatrix(ii);
            end
            IVNSPECMATRIX(1: lengthNspecmatrix(ii),ii) =  ivNspec;      
        end        
    end
     
    % Determine which column the single diffraction is in
    
    test = POTENTIALISES(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1),:)>nplanes;
       
    colweight = (1:specorder);
    nsingles = length(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1));
    if specorder > 1
        singlediffcol = uint8(sum( (test.*colweight(ones(nsingles,1),:)).').');
    else
        singlediffcol = uint8(ones(nsingles,1));        
    end
        
    startindicessinglediff = zeros(specorder,1);
    endindicessinglediff = zeros(specorder,1);
    startindicessinglediff(1) = 1;
    for ii = 1:specorder-1
        iv = find(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1)<=endindices(ii));
        if ~isempty(iv)
            endindicessinglediff(ii) = iv(length(iv));
        else
            endindicessinglediff(ii) = 0;
        end
        startindicessinglediff(ii+1) = endindicessinglediff(ii)+1;
    end
    endindicessinglediff(specorder) = lengthNdiffmatrix(1);
else
    ivss = uint32([2:endindices(length(endindices))].');
    IVNSPECMATRIX = zeros(2,specorder,'uint32');
    lengthNspecmatrix = zeros(1,specorder);
    nrowsspec = 2;

    for ii = 1:specorder
        ivreftoNspec = uint32(find(REFLORDER(ivss)==ii));

        ivNspec = ivss(ivreftoNspec);
        lengthNspecmatrix(ii) = length(ivNspec);
        if lengthNspecmatrix(ii) > nrowsspec
           IVNSPECMATRIX = [IVNSPECMATRIX;zeros(lengthNspecmatrix(ii)-nrowsspec,specorder,'uint32')];
           nrowsspec = lengthNspecmatrix(ii);
        end
        IVNSPECMATRIX(1: lengthNspecmatrix(ii),ii) =  ivNspec;      
    
    end
    
    IVNDIFFMATRIX = [];
    lengthNdiffmatrix = [];

    singlediffcol = [];
    startindicessinglediff = [];
    endindicessinglediff = [];
end

end

%-------------------------------------------------------
% Create a pointer list so that it is possible to find a combination
% that ends with edge-plane2, given an edge-plane1-plane2 combination
%
% NB!! The list points to the original POTENTIALISES (and related lists)

if difforder > 0 && specorder > 1

    if showtext >= 3
        disp('       Building a pointer list for image receiver combinations')    
    end

    % First select only those that have a single diff combo and that
    % start with the diffraction
    
    ndecimaldivider = (nplanes+2);
        
% % %     PointertoIRcombs = sparse(zeros(ndecimaldivider^(specorder-1),1));
    PointertoIRcombs = sparse(zeros(1000,1));
    IRoriginsfrom = zeros(size(ORIGINSFROM),'uint32');
    
	for ii = 1:specorder-1
        if showtext >= 3
            disp(['          order ',int2str(ii)])    
        end
                
        if ii == 1
            ivrange = uint32([startindicessinglediff(2):endindicessinglediff(2)]);
%            masterivlistselect = ivsinglediff(ivrange);
            masterivlistselect = IVNDIFFMATRIX(ivrange,1);
            masterivlistselect = masterivlistselect(singlediffcol(ivrange)==1);

            % IRoriginsfrom should be given the value zero for these
            % first-order specular combos (after one diff) so we don't
            % bother assigning any value
            
%             A = POTENTIALISES(masterivlistselect,2); 
%             [B,IA,JA] = unique(A);
            [B,IA,~] = unique(POTENTIALISES(masterivlistselect,2));
            listlength = max(B);
            if listlength > length(PointertoIRcombs)
               disp(['Expanding to length ',int2str(listlength)])
               PointertoIRcombs = [PointertoIRcombs;sparse(zeros(listlength-length(PointertoIRcombs),1))    ]; 
            end
            PointertoIRcombs(B) = masterivlistselect(IA);
            
            % Now we check if there any active wall numbers that don't
            % occur in an edge-spec combination. For those combinations
            % we can point to a specular combos instead, that ends with the
            % same plane number.
            
            wallist = POTENTIALISES(IVNSPECMATRIX(1:lengthNspecmatrix(1),1),1);
            ivnotincludedyet = find(~ismember(wallist,B));
            if ~isempty(ivnotincludedyet)
                PointertoIRcombs(wallist(ivnotincludedyet)) = IVNSPECMATRIX(ivnotincludedyet,1);
            end
            
        elseif ii >= 2
            ivrange = uint32(startindicessinglediff(ii+1):endindicessinglediff(ii+1));
            masterivlistselect = IVNDIFFMATRIX(ivrange,1);
            masterivlistselect = masterivlistselect(singlediffcol(ivrange)==1);
 
            ivlist = 0;
            for jj = 3:ii+1
                ivlist = ivlist + double(POTENTIALISES(masterivlistselect,jj))*ndecimaldivider^(ii+1-jj);
            end
            
            IRoriginsfrom(masterivlistselect) = uint32(full(PointertoIRcombs(ivlist)));

            ivlist = ivlist + double(POTENTIALISES(masterivlistselect,2))*ndecimaldivider^(ii-1);
            A = uint32(ivlist);
            [B,IA,~] = unique(A);
            PointertoIRcombs(B) = masterivlistselect(IA);

            % Now we check if there any active wall numbers that don't
            % occur as spec1 in an edge-spec1-spec2-spec3 combination.

            wallist = 0;
            for jj = 1:ii
                wallist = wallist + double(POTENTIALISES(IVNSPECMATRIX(1:lengthNspecmatrix(ii),ii),jj))*ndecimaldivider^(ii-jj);
            end

            ivnotincludedyet = find(~ismember(wallist,B));

            if ~isempty(ivnotincludedyet)
                PointertoIRcombs(wallist(ivnotincludedyet)) = IVNSPECMATRIX(ivnotincludedyet,ii);
            end
            
        end
	
	end
        
else
    ndecimaldivider = 0;
    PointertoIRcombs = [];
    IRoriginsfrom = [];
end

%-------------------------------------------------------
% Format the output data

maxval = max(IRoriginsfrom);
if maxval < 256
    IRoriginsfrom = uint8(IRoriginsfrom);
elseif maxval < 65536
    IRoriginsfrom = uint16(IRoriginsfrom);    
end
maxval = max(ORIGINSFROM);
if maxval < 256
    ORIGINSFROM = uint8(ORIGINSFROM);
elseif maxval < 65536
    ORIGINSFROM = uint16(ORIGINSFROM);    
end

    