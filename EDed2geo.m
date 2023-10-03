function [edgetoedgedata,outpar,existingfilename] = EDed2geo(edgedata,planedata,Sdata,...
    Rdata,specorder,difforder,EDversionnumber,nedgesubs,ndiff2batches,inpar)
% EDed2geo - Calculates 2nd- and higher-order edge-related geom. parameters.
% EDed2geo calculates some plane- and edge-related geometrical parameters,
% based on corners and planes from a EDreadcad, but only data that is 
% independent of the source and receiver.
% 
% A version 1 of this function was used up to v 0.221 of EDtoolbox
% and version 2 after that.
%
% Input parameters:
%	edgedata
%   planedata
%   Sdata                   Struct: needed if some edges could be redundant
%                           to compute data for.
%   Rdata
%   specorder               The highest order of any reflection kind (specular and/or diffraction).
%	difforder 	            The highest order of diffraction. If it is 0 or 1 then the parameter
%							edgeseespartialedge is not calculated. Default: 1
%   EDversionnumber
%   nedgesubs (optional)            Default: 2
%   ndiff2batches (optional)        Default: 1
%   inpar                   This parameter is different for version 1 and
%                           version 2 of this function. 
%       v1 of inpar should be showtext (optional)
%                           0 -> no text displayed. Value > 0 -> text displayed.
%       v2 of inpar should be filehandlingparameters (obligatory)
%                           filehandlingparameters is a struct which
%                           contains the field showtext.
% 
% Output parameters:
%   edgetoedgedata          struct with fields:
%       edgeseespartialedge     A sparse matrix, [nedges,nedges], with the
%                           values 0, pos. integer or neg. integer:
%                           0           Two edges can not see each other (or one
%                                       of the edges is inactive)
%                           pos.int.    A value from 1 to 2^nedgesubs-1,
%                                       indicating how much of the edges
%                                       see each other. A pos. value
%                                       indicates that the edge-to-edge
%                                       path runs along a plane.
%                                       NB! A complete visibility test is
%                                       not made; subsegment 1 on edge 1 is only
%                                       checked against subsegment 1 on
%                                       edge 2, not against all other
%                                       subsegments of edge 2.
%                           neg.int.    Same as for the pos.int. except
%                                       that the edge-to-edge path does not
%                                       run along a plane.
%       edgealignedwithedge     A sparse matrix, [nedges,nedges], with the value 1
%                           or 0 stating whether an (active) edge is aligned
%                           with another (active) edge or not.
%       edgeperptoplane         A sparse matrix, [nplanes,nedges], with the value 1
%                           or 0 stating whether an (active) edge is perpendicular
%                           to an active plane or not.
%       edgeplaneperptoplane1   A sparse matrix, [nplanes,nedges], with the value 1
%                           or 0 stating whether an (active) edge has one
%                           of its two defining planes perpendicular to
%                           another plane with which it shares an edge.
%       edgeplaneperptoplane2   A sparse matrix, [nplanes,nedges], with the value 1
%                           or 0 stating whether an (active) edge has one
%                           of its two defining planes perpendicular to
%                           another plane which has a flat edge (i.e., a 180 degree edge).
%       reftoshortlistE         
%       re1sho, re2sho
%       thetae1sho, thetae2sho
%       ze1sho, ze2sho
%       examplecombE	
%   outpar                  This parameter is different for version 1 and
%                           version 2 of this function. 
%       v1 outpar = EDinputdatahash       
%                           This is a string of characters which is
%                           uniquely representing the input data.
%                           An existing result file with the same value of
%                           this EDinputdatahash will be reused.
%       v2 outpar = elapsedtimeed2geo
%                           This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimeedgeo has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        For v2 of this function, if an existing file was
%                           found that was reused, then the reused file
%                           name is given here. If no existing file could be 
%                           reused then this variable is empty. For v1 of
%                           this function, this variable is also empty.
%
% Uses the functions EDcoordtrans2, EDinfrontofplane, EDcompress7,EDrecycleresultfiles
%                    EDcheckobstr_edgetoedge EDgetedgepoints in EDtoolbox
% Uses the function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Sep. 2023
%
% [edgetoedgedata,outpar,existingfilename] = EDed2geo(edgedata,planedata,Sdata,...
%    Rdata,specorder,difforder,EDversionnumber,nedgesubs,ndiff2batches,inpar)

% 27.5.2011  Stable version
% 13.3.2013  Fixed a problem with thin folded planes. Two edges with the
%            same start and end corners were marked as seeing each other.
% 28.10.2013 Empty indentingedgepairs caused problems
% 12.11.2014 When checking for aligned edges; the precision was rounded off
%            to 6 decimals.
% 1 Dec. 2014   Added the new output parameter edgerelatedcoordsysmatrices
% 3 Dec. 2014 Added the new output struct edgetoedgedata
% 31 March 2015 Added '' in the file saving, to handle file and directory
%               names with blanks
% 6 April 2017 Fixed a bug: planenvecs was not extracted from planedata
% 4 October 2017 Fixed around line 570, to work also for cases when
%                edgealignedwith edge is empty.
% % 21 October 2017 bug: the edgetoedgestrength was not used, so
            % it seemed like all edge-to-edge paths ran along planes.
% 27 Nov. 2017 Copied from ESIE2toolbox. Removed all file loading and
% saving
% 28 Nov. 2017 Introduced the non-global showtext parameter. Cleaned the
% code a bit. Removed the parameter planeseesplane to the EDcheckobstr_edgetoedge
% 29 Nov. 2017 Changed call from ESIE2getedgepoints to ED
% 15 Jan. 2018 Little modification; edge-to-edge obstruction test was
% active before, for a single plate.
% 8 Feb 2018 Introduced the EDinputdatahash
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.

t00 = clock;

if nargin < 10  % Must be the old version
	functionversion = 1;
	showtext = 0;
    if nargin < 9
        ndiff2batches = 1;
        if nargin < 8
            nedgesubs = 2;
        end
    end
else % nargin = 10 -> could be the old or new version
	if isstruct(inpar)
		functionversion = 2;
		filehandlingparameters = inpar;
        showtext = filehandlingparameters.showtext;
	else
		functionversion = 1;
		showtext = inpar;
	end
    if isempty(ndiff2batches)
        ndiff2batches = 1;
    end
    if isempty(nedgesubs)
        nedgesubs = 2;
    end
end

EDinputdatastruct = struct('planedata',planedata,'edgedata',edgedata,...
    'Sdata',Sdata,'Rdata',Rdata,'specorder',specorder,'difforder',difforder,...
    'nedgesubs',nedgesubs,'ndiff2batches',ndiff2batches,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

% geomacc = 1e-10;

if functionversion == 1
	outpar = EDinputdatahash;
	existingfilename = '';
elseif functionversion == 2
    %---------------------------------------------------------------
    % Sort out the file business: can an existing file be used?
    % Then copy the existing file to a new copy. Should the data be saved in a file? 
    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
        existingfilename = '';
    else
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_ed2data',EDinputdatahash);
    end
    
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_ed2data.mat'];
    
    if foundmatch == 1
        eval(['load(''',existingfilename,''')'])
        if ~strcmp(existingfilename,desiredname)
            copyfile(existingfilename,desiredname);
        end
        elapsedtimeed2geo_new = etime(clock,t00);
        elapsedtimeed2geo = [elapsedtimeed2geo_new elapsedtimeed2geo];
        outpar = elapsedtimeed2geo;
        return
    end
end

edgestartcoordsnudge =   edgedata.edgestartcoordsnudge;
edgeendcoordsnudge =    edgedata.edgeendcoordsnudge;
edgenormvecs =    edgedata.edgenormvecs;    

nplanes = size(planedata.planecorners,1);
nedges = size(edgedata.edgecorners,1);

zerosvec1 = zeros(nplanes,1);

planenvecs = planedata.planeeqs(:,1:3);

%---------------------------------------------------------------
% We check right away which edges that can not be seen by either the source
% or the receiver. We don't need to check edgeseesedge for combinations
% that involve any of these.
%
% Note that this is true only for difforder = 2! If difforder >= 3, then we
% can not exclude any edges based on this! (It could be possible to exclude
% edges, if we managed to "unwind" the
% edge-sees-an-edge-which-sees-an-edge-which-can-be-seen-by-the-source.

% PS: corrected 19.1.2012: for the inteq method we can not switch off any
% edges!
% % % % if difforder == 2,
% % % %     edgesthatareseen = any( [vispartedgesfroms vispartedgesfromr].' );
% % % %     edgesnottoworryabout = [1:nedges];
% % % %     edgesnottoworryabout(edgesthatareseen) = [];
% % % %     clear vispartedgesfroms vispartedgesfromr visplanesfroms visplanesfromr
% % % % else
    edgesnottoworryabout = [];    
% % % % end

%---------------------------------------------------------------
%
% edgeperptoplane
%
%---------------------------------------------------------------
% Make matrices over:
%   1. Which edges are perpendicular to planes ("edgeperptoplane")
%   2. Which edges can give combinations with edge1-plane-edge1
%      where plane is perpendicular to one of the two planes that
%      defines the edge ("edgeplaneperptoplane1")
%   3. Which edges can give combinations with edge1-plane-edge1
%      where plane is a part of a (unnecessarily) divided plane
%      ("edgeplaneperptoplane2")
%
% NB! Only active edges and active planes can get
% such combinations.
%
% These are needed for diff-spec-diff combos, so only when
% difforder >= 2 (and specorder >= 3).

edgeperptoplane = [];
edgeplaneperptoplane1 = [];
edgeplaneperptoplane2 = [];
if specorder >= 3
    edgeperptoplane = sparse(zeros(nplanes,nedges));
    edgeplaneperptoplane1 = sparse(zeros(nplanes,nedges));
    edgeplaneperptoplane2 = sparse(zeros(nplanes,nedges));    
    if showtext >= 4
		disp('         checking which edges are perpendicular to planes ...')
	end
    
    iv = find(edgedata.edgeseesplane == 1);
    if ~isempty(iv)
        edgelist = floor((iv-1)/nplanes)+1;
        planelist = iv - (edgelist-1)*nplanes;

        % Check first which edges have vectors that are parallel to plane
        % normal vectors. These are "edgeperptoplane" cases. Clear them
        % afterwards from the edgelist and planelist.
        
        vec1 = edgenormvecs(edgelist,:);
        vec2 = planenvecs(planelist,:);
        ivperp = find( abs(  sum( (vec1.*vec2).' )  )==1);        
        if ~isempty(ivperp)
            edgeperptoplane(iv(ivperp)) = ones(size(ivperp));    
            clear iv
            edgelist(ivperp) = [];
            planelist(ivperp) = [];
            vec2(ivperp,:) = [];
        else
            clear iv    
        end

        % Now check which edges have plane1 or plane2 perpendicular to
        % other planes.

        if ~isempty(planelist)
            vec1 = planenvecs(edgedata.planesatedge(edgelist,1),:);
            ivplaneperp1 = find( abs(  sum( (vec1.*vec2).' )  )==0);
            vec1 = planenvecs(edgedata.planesatedge(edgelist,2),:);
            ivplaneperp2 = find( abs(  sum( (vec1.*vec2).' )  )==0);
            ivplaneperp = [ivplaneperp1 ivplaneperp2];
            if ~isempty(ivplaneperp)
                ivplaneperp = unique(ivplaneperp);
                edgelist = edgelist(ivplaneperp);
                planelist = planelist(ivplaneperp);
                edgeplane1 = edgedata.planesatedge(edgelist,1);
                edgeplane2 = edgedata.planesatedge(edgelist,2);

                % Check which of these other perpendicular planes that connnect 
                % to one of the edge planes via a 90 deg. corner.                

                ivtoclear = [];
                [ivfoundcombos,loc] = ismember([planelist edgeplane1],edgedata.planesatedge,'rows');
                ivperp = find(ivfoundcombos);
                if ~isempty(ivperp)
                    edgenumb = loc(ivperp);
                    ivfinalcomb = ivperp(closwedangvec(edgenumb)==3*pi/2);
                    if ~isempty(ivfinalcomb)
                        ivreftomatrix = planelist(ivfinalcomb) + (edgelist(ivfinalcomb)-1)*nplanes;
                        edgeplaneperptoplane1(ivreftomatrix) = ones(size(ivreftomatrix));
                        ivtoclear = ivfinalcomb;
                    end
                end
                [ivfoundcombos,loc] = ismember([edgeplane1 planelist],edgedata.planesatedge,'rows');
                ivperp = find(ivfoundcombos);
                if ~isempty(ivperp)
                    edgenumb = loc(ivperp);
                    ivfinalcomb = ivperp(closwedangvec(edgenumb)==3*pi/2);
                    if ~isempty(ivfinalcomb)
                        ivreftomatrix = planelist(ivfinalcomb) + (edgelist(ivfinalcomb)-1)*nplanes;
                        edgeplaneperptoplane1(ivreftomatrix) = ones(size(ivreftomatrix));

                        %%%%%%%%% Check this 050112
                        ivtoclear = ivtoclear.';
                        ivfinalcomb = ivfinalcomb.';
                        ivtoclear = [ivtoclear ivfinalcomb];
                    end
                end
                
                [ivfoundcombos,loc] = ismember([planelist edgeplane2],edgedata.planesatedge,'rows');
                ivperp = find(ivfoundcombos);
                if ~isempty(ivperp)
                    edgenumb = loc(ivperp);
                    ivfinalcomb = ivperp(closwedangvec(edgenumb)==3*pi/2);
                    if ~isempty(ivfinalcomb)
                        ivreftomatrix = planelist(ivfinalcomb) + (edgelist(ivfinalcomb)-1)*nplanes;
                        edgeplaneperptoplane1(ivreftomatrix) = ones(size(ivreftomatrix));

                        %%%%%%%%% Check this 050112
                        if size(ivtoclear,2) == size(ivfinalcomb,2)
                            ivtoclear = [ivtoclear;ivfinalcomb];
                        elseif size(ivtoclear,1) == size(ivfinalcomb,1)                    
                            ivfinalcomb = ivfinalcomb.';
                            ivtoclear = [ivtoclear ivfinalcomb];
                        elseif isempty(ivtoclear)
                            error('ERROR: Unexpected error involving sizes of ivtoclear and ivfinalcomb')    
                        end
                    end
                end
                
                
                [ivfoundcombos,loc] = ismember([edgeplane2 planelist],edgedata.planesatedge,'rows');
                ivperp = find(ivfoundcombos);
                if ~isempty(ivperp)
                    edgenumb = loc(ivperp);
                    ivfinalcomb = ivperp(closwedangvec(edgenumb)==3*pi/2);
                    if ~isempty(ivfinalcomb)
                        ivreftomatrix = planelist(ivfinalcomb) + (edgelist(ivfinalcomb)-1)*nplanes;
                        edgeplaneperptoplane1(ivreftomatrix) = ones(size(ivreftomatrix));

                        %%%%%%%%% Check this 050112
                        if size(ivtoclear,2) == size(ivfinalcomb,2)
                            ivtoclear = [ivtoclear;ivfinalcomb];
                        elseif size(ivtoclear,1) == size(ivfinalcomb,1)                  
                            ivfinalcomb = ivfinalcomb.';
                            ivtoclear = [ivtoclear ivfinalcomb];
                        elseif isempty(ivtoclear)
                            error('ERROR: Unexpected error involving sizes of ivtoclear and ivfinalcomb')    
                        end
                    end
                end
                
                % Finally check which edges have plane1 or plane2 perpendicular to
                % two other planes that connnect via a 180 deg. corner.

                if ~isempty(ivtoclear)
                    ivtoclear = unique(ivtoclear);
                    edgelist(ivtoclear) = [];
                    planelist(ivtoclear) = [];
                    if ~isempty(edgelist)
                        
                        ivflatedges = find(closwedangvec==pi);
                        if ~isempty(ivflatedges)
                            planecombos = planesatedge(ivflatedges,:);
                            for ii = 1:length(ivflatedges)
                                iv1 = find(planelist == planecombos(ii,1));                                
                                if ~isempty(iv1)
                                    ivreftomatrix = planecombos(ii,1) + (edgelist(iv1)-1)*nplanes;
                                    edgeplaneperptoplane2(ivreftomatrix) = ones(size(ivreftomatrix));                                  
                                end
                                iv2 = find(planelist == planecombos(ii,2));
                                if ~isempty(iv2)
                                    ivreftomatrix = planecombos(ii,2) + (edgelist(iv2)-1)*nplanes;
                                    edgeplaneperptoplane2(ivreftomatrix) = ones(size(ivreftomatrix));                                  
                                end
                                
                            end
                        end
                        
                    end
                end
                  
            end
        end    
        
    end
       
end

%---------------------------------------------------------------
%
%			edgealignedwithedge
%
%---------------------------------------------------------------

edgealignedwithedge = [];

% Changed 050504 to specorder >= 2
% if specorder >= 3
if specorder >= 2
    edgealignedwithedge = sparse(eye(nedges));
	if showtext >= 4
		disp('         checking which edges are aligned with other edges ...')
	end
    
    listofactiveedges = [1:nedges].';
    listofactiveedges(edgedata.offedges) = [];
    vecstocheck = edgenormvecs(listofactiveedges,:);
    iv = find(vecstocheck(:,1)<0);
    vecstocheck(iv,:) = -vecstocheck(iv,:);
    iv = find(vecstocheck(:,1)==0 & vecstocheck(:,2)<0);
    vecstocheck(iv,:) = -vecstocheck(iv,:);
    iv = find(vecstocheck(:,1)==0 & vecstocheck(:,2)==0 & vecstocheck(:,3)<0);
    vecstocheck(iv,:) = -vecstocheck(iv,:);
    
    % Changed 12 Nov. 2014:
    % round off to 6 decimals
    
    vecstocheck = round(vecstocheck*1e6)/1e6;
    
    if ~isempty(listofactiveedges)
        [~,~,JJ] = unique(vecstocheck,'rows');
        N = histc(JJ,[1:length(listofactiveedges)]);
        iv = find(N>1);
        for ii = 1:length(iv)
            iv2 = JJ==iv(ii);
            edgestocheck = listofactiveedges(iv2);
            ntocheck = length(edgestocheck);
            cornernumberstocheck = edgedata.edgecorners(edgestocheck,:);
%             cornernumberstocheck = sort(cornernumberstocheck.').';
%             expandedcornermatrix = [reshape(repmat(cornernumberstocheck(:,1).',[ntocheck 1]),ntocheck^2,1) repmat(cornernumberstocheck(:,2),[ntocheck 1])]
            expandedstartcornermatrix = [reshape(repmat(cornernumberstocheck(:,1).',[ntocheck 1]),ntocheck^2,1) repmat(cornernumberstocheck(:,1),[ntocheck 1])];
            expandedendcornermatrix   = [reshape(repmat(cornernumberstocheck(:,2).',[ntocheck 1]),ntocheck^2,1) repmat(cornernumberstocheck(:,2),[ntocheck 1])];
            toexpandededgestocheck   =         repmat(edgestocheck,  [ntocheck 1]);
            fromexpandededgestocheck = reshape(repmat(edgestocheck.',[ntocheck 1]),ntocheck^2,1);
            
            % Now we don't need to check edges against themselves, or edge1
            % vs edge2 if edge2 vs edge1 has already been tested.
            
            indmat = reshape([1:ntocheck^2],ntocheck,ntocheck);
            indmat = triu(indmat,1);
            indmat = indmat(indmat);
            
%             expandedcornermatrix = expandedcornermatrix(indmat,:);
             expandedstartcornermatrix = expandedstartcornermatrix(indmat,:);
             expandedendcornermatrix   = expandedendcornermatrix(indmat,:);
            fromexpandededgestocheck = fromexpandededgestocheck(indmat,:);
            toexpandededgestocheck = toexpandededgestocheck(indmat,:);
            
            % We make the easiest check first: if two edges with the same
            % direction vector share one corner, they must be aligned.

            A = [expandedstartcornermatrix expandedendcornermatrix];
            A = sort(A.').';
            A = diff(A.').';
            iv3 = find(sum(A.'==0).');            
            validcombs = [fromexpandededgestocheck(iv3) toexpandededgestocheck(iv3)];
            if ~isempty(validcombs)
                indexvals = (validcombs(:,2)-1)*nedges + validcombs(:,1);
                edgealignedwithedge(indexvals) =  edgealignedwithedge(indexvals) + 1;     
                indexvals = (validcombs(:,1)-1)*nedges + validcombs(:,2);
                edgealignedwithedge(indexvals) =  edgealignedwithedge(indexvals) + 1;
                fromexpandededgestocheck(iv3) = [];
                toexpandededgestocheck(iv3) = [];
                expandedstartcornermatrix(iv3,:) = [];
                expandedendcornermatrix(iv3,:) = [];
            end
            
            % For the other edge-pair combinations we must check if the two
            % edges are really aligned. We do this by checking if two
            % vectors are aligned: edge1-vector and the vector from
            % edge1start to edge2start.            
                        
            direcvec = planedata.corners(expandedstartcornermatrix(:,2),:) - planedata.corners(expandedstartcornermatrix(:,1),:);
            direcveclengths = sqrt( sum(direcvec.'.^2) ).';
            direcvec = direcvec./(direcveclengths(:,ones(1,3)) + eps*10);
            test = abs(sum( (direcvec.*edgenormvecs(fromexpandededgestocheck,:)).' ));
            iv3 = find(abs(test-1)< 1e-8);
            validcombs = [fromexpandededgestocheck(iv3) toexpandededgestocheck(iv3)];
            if ~isempty(validcombs)
                indexvals = (validcombs(:,2)-1)*nedges + validcombs(:,1);
                edgealignedwithedge(indexvals) =  edgealignedwithedge(indexvals) + 1;     
                indexvals = (validcombs(:,1)-1)*nedges + validcombs(:,2);
                edgealignedwithedge(indexvals) =  edgealignedwithedge(indexvals) + 1;
                fromexpandededgestocheck(iv3) = [];
                toexpandededgestocheck(iv3) = [];
                expandedstartcornermatrix(iv3,:) = [];
                expandedendcornermatrix(iv3,:) = [];
            end            
        end
    end
    
end

%---------------------------------------------------------------
%
%			edgeseesedge
%
%---------------------------------------------------------------
%
% Make a matrix of which edges can see which edges
%
% Planes that are in front of the edge-planes might have visible edges.
% If closwedang = 0, then all other edges are potentially visible.
% If closwedang < pi, then all planes in front of plane1 or plane2 are OK
% If closwedang > pi, then all planes in front of plane1 and plane2 are OK

if showtext >= 4
	disp('         checking which edges see which edges...')
end

edgeseesedge = int8(zeros(nedges,nedges));
for ii = 1:nedges
	closwed = edgedata.closwedangvec(ii);
	if closwed == 0
		edgeseesedge(:,ii) = ones(nedges,1);
	else
		plane1 = edgedata.planesatedge(ii,1);
		plane2 = edgedata.planesatedge(ii,2);	
		okplanelist1 = find(planedata.planeseesplane(:,plane1)==1);
		okplanelist2 = find(planedata.planeseesplane(:,plane2)==1);		
		
		if closwed < pi
			okplanes = [okplanelist1;okplanelist2];
			okplanes = sort(okplanes);
		else
			okplanes = zerosvec1;
			okplanes(okplanelist1) = okplanes(okplanelist1)+1;
			okplanes(okplanelist2) = okplanes(okplanelist2)+1;
			okplanes = find(okplanes==2);
		end
		
		okedges = edgedata.edgesatplane(okplanes,:);
		[n1,n2] = size(okedges);
		okedges = reshape(okedges,n1*n2,1);
		
		if ~isempty(okedges)
			okedges = okedges(okedges~=0);			
			edgeseesedge(okedges,ii) = int8(double(edgeseesedge(okedges,ii)) + 1); 
		end

		edgesatplane1 = edgedata.edgesatplane(plane1,:);
		edgesatplane2 = edgedata.edgesatplane(plane2,:);
		okedges = [edgesatplane1 edgesatplane2].';
		okedges = okedges(okedges~=0);
		edgeseesedge(okedges,ii) = int8(double(edgeseesedge(okedges,ii)) + 1); 
	
	end	   	% else ..... if closwed == 0

end % for ii = 1:nedges,

% Fixed 4 Oct. 2017: sometimes edgealignedwithedge is empty
if ~isempty(edgealignedwithedge)
    % Mask out all edge pairs that are aligned
    edgeseesedge = double(edgeseesedge) - double(edgealignedwithedge);
end

% Make sure the edge-to-same-edge combos are masked out
% 	edgeseesedge = sign(edgeseesedge.*(edgeseesedge.')).*(1-eye(nedges));
edgeseesedge = int8( ( (edgeseesedge & (edgeseesedge.'))>0).*(1-eye(nedges)));

% Mask out all edges that should be switched off

mask = ones(nedges,nedges);
mask(edgedata.offedges,:) = mask(edgedata.offedges,:)*0;
mask(:,edgedata.offedges.') = mask(:,edgedata.offedges.')*0;

edgeseesedge = edgeseesedge & mask;

clear mask	

% For non-convex planes, there might be some edges that can never see each
% other.

% NB!!!!! Double check if the test below should be reversed for interior problems!!!
% PS 2011-05-27

ivplaneswindents = find(planedata.planehasindents);
if ~isempty(ivplaneswindents)
   uppertri_planeseesplane = triu(planedata.planeseesplane);
   for ii = 1:length(ivplaneswindents)
        indpl = ivplaneswindents(ii);
        disp(['Plane ',int2str(indpl),' has ind'])
        edgelist = edgesatplane(indpl,:);
        connplanes = planesatedge(edgelist,:);
        connplanes = reshape(connplanes,size(connplanes,1)*size(connplanes,2),1);
        connplanes = setdiff(connplanes,indpl);
        for jj = 1:length(connplanes)

            % First find connected plane-pairs that see each other. Their
            % edges can not see each other.
            planepairs = connplanes(uppertri_planeseesplane(connplanes(jj),connplanes));
            if ~isempty(planepairs)
                for kk = 1:length(planepairs)
                    edge1 = intersect(edgesatplane(connplanes(jj),:),edgesatplane(indpl,:));
                    edge2 = intersect(edgesatplane(planepairs(kk),:),edgesatplane(indpl,:));
                    edgeseesedge(edge1,edge2) = 0;
                    edgeseesedge(edge2,edge1) = 0;
                    
                end
            end            
        end
        
        % Second, find connected plane-pairs that have plane normal
        % vectors in the same direction. Their edges can never see each
        % other.
        A1 = planenvecs(connplanes,:);
        [A2,A2sort] = sortrows(A1);
        findnodiffs = sum(abs(diff(A2)).').';
        if any(findnodiffs==0)
            ivnodiffs = find(findnodiffs==0);
            ivnodiffs = ivnodiffs(:);
            planepairs = connplanes([A2sort(ivnodiffs) A2sort(ivnodiffs+1)]);
            if ~isempty(planepairs)
                for kk = 1:length(planepairs)
                    edge1 = intersect(edgesatplane(planepairs(kk,1),:),edgesatplane(indpl,:));
                    edge2 = intersect(edgesatplane(planepairs(kk,2),:),edgesatplane(indpl,:));
                    edgeseesedge(edge1,edge2) = 0;
                    edgeseesedge(edge2,edge1) = 0;
                    
                end
            end            
%             findnodiffs = connplanes(A2sort(findnodiffs))
        end    
   end
end

%-----------------------------------------------------------
% Check the planes with indents: some edges can not see each other
%
% 12 Nov. 2014

iv = find(planedata.planehasindents);

for ii = 1:length(iv)
   jj = iv(ii);
   disp(['Plane with indents: ',int2str(jj)])
   
   % Find the in-plane edge normal vectors
   edgenormvecs = zeros(nedges,3);
   
   for kk = 1:ncornersperplanevec(jj)
      co1 = planecorners(jj,kk);
      if kk+1<=ncornersperplanevec(jj)
          co2 = planecorners(jj,kk+1);
      else
         co2 = planecorners(jj,1); 
      end
      alongedgevec = corners(co2,:)-corners(co1,:);
      iv2 = find(edgecorners(:,1)==co1 & edgecorners(:,2)==co2);
      if isempty(iv2)
          iv2 = find(edgecorners(:,1)==co2 & edgecorners(:,2)==co1);          
      end
      alongedgevec = alongedgevec/norm(alongedgevec);      
      edgenormvecs(iv2,:) = cross(-alongedgevec,planenvecs(jj,:));      
   end

   % Now we check if two edges, A and B, in each edge pair are behind
   % or in front of each other (We check the two edge points (a little off the endpoints))
   % If (A is in front of B) & (B is behind A),then
   %    they can never see each other. 
   % If (A is behind B) & (B is behind A) & (wedge angle for A > pi) &
   %    (wedge angle for B > pi) & no obstruction (from midpoint to
   %    midpoint).
   % If (A is in front of B) & (B is in front of A) & no obstruction, then
   %    they can see each other.
    
   edgestartcoordsnudgemod = edgestartcoords + (edgestartcoordsnudge-edgestartcoords)*1e6;
   edgeendcoordsnudgemod   = edgeendcoords   + (edgeendcoordsnudge-edgeendcoords)*1e6;
   
   listofedges = sort(edgesatplane(jj,:));
   nedges_at_plane = length(listofedges);
   edgetoedgerelation = zeros(nedges_at_plane,nedges_at_plane);
   for kk = 1:nedges_at_plane
      edgetocheck = listofedges(kk);
      disp(['   Edge number ',int2str(edgetocheck)])
      cornerstocheck = edgecorners(listofedges,1);
      pointinfront1 = EDinfrontofplane(edgestartcoordsnudgemod(listofedges,:),edgenormvecs(edgetocheck,:),corners(edgecorners(edgetocheck,1),:),corners(edgecorners(edgetocheck,2),:));
      pointinfront2 = EDinfrontofplane(edgeendcoordsnudgemod(listofedges,:),edgenormvecs(edgetocheck,:),corners(edgecorners(edgetocheck,1),:),corners(edgecorners(edgetocheck,2),:));

      edgetoedgerelation(kk,:) = sign(pointinfront1+pointinfront2).';
   end
         
   edgescanobstruct = listofedges(any(edgetoedgerelation.'==-1));
   
   % In edgetoedgerelation, a 0 means two edges can never see each other
   % a 1 means that they can see each other if there is no obstruction
   % and a -1 means that they can see each other if the two wedge angles
   % are larger than pi.
   combinededgetoedgerelation = sign(edgetoedgerelation + edgetoedgerelation.');
   
   % Step 1: check all pairs with -1
   for kk = 1:nedges_at_plane
       edge1 = listofedges(kk);
       iv = find(combinededgetoedgerelation(kk,:)==-1);
       for ll=iv
           edge2 = listofedges(ll);
          if closwedangvec(edge1)<pi && closwedangvec(edge2)<pi
              combinededgetoedgerelation(kk,ll) = 1;
              combinededgetoedgerelation(ll,kk) = 1;              
          end           
       end
   end

   % Step 2: for all edgepairs with edgetoedgerelation = 1, check if
   % edgemidpoint-to-edgemidpoint crosses any of the other edges. If yes:
   % obstruction.
   %
   % NB! Only very few edges can potentially obstruct! Many of the edges
   % will have all other edges behind themselves (or aligned with themselves). 
   % They can not obstruct.
   % 
%    combinededgetoedgerelation
   
   for kk = 1:nedges_at_plane
       edge1 = listofedges(kk);
       iv = find(combinededgetoedgerelation(kk,:)==1);
       for ll=iv
           edge2 = listofedges(ll);
           
           edgestocheck = setdiff(edgescanobstruct,[edge1 edge2]);
           
%            disp(['Will check ',int2str(edge1),' to ',int2str(edge2),' against'])
           for edgeinbetween = edgestocheck
                mm = find(listofedges == edgeinbetween);
%                 disp(['      Edgeinbetween = ',int2str(edgeinbetween)])                
                if edgetoedgerelation(mm,kk) == edgetoedgerelation(mm,ll)
%                    disp(['Edges ',int2str(edge1),' and ',int2str(edge2),' are on the same side']) 
                   combinededgetoedgerelation(kk,ll) = 2;
                   combinededgetoedgerelation(ll,kk) = 2;
                else
                   disp(['      Need to check edge ',int2str(edge1),' and ',int2str(edge2),' vs edge inbetween: ',int2str(edgeinbetween)]) 
                end
           
           end
       end
   end
      
   
end



%-----------------------------------------------------------
% Go through all edgeseesedge combinations. If two edges see
% each other without belonging to the same plane, there should
% be a gain factor of 2.
%
% (12.11.2014) Also if they belong to a plane with indents, they could 
%              have a gain factor of 2.

if showtext >= 4
	disp('         checking which edge-to-edge paths that run along planes')
end

edgetoedgestrength = 2*edgeseesedge;
for ii = 1:nedges
	if sum(edgetoedgestrength(:,ii)) ~= 0
		plane1 = edgedata.planesatedge(ii,1);
		plane2 = edgedata.planesatedge(ii,2);	
		edgesatplane1 = edgedata.edgesatplane(plane1,:);
		edgesatplane2 = edgedata.edgesatplane(plane2,:);
		sameplaneedges = [edgesatplane1 edgesatplane2].';
		sameplaneedges = sameplaneedges(sameplaneedges~=0);
		edgetoedgestrength(sameplaneedges,ii) = edgetoedgestrength(sameplaneedges,ii)*0 + 1;  
	end % % if sum(edge
end

% an edge can not contribute to itself

edgetoedgestrength = edgetoedgestrength.*(1-eye(nedges));

% The edges belonging to the same planes, but that are inactive (90 degrees etc)
% must be switched off here too.
% 011021 This is redundant.
 edgetoedgestrength = edgetoedgestrength.*(edgeseesedge>0);
edgetoedgestrength = int8(edgetoedgestrength);

%----------------------------------------------------------------------------
%
%		CHECK THIN FOLDED PLANES
%
%----------------------------------------------------------------------------

% Check those edges that have identical corner pairs

rowdiff = [edgedata.edgecorners(2:end,1)-edgedata.edgecorners(1:end-1,1) edgedata.edgecorners(2:end,2)-edgedata.edgecorners(1:end-1,2)];
ivfoundedgepairs = find( sum(abs(rowdiff.')).' == 0);

for ii = 1:length(ivfoundedgepairs)
    edge1 = ivfoundedgepairs(ii);
    edge2 = edge1+1;
    edgetoedgestrength(edge1,edge2) = 0;
    edgetoedgestrength(edge2,edge1) = 0;    
end

%----------------------------------------------------------------------------
%
%		CHECK OBSTRUCTION OF EDGE-TO-EDGE PATHS
%
%----------------------------------------------------------------------------
%
% We construct two long lists: 'from-edges' and 'to-edges' for which
% obstruction should be tested. These lists are made up of all the
% combinations in edgeseesedge that have a value 1. 
% NB! We use symmetry - if edge 1 can see edge 2, then edge 2 can also see
% edge 1.
%
% Both edges should be subdivided into nedgesubs segments. So, the long
% lists must be expanded so that each from-edge/to-edge combination is
% replaced by nedgesubs^2 as many.
% For efficiency, we make a shortlist of the actual edge subdivisions and
% we then expand long lists with pointers to this shortlist.
%
% The final matrix, edgeseespartialedge, will have an integer value which
% tells whether two edges see each other:
%
%   0               Edges can not see each other, or one of the edges is inactive
%                   (an inactive edge has an integer wedge index, or has one plane
%                   which is TOTABS)
%   2^nedgesubs-1   Two edges see each other completely.
%   Other integers  Two edges see each other partly.
%   Neg. integer    As above, but in addition, the edges do not belong
%                   to the same plane.
%
% We check in the visedgesfroms and visedgesfromr lists to check which
% combinations we don't need to check.

if showtext >= 3
	disp(['   Checking for obstructing planes between edges and edges'])
    disp(' ')
    disp('NB!!! The value of nedgesubs is temporarily set to 1 for the edge-to-edge vis. test')
end

obstructtestneeded = (sum(planedata.canplaneobstruct)~=0 && strcmp(planedata.modeltype,'convex_ext')==0 & strcmp(planedata.modeltype,'singleplate')==0);

maxedgetoedgevisvalue = 2^(2*nedgesubs)-1;

if obstructtestneeded
        
    nedgesubsorig = nedgesubs;
    nedgesubs = 1;
    edgeseesedge = triu(edgeseesedge);
    
    iv = full(find(edgeseesedge~=0));
    
    if ~isempty(iv)
        fromedges = floor(iv/nedges)+1;
        toedges = iv - (fromedges-1)*nedges;
        ncombs = length(fromedges);
        
        ivcancel = find(ismember(fromedges,edgesnottoworryabout));
        fromedges(ivcancel) = [];
        toedges(ivcancel) = [];

        ivcancel = find(ismember(toedges,edgesnottoworryabout));        
        toedges(ivcancel) = [];
        fromedges(ivcancel) = [];
        clear ivcancel

        % Some of these fromedges-toedges combos might involve two
        % neighboring edges that form an indent. Those should be
        % removed before we start checking for obstructions.
        %
        % We check the [fromedges toedges] matrix against the
        % indentingedgepairs matrix.

        % Change 28.10.2013: empty indentingedgepairs caused problems
        if ~isempty(edgedata.indentingedgepairs)
            [~,loc]= ismember(edgedata.indentingedgepairs,sort([fromedges toedges].').','rows');

            if ~isempty(loc)
                ivmatch = find(loc);
               fromedges(loc(ivmatch)) = [];
               toedges(loc(ivmatch)) = [];
            end
            ncombs = length(fromedges);
        end


        % The previous version did not use 'rows' - but that must be wrong?
        % No! Because we want to subdivide each edge only once ,not both
        % edges in a pair!
         [uniqueedges,~,junique] = unique([fromedges toedges]);
                  
%         [uniqueedges,iunique,junique] = unique([fromedges toedges],'rows')
%         pause
        clear fromedges toedges
        
        % The lists edgesubcoords are the shortlists.
        
        [edgesubcoords,edgeweightlist,edgenumberlist] = EDgetedgepoints(...
            edgedata.edgestartcoords(uniqueedges,:),edgedata.edgeendcoords(uniqueedges,:),nedgesubs,1);

        % The two lists below contain pointers to the first edge segment in the
        % shortlist.
	
        fromcoords_refp1 = nedgesubs*(junique(1:ncombs)-1)+1;
        tocoords_refp1    = nedgesubs*(junique(ncombs+1:2*ncombs)-1)+1;
        clear junique
        
        % Now the original [fromedges toedges] could be recovered by:
        % uniqueedges([edgenumberlist(fromcoords_refp1) edgenumberlist(tocoords_refp1)])
        % reshape(uniqueedges([edgenumberlist(fromcoords_refp1) edgenumberlist(tocoords_refp1)]),6,2)
%         uniqueedges([edgenumberlist(fromcoords_refp1) edgenumberlist(tocoords_refp1)]);
%         npairs = length(uniqueedges([edgenumberlist(fromcoords_refp1) edgenumberlist(tocoords_refp1)]));
%         reshape(uniqueedges([edgenumberlist(fromcoords_refp1) edgenumberlist(tocoords_refp1)]),npairs/2,2)

        addmask = [0:nedgesubs-1];
        onesvec1 = ones(1,nedgesubs);
        ntot1 = ncombs*nedgesubs;
	
        % First, the from-list gets repetitions of the same values
        expandfrom_ref = fromcoords_refp1(:,onesvec1);
        clear fromcoords_refp1;
        expandfrom_ref = reshape(expandfrom_ref.',ntot1,1);
        
        % Second, the expanded from-list gets repetitions that are
        % incremented in steps of 1 (because the shortlists have the
        % edge segments placed after each other).
        expandfrom_ref = expandfrom_ref(:,onesvec1);
        expandfrom_ref = expandfrom_ref + addmask(ones(ntot1,1),:);
        expandfrom_ref = reshape(expandfrom_ref.',ncombs*nedgesubs^2,1);
        
        % Same thing for the to-list except that the first expansion
        % gives repetitions that are incremented in steps of 1.
        expandto_ref   = tocoords_refp1(:,onesvec1);
        clear tocoords_refp1;
        expandto_ref = expandto_ref + addmask(ones(ncombs,1),:);
        expandto_ref = reshape(expandto_ref.',ntot1,1);
	
        % Second expansion for the to-list: repetitions
        expandto_ref = expandto_ref(:,onesvec1);
        expandto_ref = reshape(expandto_ref.',ncombs*nedgesubs^2,1);
	
        % Now we have expanded lists of all combinations, and we can 
        % make lists with the coordinates, the edge numbers, and the weights of
        % the segments.
        
        % edgesubcoords(expandfrom_ref,:) will give all the starting points,
        % "fromcoords"
        % edgesubcoords(expandto_ref,:) will give all the ending points,
        % "tocoords"
        % uniqueedges(edgenumberlist(expandfrom_ref)) will give all
        % the starting edges, "fromedges"
        % uniqueedges(edgenumberlist(expandto_ref)) will give all
        % the ending edges, "toedges"

% [uniqueedges(edgenumberlist(expandfrom_ref))   uniqueedges(edgenumberlist(expandto_ref))]
% pause

        % We can make a simple test of which planes that could be excluded
        % from the obstruction test by checking the edgeseesplane for the
        % relevant edges. It can be multiplied by the canplaneobstruct list
        % when calling the findobstructedpaths function.
        isplaneactiveforobstruction = (sum(edgedata.edgeseesplane(:,uniqueedges).'==1)>0);
        
        
% %         save /Users/petersve/Documents/Temp/Temp2.mat
% %         pause

        global STARTPLANES ENDPLANES
        global BIGFROMCOORDS BIGTOCOORDS BIGSTARTPLANES BIGENDPLANES
% The global definition was missing - why?
        global REFTOFROMCOSHO REFTOTOCOSHO
        
%             FROMCOORDS = edgesubcoords(expandfrom_ref,:);
        FROMCOORDSSHORTLIST = edgesubcoords;
        REFTOFROMCOSHO = uint32(expandfrom_ref);
%         TOCOORDS   = edgesubcoords(expandto_ref,:);
        TOCOORDSSHORTLIST   = edgesubcoords;
        REFTOTOCOSHO = uint32(expandto_ref);
        if nedges < 256
            bigfromedge   = uint8(uniqueedges(edgenumberlist(expandfrom_ref)).');
        else
            bigfromedge   = uint16(uniqueedges(edgenumberlist(expandfrom_ref)).');
        end
        clear expandfrom_ref
        bigfromedge = bigfromedge(:);
        if nedges < 256
            bigtoedge     = uint8(uniqueedges(edgenumberlist(expandto_ref)).');
        else
            bigtoedge     = uint16(uniqueedges(edgenumberlist(expandto_ref)).');
        end
        clear expandto_ref
        bigtoedge = bigtoedge(:);
%         [bigfromedge bigtoedge]
%         pause
        %  bigfromweight = edgeweightlist((expandfrom_ref));
      %  bigtoweight   = edgeweightlist((expandto_ref));
        
        STARTPLANES = [edgedata.planesatedge(bigfromedge,1) edgedata.planesatedge(bigfromedge,2)];
        ENDPLANES   = [edgedata.planesatedge(bigtoedge,1) edgedata.planesatedge(bigtoedge,2)];        
        
%                     REFTOFROMCOSHO = REFTOFROMCOSHO(nonobstructedpaths);
%                     REFTOTOCOSHO = REFTOTOCOSHO(nonobstructedpaths);
%                     STARTPLANES = STARTPLANES(nonobstructedpaths,:);
%                     ENDPLANES = ENDPLANES(nonobstructedpaths,:);
%                     bigfromedge = bigfromedge(nonobstructedpaths);
%                     bigtoedge = bigtoedge(nonobstructedpaths);
%                     npathstocheck = size(REFTOFROMCOSHO,1);

        shouldplanebechecked = double(planedata.canplaneobstruct).*double(isplaneactiveforobstruction);
        nplanesperbatch = ceil(sum(shouldplanebechecked)/ndiff2batches);
        if nplanesperbatch > 0
            temp = cumsum(shouldplanebechecked);
            diff2batchlist = zeros(ndiff2batches,2);
            diff2batchlist(1,1) = 1;
            loopisfinished = 0;
            for ii = 1:ndiff2batches-1
                ivtemp = find(temp==nplanesperbatch*ii);
                if ~isempty(ivtemp)
                    diff2batchlist(ii,2) = ivtemp(1);
                    diff2batchlist(ii+1,1) = ivtemp(1)+1;
               else
                   loopisfinished = 1;
               end               
            end
            diff2batchlist(ii,2) = nplanes;
            ivtemp = diff2batchlist(:,1)>0;
            diff2batchlist = diff2batchlist(ivtemp,:);
            ndiff2batches = size(diff2batchlist,1);

            npathstocheck = size(REFTOFROMCOSHO,1);
            for ii = 1:ndiff2batches
                if npathstocheck > 0
                    maskedchkpla = shouldplanebechecked;
                    if ii > 1
                        maskedchkpla(1:diff2batchlist(ii,1)-1) = maskedchkpla(1:diff2batchlist(ii,1)-1)*0;    
                    end
                    if ii < ndiff2batches
                        maskedchkpla(diff2batchlist(ii,2)+1:end) = maskedchkpla(diff2batchlist(ii,2)+1:end)*0;                    
                    end
                    nplanestocheck = sum(maskedchkpla);
                    if showtext >= 3 
                        if round(ii/1)*1==ii
                            disp(['   Batch ',int2str(ii),' of ',int2str(ndiff2batches),': '])
                            disp(['   ',int2str(npathstocheck),' paths to check, ',int2str(nplanestocheck),' planes to check obstruction for'])
%                         disp(['   plane ',int2str(find(maskedchkpla))])
                        end
                    end
                
                    % We create the BIGFROMCOORDS matrices and BIGTOCOORDS matrices
                    % outside since they need to be exactly identical from batch to
                    % batch as long as nplanestocheck is the same.
        
                    npathstocheck = size(REFTOFROMCOSHO,1);
                
                    ntot = nplanestocheck*npathstocheck;
                
                    BIGFROMCOORDS = reshape(repmat(FROMCOORDSSHORTLIST(REFTOFROMCOSHO,:).',[nplanestocheck,1]),3,ntot).';
                    BIGTOCOORDS = reshape(repmat(TOCOORDSSHORTLIST(REFTOTOCOSHO,:).',[nplanestocheck,1]),3,ntot).';
                    BIGSTARTPLANES = reshape(repmat(STARTPLANES.',[nplanestocheck,1]),2,ntot).';
                    BIGENDPLANES = reshape(repmat(ENDPLANES.',[nplanestocheck,1]),2,ntot).';
                
%                     [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_edgetoedge(maskedchkpla,planedata.planeseesplane,...
%                         planedata.planeeqs,planedata.planeeqs(:,1:3),planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.rearsideplane);
                    [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_edgetoedge(maskedchkpla,...
                        planedata.planeeqs,planedata.planeeqs(:,1:3),planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,planedata.ncornersperplanevec,planedata.rearsideplane);
                                    
                    if ~isempty(edgehits)
                           % Look for the edge-to-edge combos that share a
                           % plane. When such combos have hit an edge, it
                           % must mean that they are planes with indents.
                          wasedgehitnotindent = prod(double(diff(sort([STARTPLANES(edgehits,:) ENDPLANES(edgehits,:)].')))).';
                        indentedgehits = find(wasedgehitnotindent==0);
                        if ~isempty(indentedgehits)
                          nonobstructedpaths = setdiff(nonobstructedpaths,edgehits(indentedgehits));
                        end
                        
                    end
                    
                    if ~isempty(nonobstructedpaths)
                        REFTOFROMCOSHO = REFTOFROMCOSHO(nonobstructedpaths);
                        REFTOTOCOSHO = REFTOTOCOSHO(nonobstructedpaths);
                        STARTPLANES = STARTPLANES(nonobstructedpaths,:);
                        ENDPLANES = ENDPLANES(nonobstructedpaths,:);
                        bigfromedge = bigfromedge(nonobstructedpaths);
                        bigtoedge = bigtoedge(nonobstructedpaths);
                        npathstocheck = size(REFTOFROMCOSHO,1);
                    else
                        REFTOFROMCOSHO = [];
                        REFTOTOCOSHO = [];
                        STARTPLANES = [];
                        ENDPLANES = [];
                        bigfromedge = [];
                        bigtoedge = [];
                        npathstocheck = 0;
                    end
                   
                end   %     if npathstocheck > 0,
            
            end           %  for ii = 1:ndiff2batches,
        else
            nonobstructedpaths = 1;
            
        end         %if nplanesperbatch > 0,

        
            clear REFTOFROMCOSHO REFTOTOCOSHO STARTPLANES ENDPLANES


        if ~isempty(nonobstructedpaths)
            
            % Add the segment pieces together
            % NB! We don't try to implement the correct way to do it since
            % we would need nedgesubs^2 bits to represent all edge-segment
            % to edge-segment visibility possibilities.
            % Instead we just establish whether the two edges see each
            % other fully or partly.

            visiblesegmentscounter = ones(size(bigfromedge));
            test = [bigfromedge bigtoedge];
            
            ncombs = length(bigfromedge);
            dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
            ivremove = find(dtest==1);
            
            while ~isempty(ivremove)
                visiblesegmentscounter(ivremove+1) = visiblesegmentscounter(ivremove+1) + visiblesegmentscounter(ivremove);
                visiblesegmentscounter(ivremove) = [];
                bigfromedge(ivremove) = [];
                bigtoedge(ivremove,:) = [];
                
                test = [bigfromedge bigtoedge];
                ncombs = length(bigfromedge);
                dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                ivremove = find(dtest==1);
	
            end
            
            indexvec = uint32((double(bigfromedge)-1)*nedges + double(bigtoedge));
            if maxedgetoedgevisvalue < 128
                edgeseespartialedge = int8(zeros(nedges,nedges));            
            elseif maxedgetoedgevisvalue < 32768
                edgeseespartialedge = int16(zeros(nedges,nedges));            
            else
                edgeseespartialedge = int32(zeros(nedges,nedges));            
            end                
                
            iv = find(visiblesegmentscounter==nedgesubs^2);
            edgeseespartialedge(indexvec(iv)) = maxedgetoedgevisvalue*ones(size(iv));
            iv = find(visiblesegmentscounter>0 & visiblesegmentscounter<nedgesubs^2);
            edgeseespartialedge(indexvec(iv)) = maxedgetoedgevisvalue/2*ones(size(iv));
            
            if maxedgetoedgevisvalue < 128
               edgeseespartialedge = int8(double(edgeseespartialedge) + double(triu(edgeseespartialedge)).');
           elseif maxedgetoedgevisvalue < 32768
               edgeseespartialedge = int16(double(edgeseespartialedge) + double(triu(edgeseespartialedge)).');
           else
               edgeseespartialedge = int32(double(edgeseespartialedge) + double(triu(edgeseespartialedge)).');               
           end
        else
            edgeseespartialedge = sparse(zeros(nedges,nedges));            
        end
    else
        edgeseespartialedge = [];
    end
else
    if maxedgetoedgevisvalue < 128
        edgeseespartialedge = int8((edgeseesedge*maxedgetoedgevisvalue));            
    elseif maxedgetoedgevisvalue < 32768
        edgeseespartialedge = int16((edgeseesedge*maxedgetoedgevisvalue));            
    else
        edgeseespartialedge = int32((edgeseesedge*maxedgetoedgevisvalue));            
    end                
%     edgeseespartialedge = sparse(edgeseesedge*maxedgetoedgevisvalue);
    if showtext >= 2
        disp('      No obstruction tests needed!')
    end
end   %(obstructtestneeded)


% NB!!!! 27 May 2011: the part below looks suspect. Probable error:
% expression for multfactor is probably wrong if allthinplanes and
% allplanesinsameplane are both -1.
%
% For the special case of thin planes that are all in the same plane, only
% edge pairs that belong to the same plane should be included. That is,
% when edgeseespartialedge is negative, then we should set
% edgeseespartialedge to zero.

% allthinplanes = sign(1-sum(~planeisthin));
% 
% allplanesinsameplane = 1 - sum(sum(planeseesplane==1));
% 
% multfactor = 1 - allthinplanes*allplanesinsameplane
% 
% iv = find(edgetoedgestrength==2)
% if ~isempty(iv),
%     if maxedgetoedgevisvalue < 128,
%         edgeseespartialedge(iv) = int8(-double(edgeseespartialedge(iv))*multfactor);            
%     elseif maxedgetoedgevisvalue < 32768,
%         edgeseespartialedge(iv) = int16(-double(edgeseespartialedge(iv))*multfactor);            
%     else
%         edgeseespartialedge(iv) = int32(-double(edgeseespartialedge(iv))*multfactor);            
%     end                
% end
% % % % % edgeseespartialedge(iv) = -double(edgeseespartialedge(iv));

% 13 March 2013: Reintroduced one part: when edgetoedgestrength = 0, then
% edgeseespartialedge must also be zero

edgeseespartialedge = double(edgeseespartialedge).*double(edgetoedgestrength>0);

           if maxedgetoedgevisvalue < 128
                edgeseespartialedge = int8(edgeseespartialedge);            
            elseif maxedgetoedgevisvalue < 32768
                edgeseespartialedge = int16(edgeseespartialedge);            
            else
                edgeseespartialedge = int32(edgeseespartialedge);            
           end                
 
            % 21 October 2017 bug: the edgetoedgestrength was not used, so
            % it seemed like all edge-to-edge paths ran along planes.
            
iv = find(edgetoedgestrength==2);
edgeseespartialedge(iv) = -edgeseespartialedge(iv);

clear edgeseesedge edgetoedgestrength

%----------------------------------------------------------------------------
%
%		CYLINDRICAL EDGE-TO-EDGE PARAMETERS
%
%----------------------------------------------------------------------------
%
% For each edge, calculate the cylindrical coordinates of the starting
% and ending points of all other edges relative to the first edge.

if difforder >= 2 && ~isempty(edgeseespartialedge)
	
	if showtext >= 4
		disp('         E2E')
	end
	
    zerosvec4 = zeros(nedges,nedges);
	Bigre1 = (zerosvec4);
	Bigthetae1 = (zerosvec4);
	Bigze1 = (zerosvec4);
	Bigre2 = (zerosvec4);
	Bigthetae2 = (zerosvec4);
	Bigze2 = (zerosvec4);
    clear zerosvec4
    
	for edge1 = 1:nedges
		if showtext >= 4
            if round(edge1/1)*1 == edge1
    			disp(['            Edge no. ',int2str(edge1)])
            end
        end
		edge1coords = [edgedata.edgestartcoords(edge1,:);edgedata.edgeendcoords(edge1,:)];
%%%		iv = find(edgeseesedge(edge1,:)==1).';
		iv = find(edgeseespartialedge(edge1,:)~=0).';        

		if ~isempty(iv)		
			% First find the subset of edges which belong to the same plane as the
			% edge itself. These should be treated separately for higher accuracy

			% iv1 will be the edges on the reference plane. They should have theta = 0.
			
            
			refplane = edgedata.planesatedge(edge1,1);
            ncornersperplanevec = double(planedata.ncornersperplanevec);
			iv1 = edgedata.edgesatplane( refplane,1:planedata.ncornersperplanevec(  refplane ));
%%%%			iv1 = strtrim( iv1( find(iv1~= edge1)).')  % This removed edge nr 32!!
			iv1 = iv1( iv1~= edge1).';

			edgestart = edgedata.edgestartcoordsnudge(iv1,:);
			edgeend   = edgedata.edgeendcoordsnudge(iv1,:);
	
% % % 			[rs,thetas,zs,slask,slask,slask,zw1,Atrans1,xneworigo1] = ...
% % % 			ESIE2coordtrans(edgestart,[],edge1coords,edgenvecs(edge1,:));
% % % 			Bigre1(iv1,edge1) = rs;
% % % 			Bigthetae1(iv1,edge1) = 0*iv1;
% % % 			Bigze1(iv1,edge1) = zs;
% % % 
% % % 			[rs,thetas,zs,slask,slask,slask,zw1,Atrans1,xneworigo1] = ...
% % % 			ESIE2coordtrans(edgeend,[],edge1coords,edgenvecs(edge1,:));
% % % 			Bigre2(iv1,edge1) = rs;
% % % 			Bigthetae2(iv1,edge1) = 0*iv1;
% % % 			Bigze2(iv1,edge1) = zs;

            [rs,~,zs,rr,~,zr] = ...
			EDcoordtrans2(edgestart,edgeend,edge1coords,edgedata.edgenvecs(edge1,:));
			Bigre1(iv1,edge1) = rs;
			Bigthetae1(iv1,edge1) = 0*iv1;
			Bigze1(iv1,edge1) = zs;
			Bigre2(iv1,edge1) = rr;
			Bigthetae2(iv1,edge1) = 0*iv1;
			Bigze2(iv1,edge1) = zr;

% % % 			sameedgesold = ESIE2findsame(iv,iv1);
            [samevalues,iivec,~] = intersect(iv,iv1);
% Bug found 050421!!
% % %             %             sameedges = [iivec jjvec];
% % % 			if sum(sum(sameedges)) ~= 0,
% % %             	iv( sameedges(1,:).' ) = [];
% % %             end
            if ~isempty(samevalues)
                iv(iivec) = [];
            end

            % iv2 will be the edges on the non-reference plane.
			% They should have theta = 2*pi - closthetavec.

			if ~isempty(iv)

				if edgedata.planesatedge(edge1,2) > 0
					secplane = edgedata.planesatedge(edge1,2);

					iv2 = edgedata.edgesatplane( secplane,1:ncornersperplanevec(secplane));
%%%					iv2 = strtrim( iv2( find(iv2~= edge1)).' );
					iv2 = iv2( iv2~= edge1).';

					if ~isempty(iv2)
						edgestart = edgedata.edgestartcoordsnudge(iv2,:);
						edgeend   = edgedata.edgeendcoordsnudge(iv2,:);

% % % 						[rs,thetas,zs,slask,slask,slask,zw1,Atrans1,xneworigo1] = ...
% % % 						ESIE2coordtrans(edgestart,[],edge1coords,edgenvecs(edge1,:));
% % % 						Bigre1(iv2,edge1) = rs;
% % % 						Bigthetae1(iv2,edge1) = (2*pi-closwedangvec(edge1))*ones(size(iv2));
% % % 						Bigze1(iv2,edge1) = zs;
% % % 
% % % 						[rs,thetas,zs,slask,slask,slask,zw1,Atrans1,xneworigo1] = ...
% % % 						ESIE2coordtrans(edgeend,[],edge1coords,edgenvecs(edge1,:));
% % % 						Bigre2(iv2,edge1) = rs;
% % % 						Bigthetae2(iv2,edge1) = (2*pi-closwedangvec(edge1))*ones(size(iv2));
% % % 						Bigze2(iv2,edge1) = zs;

                        [rs,~,zs,rr,~,zr] = ...
						EDcoordtrans2(edgestart,edgeend,edge1coords,edgedata.edgenvecs(edge1,:));
						Bigre1(iv2,edge1) = rs;
						Bigthetae1(iv2,edge1) = (2*pi-edgedata.closwedangvec(edge1))*ones(size(iv2));
						Bigze1(iv2,edge1) = zs;
						Bigre2(iv2,edge1) = rr;
						Bigthetae2(iv2,edge1) = (2*pi-edgedata.closwedangvec(edge1))*ones(size(iv2));
						Bigze2(iv2,edge1) = zr;

%						sameedges = ESIE2findsame(iv,iv2);
% Bug found 050421!!!
% % % %                         [samevalues,iivec,jjvec] = intersect(iv,iv1);
% % % %                         sameedges = [iivec;jjvec];
% % % % 						if sum(sum(sameedges)) ~= 0,
% % % % 							iv( sameedges(1,:).' ) = [];
% % % %                         end

                        [samevalues,iivec,~] = intersect(iv,iv2);                        
                        if ~isempty(samevalues)
                            iv(iivec) = [];
                        end

                        
					end
				end
			end
			if ~isempty(iv)

				% Move the edge coordinates to be checked a short distance away
	
				edgestart = edgestartcoordsnudge(iv,:);
				edgeend   = edgeendcoordsnudge(iv,:);
%				edgestart = edgestart + geomacc*(edgeend-edgestart);
%				edgeend   = edgeend   - geomacc*(edgeend-edgestart);

                [rs,thetas,zs,rr,thetar,zr] = ...
				EDcoordtrans2(edgestart,edgeend,edge1coords,edgedata.edgenvecs(edge1,:));							
				Bigre1(iv,edge1) = rs;
				Bigthetae1(iv,edge1) = thetas;
				Bigze1(iv,edge1) = zs;
				Bigre2(iv,edge1) = rr;
				Bigthetae2(iv,edge1) = thetar;
				Bigze2(iv,edge1) = zr;

			end
		end
    end

    %-----------------------------------------------------------
    % Go through all edges that are in-plane with each other, that is, two
    % edges have at least one plane each that are in-plane with each other.

    % First we identify all planes that have at least one more co-planar
    % plane
    
    planehascoplanar = any(planedata.planeseesplane == -1);
    
    % Then we go through the list of edges and every edge with a connected
    % plane that has another co-planar plane must also have potential
    % in-plane edges.

    % For each edge, we check if there are other edges (that don't belong to the same plane!!)
    % with thetaangle = 0 or thetaangle = 2*pi. If there are other such
    % edges, those edge-to-edge combinations should be shut off.
    %
    % After all such edge-to-edge combinations have been shut off, we still
    % need another pass to cancel edge-to-edge paths that pass entirely
    % across other edges.
    
    ivedges = 1:nedges;
    ivedges(edgedata.offedges) = [];
    
    for ii = ivedges
       if any( planehascoplanar(edgedata.planesatedge(1,:)) )
          
           ivcancelcombs = find( (edgeseespartialedge(:,ii) == -maxedgetoedgevisvalue) & ( (Bigthetae1(:,ii) == 0) | (Bigthetae1(:,ii) == 2*pi) ) );           
           if ~isempty(ivcancelcombs)
               edgeseespartialedge(ivcancelcombs,ii) = 0;
               edgeseespartialedge(ii,ivcancelcombs.') = 0;
           end
           
       end
        
    end
    
    % The only edge-to-edge combinations that we can easily discard are
    % those for which edges are aligned with each other (same zstart and
    % zend): only the closest of those should be kept!
    
    for ii = ivedges
       if any( planehascoplanar(edgedata.planesatedge(1,:)) )

           
           ivotheredges = find( (edgeseespartialedge(:,ii) == -maxedgetoedgevisvalue) & ((Bigthetae1(:,ii) == pi)) );
           if ~isempty(ivotheredges)
               
               nudgeval = 1e-10*edgelengthvec(ii)*100;
                ivselectedges = find( abs(Bigze1(ivotheredges,ii))<nudgeval | abs(Bigze2(ivotheredges,ii))<nudgeval);
                if length(ivselectedges) > 1
                   disp(['Edge no. ',int2str(ii)])
                    
%                     ivotheredges(ivselectedges)
                    lengthok = abs([Bigze1(ivotheredges(ivselectedges),ii)  Bigze2(ivotheredges(ivselectedges),ii)] - edgelengthvec(ii)) < nudgeval;
                    ivstillok = find(any(lengthok.'));
                    
                    if ~isempty(ivstillok)
                        ivotheredges = ivotheredges(ivselectedges(ivstillok));
                    else
                       ivotheredges = []; 
                    end

                    meanradialdist = mean( [Bigre1(ivotheredges,ii) Bigre2(ivotheredges,ii)].'    );
                    
                    ivshutoff = ivotheredges;
                    noshutoff = meanradialdist == min(meanradialdist);
                    ivshutoff(noshutoff) = [];

                    if ~isempty(ivshutoff)
                       edgeseespartialedge(ivshutoff,ii) = 0; 
                       edgeseespartialedge(ii,ivshutoff) = 0; 
                    end
                    
                end           
            
           end
           
           
       end
        
    end
    
    
	%-----------------------------------------------------------------------
	% Find identical combinations

% %     [B,I,J] = unique([Bigre1 Bigre2 Bigthetae1 Bigthetae2 Bigze1 Bigze2 edgelengthvec(:,ones(1,nedges))],'rows');
% %   
% %     pause
   
    if showtext >= 3
        disp('   Looking for identical combinations')    
    end
	[reftoshortlistE,re1sho,re2sho,thetae1sho,thetae2sho,...
	ze1sho,ze2sho,edgelengthsho,examplecombE] = EDcompress7(Bigre1,Bigre2,...
	Bigthetae1,Bigthetae2,Bigze1,Bigze2, edgedata.edgelengthvec(:,ones(1,nedges)).');    

else
	reftoshortlistE = [];
	re1sho = [];
	re2sho = [];
	thetae1sho = [];
	thetae2sho = [];
	ze1sho = [];
	ze2sho = [];
	examplecombE = [];
end

%----------------------------------------------------------------------------
%
%		SAVE THE VARIABLES
%
%----------------------------------------------------------------------------
%
% Also the variables from cadgeo??
% Yes: planeeqs planenvecs ncornersperplanevec planecorners corners
%      minvals maxvals

edgetoedgedata = struct('reftoshortlistE',reftoshortlistE,'re1sho',re1sho,...
    're2sho',re2sho,'thetae1sho',thetae1sho,'thetae2sho',thetae2sho,...
    'ze1sho',ze1sho,'ze2sho',ze2sho,'examplecombE',examplecombE,...
    'edgeseespartialedge',edgeseespartialedge,'edgealignedwithedge',edgealignedwithedge,...
    'edgeperptoplane',edgeperptoplane,'edgeplaneperptoplane1',edgeplaneperptoplane1,...
    'edgeplaneperptoplane2',edgeplaneperptoplane2);

if functionversion == 2
	elapsedtimeed2geo = etime(clock,t00);
    outpar = elapsedtimeed2geo;

	if filehandlingparameters.saveed2datafile == 1
    	eval(['save(''',desiredname,''',''planedata'',''edgedata'',''edgetoedgedata'',''EDinputdatahash'',''elapsedtimeed2geo'');'])
	end
end



