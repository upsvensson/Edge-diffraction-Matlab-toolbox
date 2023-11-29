function [outputstruct,elapsedtimeSRgeo,existingfilename] = EDSorRgeo...
    (planedata,edgedata,inputstruct,typeofcoords,EDversionnumber,filehandlingparameters)
% EDSorRgeo - Calculates some source- or receiver-related geometrical parameters.
% Calculates some source- or receiver-related geometrical parameters,
% based on corners,edges and planes in an eddata-file, and on a list of
% source/receiver coordinates. The output is returned as an expanded version
% of the input struct. For piston sources only one point of the piston is
% checked (the center point). The position is nudged a little bit in the
% direction of the plane normal vector.
%
% From v0.4 of the EDtoolbox, the input parameter list changed.
%
% Input parameters:
%   planedata               struct
%   edgedata                struct
%	inputstruct             struct (Sdata or Rdata) where the fields
%                           .coordinates, .nedgesubs, .sourcetype and 
%                           .pistonplanes (if Sdata) are the only ones used by this
%                           function. The function adds new fields to this existing struct.
%                           .nedgesubs: the number of subdivisions that each edge will be
%							subdivided into for visibility/obstruction tests. Default: 2
%							NB! When nedgesubs = 2, the two end points will be checked.
%   typeofcoords            'S' or 'R' - specifying if the point coordinates are sources
%                           or receivers. This determines what the output data in the output
%                           file will be called.
%   EDversionnumber
%   filehandlingparameters  a struct which contains the field showtext.
%
% Output parameters:
%   outputstruct            struct with fields
%       coordinates         (renamed) copy of the input parameter 'pointcoords'
%       visplanesfroms/visplanesfromr           Matrix, [nplanes,nsources]/[nplanes/nreceivers]
%                                           with integer values 0-5:
%                           0 means S/R behind a plane which is reflective or totabs
%                           1 means S/R is aligned with a plane, but outside it
%                           2 means S/R is in front of a plane which is reflective
%                           3 means S/R is in front of a plane which is totabs
%                           4 means S/R is inside a plane which is reflective
%                           5 means S/R is inside a plane which is totabs
%       vispartedgesfroms/vispartedgesfromr     Matrix, [nedges,nsources]/[nedges/nreceivers]
%                                           with integer values 0-2^nedgesubs-1 indicating
%                                           part-visibility of the edge.
%       soutsidemodel/routsidemodel             List, [nsources,1]/[nreceivers,1], with values
%                                           1 or 0 indicating whether the S/R is outside (1)
%                                           the model or not (0).
%                                           NB! Only simple tests are done, so an S/R could still
%                                           be outside the model even if the value here is 0.
%       reftoshortlistS/reftoshortlistR         Matrix, [nedges,nreceivers] or
%                                           [nedges,nsources] with pointers
%                                           to the three shortlists below.
%       rSsho/rRsho                             A shortlist with the unique
%                                           values of radius for sources/receivers relative to the the edges.
%       thetaSsho/thetaRsho
%       zSsho/zRsho
%   elapsedtimeSRgeo        This tells how long time was used inside this
%                           function. If an existing file was reused, then
%                           elapsedtimeSgeo has a second value which tells
%                           how much time was used for the existing file.
%   existingfilename        If an existing file was found that was reused,
%                           then the reused file name is given here. If no 
%                           existing file could be reused then this variable
%                           is empty. 
%
% NB! The text on the screeen, and in the code refers to 'R' or 'receivers'
% but it should be S or R.
%
% Uses the functions EDinfrontofplane, EDpoinpla, EDcompress3, EDcoordtrans1
% EDgetedgepoints, EDcheckobstr_pointtoedge, EDrecycleresultfiles
% from EDtoolbox.
% Uses the function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 29 Nov. 2023
%
% [outputstruct,elapsedtimeSRgeo,existingfilename] = EDSorRgeo(planedata,...
% edgedata,inputstruct,typeofcoords,EDversionnumber,filehandlingparameters);

%  1 June 2006 Functioning version
%  2 Dec 2014  Introduced new output parameters in the file: reftoshortlistR,rRsho,thetaRsho,zRsho
%  3 Dec 2014  Changed to two struct input parameters instead of inputfiles.
%  3 Dec 2014  Also changed to an output struct in addition to the file saving.
% 31 March 2015 Added '' in the file saving, to handle file and directory
%               names with blanks
% 1 Feb 2016 SM - Added 'rownumb' to planedata.planeeqs(rownumb,1:3) when EDpoinpla function is called
% 1 Feb 2016 SM - Changed 'edgeatplane' by 'edgedata.edgesatplane' and
%               'planehasindents' by ' planedata.planehasindents'
% 10 Oct. 2017 - Added the fields .vispartedgesfroms_start &
%               .vispartedgesfroms_end
% 27 Nov. 2017 Copied from ESIE2toolbox. Removed the input parameter
%               difforder, and removed the file saving.
% 17 Jan 2018 Turned off the check if an S/R is very close to a thin plane.
% 17 Jan 2018 Turned off some text printouts
% 8 Feb 2018 Introduced the EDinputdatahash
% 15 Mar 2018 Fixed a little bug for cases where some plane was TOTABS.
% 18 June 2018 Introduced the variable geomacc which is passed on to
% EDinfrontofplane and to EDpoinpla.
% 20 Jan 2021 Changed the format of .vispartedgesfroms_start to uint8 etc
% 20 Jan 2021 The obstruction test is skipped altogether when a convex
% model is used.
% 14 Mar 2021 Update to change of EDpoinpla
% 28 Sep 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 27 Oct. 2023 Renamed the field in the output struct from 'sources'/'receivers'
% to 'coordinates'. Also, changed so that the entire Sdata/Rdata struct is
% taken as input, so that it can be expanded.
% 29 Oct. 2023 nedgesubs was made into a field in the inputstruct
% 30 Oct. 2023 Fine-tuned the EDinputdatahash. Adjusts the piston center
% point a little bit for the visibility test; then nudges it back.
% 6 Nov. 2023 Added the piston coordinates and piston gauss order to the
% EDinputdatahash. And the source amplitudes.
% 10 Nov. 2023 Added a check which reduced the size of rrsho, thetarsho,
% zrsho. Previously, the edge-related coordinates were calculated for all
% S/R vs. all edges. Now a multiplication with (vispartedgesfromSR>0)
% reduces the number of values significantly, which makes the subesequent
% compression generate much shorter shortlists.
% 12 Nov. 2023 Introduced piston corner coordinates
% 20 Nov. 2023 Another reduction of the size of rrsho etc.
% 29 Nov. 2023 Adapted to the name change of the field pistongausspoints to
% pistongaussorder

t00 = clock;
geomacc = 1e-9;
nudgedist = 1e-5;   % for piston sources

showtext = filehandlingparameters.showtext;

nedgesubs = inputstruct.nedgesubs;

typeofcoords = lower(typeofcoords(1));
if typeofcoords~='r' && typeofcoords~='s'
    error('ERROR: The input parameter typeofcoords must have the value S or R')    
end

if typeofcoords == 's'
    EDinputdatastruct = struct('corners',planedata.corners,'planecorners',...
        planedata.planecorners,'offedges',edgedata.offedges,...
        'coordinates',inputstruct.coordinates,...
        'amplitudes',inputstruct.sourceamplitudes,...
        'pistoncoordinates',inputstruct.pistoncornercoordinates,...
        'pistongaussorder',inputstruct.pistongaussorder,...
        'nedgesubs',inputstruct.nedgesubs,'EDversionnumber',EDversionnumber);
else
    EDinputdatastruct = struct('corners',planedata.corners,'planecorners',...
        planedata.planecorners,'offedges',edgedata.offedges,...
        'coordinates',inputstruct.coordinates,...
        'nedgesubs',inputstruct.nedgesubs,'EDversionnumber',EDversionnumber);
end    
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Should the data be saved in a file? 

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else
    if typeofcoords == 's'
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_Sdata',EDinputdatahash);
    else
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_Rdata',EDinputdatahash);
    end    
end

if typeofcoords == 's'
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_Sdata.mat'];
else
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_Rdata.mat'];
end

if foundmatch == 1
    eval(['load(''',existingfilename,''')'])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    if typeofcoords == 's'
        outputstruct = Sdata;
        elapsedtimeSgeo_new = etime(clock,t00);
        elapsedtimeSRgeo = [elapsedtimeSgeo_new elapsedtimeSgeo];
    else
        outputstruct = Rdata;
        elapsedtimeRgeo_new = etime(clock,t00);
        elapsedtimeSRgeo = [elapsedtimeRgeo_new elapsedtimeRgeo];
    end
    return
end

%---------------------------------------------------------------
% No file was found, so this function should be run.

% ncorners = size(planedata.corners,1);
nplanes = size(planedata.planecorners,1);
maxncornersperplane = double(max(planedata.ncornersperplanevec));
nedges = size(edgedata.edgecorners,1);
nreceivers = size(inputstruct.coordinates,1);

n2 = size(planedata.canplaneobstruct,2);
% if n2>1
% %     canplaneobstruct = planedata.canplaneobstruct.';
% end

onesvec1R = ones(1,nreceivers,'uint8');
onesvec1ES = ones(1,nedgesubs);

% If a pistonsource was defined, then we need to nudge the center position
% away from the plane.

if typeofcoords == 's'
    if strcmp(inputstruct.sourcetype,'polygonpiston') == 1
        disp('Nudging!!')
        npistons = size(inputstruct.coordinates,1);
        for ii = 1:npistons
            nvec = planedata.planeeqs(inputstruct.pistonplanes(ii),1:3);
            inputstruct.coordinates(ii,:) = inputstruct.coordinates(ii,:) + ...
                nudgedist*nvec;
        end
    end
end

%##################################################################
%##################################################################
%##################################################################
%
%		PLANE RELATED PARAMETERS
%
%##################################################################

% activeplanelist = find(planedata.reflfactors ~= 0);
totabsplanelist = find(planedata.reflfactors == 0);
ntotabsplanes = length(totabsplanelist);
nthinplanes = length(find(planedata.planeisthin));

if showtext >= 3
    if lower(typeofcoords(1)) == 'r'
    	disp('         Checking visible planes from R')
    else
    	disp('         Checking visible planes from S')        
    end
end

%--------------------------------------------------------------
%
%		visplanesfromr      [nplanes,nrec]  uint8
%
%           0 means S/R behind a plane which is reflective or totabs
%           1 means S/R is aligned with a plane, but not inside the plane
%           2 means S/R is in front of a plane which is reflective
%           3 means S/R is in front of a plane which is totabs
%           4 means S/R is inside a plane which is reflective
%           5 means S/R is inside a plane which is totabs
%
% We can call EDinfrontofplane with a single call by extracting receiver
% numbers and plane numbers for the [nplanes,nsources] matrix.
%
% EDinfrontofplane returns -1 for behind, 0 for in-plane with and 1 for in
% front of so if we add the value 1 we get close to the final value we want
% to have. Another modification is to give totabs planes the value 3
% instead of 2. A last modification is to find which S/R are inside the
% plane.

iv = [1:nplanes*nreceivers].';                
if nreceivers < 256
    colnumb = uint8(ceil(iv/nplanes));             % This is the receiver number
elseif nreceivers < 65536
    colnumb = uint16(ceil(iv/nplanes));             % This is the receiver number
else
    colnumb = uint32(ceil(iv/nplanes));             % This is the receiver number    
end
if nplanes < 256
    rownumb = uint8(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
elseif nplanes < 65536
    rownumb = uint16(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
else
    rownumb = uint32(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
end
clear iv 

% The function EDinfrontofplane returns:
%   +1  if point is in front of plane
%    0  if point belongs to (infinite) plane
%   -1 if point is behind plane
% The addition of 1 leads to that visplanesfromr is:
%   +2  if point is in front of plane
%   +1  if point belongs to (infinite) plane
%    0 if point is behind plane

%%%visplanesfromr = EDinfrontofplane(pointcoords(colnumb,:),planedata.planeeqs(rownumb,1:3),...
visplanesfromr = EDinfrontofplane(inputstruct.coordinates(colnumb,:),planedata.planeeqs(rownumb,1:3),...
planedata.corners(planedata.planecorners(rownumb,1),:),planedata.corners(planedata.planecorners(rownumb,2),:),'','',geomacc) + 1;

% If any source/receiver belongs to a plane, we should check if the S/R is
% inside the finite plane or not:
%   Yes -> issue an error message, because source/receivers can not be too
%          close to any plane
%   No  -> no error message, but if difforder > 1, issue a warning message
%          that ESIEBEM is recommended.

iv_closetoplane = find(visplanesfromr==1);
if any(iv_closetoplane)
    [hitvec,edgehit,edgehitnumbers,cornerhit,cornerhitnumbers] = EDpoinpla(inputstruct.coordinates(colnumb(iv_closetoplane),:),rownumb(iv_closetoplane),...
        planedata.minvals,planedata.maxvals,planedata.planecorners,planedata.corners,...
        planedata.ncornersperplanevec,planedata.planeeqs(:,1:3),geomacc);
    ivinside = find(hitvec);
    ivedgehit = find(edgehit);
    ivcornerhit = find(cornerhit);
    ivoutside = find(hitvec==0);

    if any(ivinside)
        error(['ERROR: at least ',upper(typeofcoords),' no. ',int2str( colnumb(iv_closetoplane(ivinside(1)))),...
            ' is too close to plane no. ',int2str(  rownumb(iv_closetoplane(ivinside(1)))  )])
    end    
    if any(edgehit)
        error(['ERROR: at least ',upper(typeofcoords),' no. ',int2str( colnumb(iv_closetoplane(ivedgehit(1)))),...
            ' is too close to one edge'])
    end    
    if any(cornerhit)
        error(['ERROR: at least ',upper(typeofcoords),' no. ',int2str( colnumb(iv_closetoplane(ivedgehit(1)))),...
            ' is too close to one corner'])
    end    
    if any(ivoutside)
        error(['ERROR: at least ',upper(typeofcoords),' no. ',int2str( colnumb(iv_closetoplane(ivoutside(1)))),...
            ' is too close to the extension of plane no. ',int2str(  rownumb(iv_closetoplane(ivoutside(1)))  )])
    end
end

clear rownumb colnumb

if ntotabsplanes > 0

    % It's unclear why the first uint8 assignment doesn't work!
    % Bugfix 15 Mar 2018: the line "visplanesfromr(iv) = ..." needed a
    % "double" for both terms.
        
    colvec = (0:nreceivers-1)*nplanes;
    iv = uint32(totabsplanelist(:,onesvec1R) + colvec(ones(ntotabsplanes,1),:));    
    visplanesfromr(iv) = uint8(double(visplanesfromr(iv)) + double(visplanesfromr(iv)==2));
    clear iv colvec
    visplanesfromr = uint8(visplanesfromr);
end

% For all the planes that the S/R are aligned with,
% check if the S/R is inside or outside the plane.
%
% We can call EDpoinpla with a single call by extracting S/R coordinates
% (=colnumb) and planenumbers (=rownumb).

ivec_visR1 = find(visplanesfromr==1);
if ~isempty(ivec_visR1)
	if nreceivers < 256
        colnumb = uint8(ceil(ivec_visR1/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec_visR1/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec_visR1/nplanes));             % This is the receiver number    
	end
%	rownumb = ivec_visR1 - (double(colnumb)-1)*nplanes;       % This is the plane number
	
    if nplanes < 256
        rownumb = uint8(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
    elseif nplanes < 65536
        rownumb = uint16(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
    else
        rownumb = uint32(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
    end
    [hitvec,edgehit,edgehitnumbers,cornerhit,cornerhitnumbers] = ...
        find(EDpoinpla(inputstruct.coordinates(colnumb,:),rownumb,...
        planedata.minvals,planedata.maxvals,planedata.planecorners,...
        planedata.corners,planedata.ncornersperplanevec,planedata.planeeqs(rownumb,1:3)));
    if ~isempty(hitvec)
        visplanesfromr(ivec_visR1(hitvec)) = uint8(ones(size(hitvec))*4);
        if ntotabsplanes > 0
            insidetotabsplane = (reflfactors(rownumb(hitvec))==0);
            ivec2 = find(insidetotabsplane);
            if ~isempty(ivec2)
                visplanesfromr(ivec_visR1(hitvec(ivec2))) = uint8(ones(size(hitvec(ivec2)))*5);
            end
        end
        ivec_visR1(hitvec) = [];
    end
end

visplanesfromr = reshape(visplanesfromr,nplanes,nreceivers);

%--------------------------------------------------------------
%
%       routsidemodel       (0 or 1, size [1,nrec])
%
% We make a simple check: if the S/R is outside the big cube defined by all
% the planes, and we have an interior problem, then S/R must be outside the
% model.

routsidemodel = uint8(zeros(1,nreceivers));
% % % % if int_or_ext_model == 'i'
% % % %     ivec = find(pointcoords(:,1)<totalmodelmin(1) | pointcoords(:,2)<totalmodelmin(2) | pointcoords(:,3)<totalmodelmin(3) | ...
% % % %           pointcoords(:,1)>totalmodelmax(1) | pointcoords(:,2)>totalmodelmax(2) | pointcoords(:,3)>totalmodelmax(3));
% % % % 	if ~isempty(ivec)
% % % %         routsidemodel(ivec) = ones(size(ivec));
% % % % 	end
% % % % 	routsidemodel = uint8(routsidemodel);
% % % %     
% % % %     % If we have a convex model, then the obstruction check is turned off,
% % % %     % so we must check here if the source/receiver is inside.
% % % %     if sum(closwedangvec<pi) == 0
% % % % 
% % % %         nplanesvisible = sum(sign(visplanesfromr));
% % % %         ivec = find(nplanesvisible < nplanes);
% % % %         if ~isempty(ivec)
% % % %             routsidemodel(ivec) = ones(size(ivec));        
% % % %         end
% % % %     end
% % % % end

%##################################################################
%##################################################################
%##################################################################
%
%		EDGE RELATED PARAMETERS
%
%##################################################################

% % % histvec = [1:nedges];
% closwedlargerthanpi = edgedata.closwedangvec>pi;
% closwedsmallerthanpi = edgedata.closwedangvec<pi;

if showtext >= 3
    if lower(typeofcoords(1)) == 'r'
    	disp('         Checking which edges are seen from R')
    else
    	disp('         Checking which edges are seen from S')        
    end    
end

%--------------------------------------------------------------
%
%		visedgesfromr      [nedges,nrec]   uint8
%       
%       6    Edge belongs to a plane which is aligned with R and thin and
%            rigid (but R not inside the plane)
%       5    Edge belongs to a plane which is aligned with R and not thin 
%
% These can be derived from the cases where visplanesfromr = 1

visedgesfromr = uint8(ones(nedges,nreceivers));

if ~isempty(ivec_visR1)
	if nreceivers < 256
        colnumb = uint8(ceil(ivec_visR1/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec_visR1/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec_visR1/nplanes));             % This is the receiver number    
	end
	if nplanes < 256
        rownumb = uint8(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	end

    % Divide these lists into two categories: thin planes and non-thin
    % planes.
    
    iv2 = find(planedata.planeisthin(rownumb));
    if ~isempty(iv2)
        colnumb_thin = colnumb(iv2);            
        rownumb_thin = rownumb(iv2);
        colnumb(iv2) = [];
        rownumb(iv2) = [];
    else
        colnumb_thin = [];
        rownumb_thin = [];
    end

    if ~isempty(colnumb)
        % Select all the edges that are connected to these planes
        
        if nedges < 256
            edgeselection = uint8(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        else
            edgeselection = uint32(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        end
        if length(rownumb) > 1
            maxcols = sum(sign(sum(double(edgeselection))));
        else
            maxcols = sum(sign(double(edgeselection)));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb = colnumb(:,ones(1,maxcols));
        
        indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
        indexvec(edgeselection==0) = [];
        
        visedgesfromr(indexvec) = uint8(5*ones(size(indexvec)));
    end
    if ~isempty(colnumb_thin)
        % Select all the edges that are connected to these planes
	
%        edgeselection = edgedata.edgesatplane(rownumb_thin,1:maxncornersperplane);
        if nedges < 256
            edgeselection = uint8(edgedata.edgesatplane(rownumb_thin,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgedata.edgesatplane(rownumb_thin,1:maxncornersperplane));
        else
            edgeselection = uint32(edgedata.edgesatplane(rownumb_thin,1:maxncornersperplane));
        end
        if length(rownumb_thin) > 1
            maxcols = sum(sign(sum(edgeselection)));
        else
            maxcols = sum(sign(edgeselection));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb_thin = colnumb_thin(:,ones(1,maxcols));
 
        % The line below was a bug: 030731
%         indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
         indexvec = uint32(double(edgeselection) + (double(colnumb_thin)-1)*nedges);
          indexvec(edgeselection==0) = [];
        
        visedgesfromr(indexvec) = uint8(6*ones(size(indexvec)));
    end
end

%       4    Edge belongs to a plane which the S/R is inside,
%            and the plane is totabs
%       3    Edge belongs to a plane which the S/R is inside,
%            and the plane has indents
%       2    Edge belongs to a plane which the S/R is inside,
%            and the plane has no indents
%
% These can be derived from the cases where visplanesfromr = 4 and 5

ivec = find(visplanesfromr==5);
if ~isempty(ivec)
%     colnumb = ceil(ivec/nplanes);               % This is the receiver number
	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
%     rownumb = ivec - (double(colnumb)-1)*nplanes;       % This is the plane number
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Select all the edges that are connected to these planes
	
%    edgeselection = edgedata.edgesatplane(rownumb,1:maxncornersperplane);
    if nedges < 256
        edgeselection = uint8(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
    elseif nedges < 65536
        edgeselection = uint16(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
    else
        edgeselection = uint32(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
    end
    if length(rownumb) > 1
        maxcols = sum(sign(sum(edgeselection)));
    else
        maxcols = sum(sign(edgeselection));    
    end
    edgeselection = edgeselection(:,1:maxcols);
    colnumb = colnumb(:,ones(1,maxcols));
    
    indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
    indexvec(edgeselection==0) = [];
    
    visedgesfromr(indexvec) = uint8(4*ones(size(indexvec)));
end

ivec = find(visplanesfromr==4);    
if ~isempty(ivec)
%    colnumb = ceil(ivec/nplanes);               % This is the receiver number
	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
%    rownumb = ivec - (double(colnumb)-1)*nplanes;       % This is the plane number
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Divide these lists into two categories: planes with indents and 
    % planes without
    
    iv2 = find(planedata.planehasindents(rownumb));
    if ~isempty(iv2)
        colnumb_ind = colnumb(iv2);            
%         rownumb_ind = rownumb(iv2);
        colnumb(iv2) = [];
        rownumb(iv2) = [];
    else
        colnumb_ind = [];
%         rownumb_ind = [];
    end
    if ~isempty(colnumb)
        % Select all the edges that are connected to these planes
	
    %    edgeselection = edgedata.edgesatplane(rownumb,1:maxncornersperplane);
        if nedges < 256
            edgeselection = uint8(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        else
            edgeselection = uint32(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
        end
        if length(rownumb) > 1
            maxcols = sum(sign(sum(double(edgeselection))));
        else
            maxcols = sum(sign(double(edgeselection)));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb = colnumb(:,ones(1,maxcols));
        
        indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
        indexvec(edgeselection==0) = [];
        
        visedgesfromr(indexvec) = uint8(2*ones(size(indexvec)));
    end
    if ~isempty(colnumb_ind)
        % Select all the edges that are connected to these planes
	
        edgeselection = edgedata.edgesatplane(rownumb,1:maxncornersperplane);
        if length(rownumb) > 1
            maxcols = sum(sign(sum(edgeselection)));
        else
            maxcols = sum(sign(edgeselection));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb_ind = colnumb_ind(:,ones(1,maxcols));
        
        indexvec = uint32(edgeselection + (colnumb_ind-1)*nedges);
        indexvec(edgeselection==0) = [];
        
        visedgesfromr(indexvec) = uint8(3*ones(size(indexvec)));
    end
end

%       0   R can never see the edge because it is behind the edge-planes.
%           NB! If the closwedang > pi, then it is enough if R is behind
%           one plane. If the closwedang <= pi, then R must be behind both
%           planes in order not to see the edge.

ivec = find(visplanesfromr==0);

if ~isempty(ivec)
%    colnumb = ceil(ivec/nplanes);               % This is the receiver number
 	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
%    rownumb = ivec - (double(colnumb)-1)*nplanes;       % This is the plane number
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Select all the edges that are connected to these planes

    if nedges < 256
        edgeselection = uint8(edgedata.edgesatplane(rownumb,1:maxncornersperplane));        
    elseif nedges < 65536
        edgeselection = uint16(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
    else
        edgeselection = uint32(edgedata.edgesatplane(rownumb,1:maxncornersperplane));
    end
    if length(rownumb) > 1
        maxcols = sum(sign(sum(double(edgeselection))));
    else
        maxcols = sum(sign(double(edgeselection)));    
    end
    edgeselection = edgeselection(:,1:maxcols);
    colnumb = colnumb(:,ones(1,maxcols));
    ntot2 = length(rownumb)*maxcols;
    
    indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
    indexvec = reshape(indexvec,ntot2,1);
    edgeselection = reshape(edgeselection,ntot2,1);

    iv2 = find(edgeselection==0);    
    indexvec(iv2) = [];
    edgeselection(iv2) = [];
    ntot2 = ntot2 - length(iv2);
    
    [indexvec,sortvec] = sort(indexvec);
    edgeselection = edgeselection(sortvec);
    clear sortvec
    
    markdoubles = (edgeselection(1:ntot2-1)==edgeselection(2:ntot2));
    markdoubles = sign([markdoubles;0] + [0;markdoubles]);
    iv2 = uint32(find(markdoubles));
    clear markdoubles
    visedgesfromr(indexvec(iv2)) = 0;
    indexvec(iv2) = [];
    edgeselection(iv2) = [];
    iv2 = uint32(find(edgedata.closwedangvec(edgeselection)<=pi));
    indexvec(iv2) = [];
    visedgesfromr(indexvec) = 0;
    clear indexvec
end

clear edgeselection iv2

%##################################################################
%
% Mask out edges that should be switched off

mask = (ones(nedges,1));
mask(edgedata.offedges) = mask(edgedata.offedges) & 0;

% % % % visedgesfroms = uint8(double(visedgesfroms).*mask(:,onesvec1S));

visedgesfromr = uint8(double(visedgesfromr).*mask(:,onesvec1R));

clear onesvec1R

%##################################################################
%##################################################################
%##################################################################
%
%		CHECK OBSTRUCTION OF R-TO-EDGE PATHS
% 
% This check is not needed if planedata.modeltype = 'convex_ext'
%
%   vispartedgesfromr   0-2^nedgesubs-1, [nedges,nrec]
%                       telling which segments of an edge that
%                       are visible
%
% A. If visedgesfromr = 6 then check aligned-plane obstructions
%           If OK,     set visedgesfromr = 1 (check ordinary obst.)
%           If not OK, set visedgesfromr = 0
%
% B. If visedgesfromr = 3 then check within-plane obstructions
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = 0-full
%
% C. If visedgesfromr = 2 then
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = Full
%
% D. If visedgesfromr = 1 then check ordinary obstructions
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = 0-full

% edgedividervec = [0:nedgesubs-1].';
% weightvec = 2.^edgedividervec;
maxvisibilityval = 2^nedgesubs-1;

if nedgesubs<8
    vispartedgesfromr = zeros(size(visedgesfromr),'uint8');
elseif nedgesubs<16
    vispartedgesfromr = zeros(size(visedgesfromr),'uint16');
else
    vispartedgesfromr = zeros(size(visedgesfromr),'uint32');
end

if strcmp(planedata.modeltype,'convex_ext') == 1
    iv = uint32(find(visedgesfromr==1));
    vispartedgesfromr(iv) = maxvisibilityval;
else
    
    if showtext >= 3
        if lower(typeofcoords(1)) == 'r'
            disp('         Checking for obstructing planes between R and edges')
        else
            disp('         Checking for obstructing planes between S and edges')
        end
    end

    iv = uint32(find(visedgesfromr>=5));
    if ~isempty(iv)
        if showtext > 1
            disp('      We check aligned-edge obstructions')
        end

        % All edges must be checked vs all other edges?
        % Make subdivision of edges. If a line from R to edge segments
        % pass through a plane that is constructed by an edge and perpendicular to
        % the studied plane, then that edge is obstructing. 

        combnumbers = double(iv);
        recnumbers = ceil(combnumbers/nedges);
    %     if nedges < 256
    % % 		edgenumbers = uint8(combnumbers - (recnumbers-1)*nedges);
    %     elseif nedges < 65536
    % % 		edgenumbers = uint16(combnumbers - (recnumbers-1)*nedges);
    %     else
    % %    		edgenumbers = uint32(combnumbers - (recnumbers-1)*nedges);
    %     end
        nedgesperreceiver = histc(recnumbers,(1:nreceivers));
        iv1 = find(nedgesperreceiver > 2, 1);
        if isempty(iv1)
            if showtext > 2
                disp('   No aligned-edges can obscure each other')    
            end
            visedgesfromr(iv)    = ones(size(iv));
        else
            error('ERROR: Obstruction check of aligned-edges not implemented yet!')        
        end
     end

    iv = uint32(find(visedgesfromr==3));
    if ~isempty(iv)
        if showtext > 2
            disp('      We check within-plane obstructions')
        end
    %     visedgesfromr(iv)      = zeros(size(iv));
        error('ERROR: Obstruction check of within-same-plane-edges not implemented yet!')        
        %    vispartedgesfromr(iv) = 0 - full;
    end

    iv = uint32(find(visedgesfromr==2));
    if ~isempty(iv)
        if showtext > 2
            disp('      Edge fully visible')
        end
        visedgesfromr(iv)      = zeros(size(iv));
        vispartedgesfromr(iv)  =  maxvisibilityval*ones(size(iv));

    end

    iv = uint32(find(visedgesfromr==4));
    if ~isempty(iv)
        if showtext > 2
            disp('      Edge in totabs plane')
        end
        visedgesfromr(iv)      = zeros(size(iv));
        vispartedgesfromr(iv)  =  zeros(size(iv));

    end

    % Below is the main part of the obstruction test. In the matrix
    % visedgesfromr, values 1 indicate that edge-receiver combos should be
    % tested. visedgesfromr has the size [nedges,nreceivers].

    iv = uint32(find(visedgesfromr==1));
    if ~isempty(iv)

    %     visedgesfromr(iv) = zeros(size(iv));
        combnumbers = double(iv);
        clear iv

        % combnumbers is a list of the "combined" values for edge and receiver. 

        if nreceivers < 256
            recnumbers = uint8(ceil(combnumbers/nedges));
        elseif nreceivers < 65536
            recnumbers = uint16(ceil(combnumbers/nedges));
        else
            recnumbers = uint32(ceil(combnumbers/nedges));
        end
        if nedges < 256
            edgenumbers = uint8(combnumbers - (double(recnumbers)-1)*nedges);
        elseif nedges < 65536
            edgenumbers = uint16(combnumbers - (double(recnumbers)-1)*nedges);            
        else
            edgenumbers = uint32(combnumbers - (double(recnumbers)-1)*nedges);            
        end
        combnumbers = uint32(combnumbers);
        ncombs = length(recnumbers);
    %     ntot = ncombs*nplanes;

        % Mark, in a big matrix, which planes can at all obstruct.
        % Planes that can obstruct must be seen by the receiver or the edge but
        % not both - because then they would be on the same side!
        % NB! Also totabs planes can obstruct.
        % The big matrix will have nedges*nplanes rows. The vector planenumb
        % will have values such as [1 1 1 1... 2 2 2 2 .....]

        iv = [1:ncombs*nplanes].';                
        % Plane number is given by the col no.
        if nplanes < 256
            planenumb = uint8(ceil(iv/ncombs));
        elseif nplanes < 65536
            planenumb = uint16(ceil(iv/ncombs));
        else
            planenumb = uint32(ceil(iv/ncombs));        
        end
        indexvec = planenumb;

        % The rownumber will simply have values [1 2 3... N 1 2 3 ... N 1 2 3
        % ... N] where N is the number of combinations.
        if ncombs < 256
            rownumber = uint8(iv - (double(planenumb)-1)*ncombs);
        elseif ncombs < 65536
            rownumber = uint16(iv - (double(planenumb)-1)*ncombs);
        else
            rownumber = uint32(iv - (double(planenumb)-1)*ncombs);
        end
        clear planenumb iv 

        % The peak in the memory need might be after the two indexvec lines
        % below

        % indexvec2 will point to the right location in an [nplanes,nedges] matrix
        % indexvec will point to the right location in an [nplanes,nreceivers] matrix
        indexvec2 = uint32(double(indexvec) + (double(edgenumbers(rownumber))-1)*nplanes);
        indexvec = uint32(double(indexvec) + (double(recnumbers(rownumber))-1)*nplanes);
        clear rownumber

        % The big matrix checkplane, size [1,ncombs*nplanes], will contain the
        % value 1 if a plane should be checked for obstruction.
        % The index order of checkplane is:
        % [Plane1&comb1 Plane1&comb2 Plane1&comb3 ...] where comb1 is the first
        % combination of receiver and edge in the recnumbers & edgenumbers
        % lists.
        % Only planes that are seen by the S/R *or* the edge, but not both
        % should be checked for obstruction!!
        % We remove combinations where:
        %   ... S/R is aligned with a plane (visplanesfromr == 1)
        %   ... edge belongs to a plane or is aligned with a plane (edgeseesplane <0)
        %   ... S/R is behind and edge is behind a plane (visplanesfromr == 0 & edgeseesplane == 0)
        %   ... S/R is in front of and edge is in front of a plane ( (visplanesfromr == 2 | visplanesfromr == 3) & edgeseesplane == 1)
        %
        % Comment 050116 (PS):
        %   We would actually need a matrix called planeseesedge. Here we use
        %   edgeseesplane instead, and it is not completely clear if that works
        %   fine. 

        checkplane = (visplanesfromr(indexvec)~=1) & (edgedata.edgeseesplane(indexvec2)>=0) & not(visplanesfromr(indexvec)==0 & edgedata.edgeseesplane(indexvec2)==0 ) & not( (visplanesfromr(indexvec)==2 | visplanesfromr(indexvec)==3) & edgedata.edgeseesplane(indexvec2)==1 );
        if size(checkplane,1) > size(checkplane,2)
            checkplane = checkplane.';    
        end

        if size(checkplane,1) == 1 || size(checkplane,2) == 1
            checkplane = reshape(checkplane,ncombs,nplanes);
        end

        clear indexvec indexvec2 edgeseesplane

        % If there are some R-edge combos that have no planes to check
        % obstruction for, we can mark those combos ("combsnottocheck") as fully visible.
        %
        % The remaining ones ("combstocheck") still need to be checked.

        n1 = size(checkplane,1);
        if n1 < 65536
            combstocheck = uint16(find( sum(checkplane.') ));        
        elseif n1 < 4e9
            combstocheck = uint32(find( sum(checkplane.') ));                    
        end
        combsnottocheck = (1:ncombs);
        combsnottocheck(combstocheck) = [];
        if ~isempty(combsnottocheck)
            vispartedgesfromr(combnumbers(combsnottocheck)) = maxvisibilityval*ones(size(combsnottocheck));
        end

        if ~isempty(combstocheck)

            checkplane = checkplane(combstocheck,:);
            recnumbers = recnumbers(combstocheck);
            maxrec = max(recnumbers);
            if maxrec < 256
                recnumbers = uint8(recnumbers);        
            elseif maxrec < 65536
                recnumbers = uint16(recnumbers);                
            else
                recnumbers = uint32(recnumbers);                
            end
            edgenumbers = edgenumbers(combstocheck);
    %         combnumbers = combnumbers(combstocheck);
            ncombs = length(combstocheck);

            % Now, checkplane is a matrix of size [ncombs,nplanes] where each
            % row corresponds to one path (from one receiver to one edge) that needs
            % an obstruction check. For that row, checkplane has the value 1 for
            % each plane that needs to be checked.
            %    
            % Expand all edges into their edge segment/subdivision
            % We treat all the little edge subdivisions as separate edges

            nposs = ncombs*nedgesubs;        

            expandedrecnumbers = recnumbers(:,onesvec1ES);
            clear recnumbers
            expandedrecnumbers = reshape(expandedrecnumbers.',nposs,1);
    %         expandedcombnumbers = combnumbers(:,onesvec1ES);
            clear combnumbers
    %         expandedcombnumbers = reshape(expandedcombnumbers.',nposs,1);
            if nposs < 65536
    %             okcombs = zeros(size(expandedrecnumbers),'uint16');
            elseif nposs < 4e9
    %             okcombs = zeros(size(expandedrecnumbers),'uint32');                
            end

             [tocoords,expandededgeweightlist,expandededgenumbers] = EDgetedgepoints(edgedata.edgestartcoords(edgenumbers,:),...
                 edgedata.edgeendcoords(edgenumbers,:),nedgesubs,1);
             expandededgenumbers = edgenumbers(expandededgenumbers);

%%%            [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDcheckobstr_pointtoedge(pointcoords,expandedrecnumbers,tocoords,reshape(repmat(checkplane.',[nedgesubs,1]),nplanes,nposs).',planedata);
            [nonobstructedpaths,nobstructions,edgehits,cornerhits] = ...
                EDcheckobstr_pointtoedge(inputstruct.coordinates,expandedrecnumbers,...
                tocoords,reshape(repmat(checkplane.',[nedgesubs,1]),nplanes,nposs).',planedata);

    %         expandedcombnumbers = expandedcombnumbers(nonobstructedpaths);
            expandededgeweightlist = expandededgeweightlist(nonobstructedpaths);
            expandedrecnumbers = expandedrecnumbers(nonobstructedpaths);
            expandededgenumbers = expandededgenumbers(nonobstructedpaths);

            % Pack all non-obstructed edge segments together and add their weights together

            test = [double(expandedrecnumbers) double(expandededgenumbers)];
            ncombs = length(expandedrecnumbers);
            dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
            ivremove = find(dtest==1);

            while ~isempty(ivremove)
                expandededgeweightlist(ivremove+1) = double(expandededgeweightlist(ivremove+1)) + double(expandededgeweightlist(ivremove));
                expandededgeweightlist(ivremove) = [];
                expandedrecnumbers(ivremove) = [];
                expandededgenumbers(ivremove,:) = [];

                test = [double(expandedrecnumbers) double(expandededgenumbers)];
                ncombs = length(expandedrecnumbers);
                dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                ivremove = find(dtest==1);
            end

            indexvec = uint32(nedges*(double(expandedrecnumbers)-1)+double(expandededgenumbers));
            vispartedgesfromr(indexvec)  =  expandededgeweightlist;
            clear indexvec                
        end
    end
end

%##################################################################
%
% Compute the edge-related receiver coordinates

% if nreceivers*nedges < 256
%     reftoshortlistR = zeros(nedges,nreceivers,'uint8');
% elseif nreceivers*nedges < 65536
%     reftoshortlistR = zeros(nedges,nreceivers,'uint16');
% else
%     reftoshortlistR = zeros(nedges,nreceivers,'uint32');    
% end

subcoordinates = 0;
if typeofcoords == 's'
    if strcmp(inputstruct.sourcetype,'polygonpiston') == 1
        subcoordinates = 1;
    end
end

if subcoordinates == 0
    rRcomplete     = zeros(nedges,nreceivers);
    thetaRcomplete = zeros(nedges,nreceivers);
    zRcomplete     = zeros(nedges,nreceivers);
    
    for ii = 1:nedges
        % New 20 Nov. 2023: if one edge doesn't see a single S/R, then
        % don't compute the edge-related coordinates for that edge
        if any(vispartedgesfromr(ii,:))
            [rR,thetaR,zR] = EDcoordtrans1(inputstruct.coordinates,...
                [edgedata.edgestartcoords(ii,:);edgedata.edgeendcoords(ii,:)],...
                edgedata.edgenvecs(ii,:),reshape(edgedata.edgerelatedcoordsysmatrices(ii,:),3,3));
            rRcomplete(ii,:) = rR.';
            thetaRcomplete(ii,:) = thetaR.';
            zRcomplete(ii,:) = zR.';    
        end
    end
    
    % New 10 Nov. 2023:
    rRcomplete = rRcomplete.*(vispartedgesfromr>0);
    thetaRcomplete = thetaRcomplete.*(vispartedgesfromr>0);
    zRcomplete = zRcomplete.*(vispartedgesfromr>0);
    
    [reftoshortlistR,rRsho,thetaRsho,zRsho,~] = ...
    EDcompress3(rRcomplete,thetaRcomplete,zRcomplete);
else
    nsubcoordinates = size(inputstruct.pistongausscoordinates{1},1)
    rRcomplete     = zeros(nedges,nreceivers,nsubcoordinates);
    thetaRcomplete = zeros(nedges,nreceivers,nsubcoordinates);
    zRcomplete     = zeros(nedges,nreceivers,nsubcoordinates);
    for ii = 1:nedges
        for jj = 1:nreceivers
            if vispartedgesfromr(ii,jj) > 0
                [rR,thetaR,zR] = EDcoordtrans1(inputstruct.pistongausscoordinates{jj},...
                    [edgedata.edgestartcoords(ii,:);edgedata.edgeendcoords(ii,:)],...
                    edgedata.edgenvecs(ii,:),reshape(edgedata.edgerelatedcoordsysmatrices(ii,:),3,3));
                rRcomplete(ii,jj,:) = rR.';
                thetaRcomplete(ii,jj,:) = thetaR.';
                zRcomplete(ii,jj,:) = zR.';  
            end
        end
    end
    
    [reftoshortlistR,rRsho,thetaRsho,zRsho,~] = ...
    EDcompress3(rRcomplete,thetaRcomplete,zRcomplete);

end
%----------------------------------------------------------------------------
%
%		STORE THE VARIABLES IN THE STRUCT
%
%--------------------------------------------------------------------------

if nedgesubs<8
    vispartedgesfromr_start = zeros(size(vispartedgesfromr),'uint8');
    vispartedgesfromr_end   = ones(size(vispartedgesfromr),'uint8');
elseif nedgesubs<16
    vispartedgesfromr_start = zeros(size(vispartedgesfromr),'uint16');
    vispartedgesfromr_end   = ones(size(vispartedgesfromr),'uint16');
else
    vispartedgesfromr_start = zeros(size(vispartedgesfromr),'uint32');
    vispartedgesfromr_end   = ones(size(vispartedgesfromr),'uint32');
end

if typeofcoords=='r'

    outputstruct = inputstruct;
    outputstruct.visplanesfromr = visplanesfromr;
    outputstruct.vispartedgesfromr = vispartedgesfromr;
    outputstruct.routsidemodel = routsidemodel;
    outputstruct.vispartedgesfromr_start = vispartedgesfromr_start;
    outputstruct.vispartedgesfromr_end = vispartedgesfromr_end;
    outputstruct.reftoshortlistR = reftoshortlistR;
    outputstruct.rRsho = rRsho;
    outputstruct.thetaRsho = thetaRsho;
    outputstruct.zRsho = zRsho;

    elapsedtimeRgeo = etime(clock,t00);
    elapsedtimeSRgeo = elapsedtimeRgeo;

    if filehandlingparameters.saveSRdatafiles == 1
        Rdata = outputstruct;
    	eval(['save(''',desiredname,''',''Rdata'',''EDinputdatahash'',''elapsedtimeRgeo'');'])
    end   
    
else

    if strcmp(inputstruct.sourcetype,'polygonpiston') == 1
        disp('Nudging back!!')
        npistons = size(inputstruct.coordinates,1);
        for ii = 1:npistons
            nvec = planedata.planeeqs(inputstruct.pistonplanes(ii),1:3);
            inputstruct.coordinates(ii,:) = inputstruct.coordinates(ii,:) - ...
                nudgedist*nvec;
        end

        % Create rotation matrices, one for each pistonedge
        nmaxcornernumbersperpiston = size(inputstruct.pistoncornernumbers,2);
        nrotationmatrices = numel(inputstruct.pistoncornernumbers);
        rotationmatrices = cell(nrotationmatrices,1);
        reftorotationmatrices = [1:nrotationmatrices].';
        reftorotationmatrices = reshape(reftorotationmatrices,nmaxcornernumbersperpiston,npistons).';
        for ii = 1:npistons
            nvec = planedata.planeeqs(inputstruct.pistonplanes(ii),1:3);
            for jj = 1:nmaxcornernumbersperpiston
                startnumber = inputstruct.pistoncornernumbers(ii,jj);
                if startnumber > 0
                    if jj == nmaxcornernumbersperpiston
                        endnumber = inputstruct.pistoncornernumbers(ii,1);
                    else
                        endnumber = inputstruct.pistoncornernumbers(ii,jj+1);
                        if endnumber == 0
                            endnumber = inputstruct.pistoncornernumbers(ii,1);
                        end
                    end
%                     disp(['ii = ',int2str(ii),', jj = ',int2str(jj),' start no. = ',int2str(startnumber),' end no. = ',int2str(endnumber)])
%                     disp(['   Rot. matrix no. ',int2str(reftorotationmatrices(ii,jj))])
                    pistonedgevec = inputstruct.pistoncornercoordinates(endnumber,:) - ...
                        inputstruct.pistoncornercoordinates(startnumber,:);
                    pistonedgevec = pistonedgevec/norm(pistonedgevec);
                    yvec = cross(nvec,pistonedgevec);
                    singlerotmatrix = [0 1 0;0 0 1;1 0 0]*inv([nvec(:) pistonedgevec(:) yvec(:)]);
                    rotationmatrices{reftorotationmatrices(ii,jj)} = singlerotmatrix;
                else
                end
            end
        end
        inputstruct.pistonrotationmatrices = rotationmatrices;
        inputstruct.reftorotationmatrices = reftorotationmatrices;
    end

    outputstruct = inputstruct;
    outputstruct.visplanesfroms = visplanesfromr;
    outputstruct.vispartedgesfroms = vispartedgesfromr;
    outputstruct.soutsidemodel = routsidemodel;
    outputstruct.vispartedgesfroms_start = vispartedgesfromr_start;
    outputstruct.vispartedgesfroms_end = vispartedgesfromr_end;
    outputstruct.reftoshortlistS = reftoshortlistR;
    outputstruct.rSsho = rRsho;
    outputstruct.thetaSsho = thetaRsho;
    outputstruct.zSsho = zRsho;

    elapsedtimeSgeo = etime(clock,t00);
    elapsedtimeSRgeo = elapsedtimeSgeo;
    
    if filehandlingparameters.saveSRdatafiles == 1
        Sdata = outputstruct;
    	eval(['save(''',desiredname,''',''Sdata'',''EDinputdatahash'',''elapsedtimeSgeo'');'])
    end   
    
end

