function [tfdirect,tfgeom,tfdiff,timingdata,elapsedtimemaketfs,existingfilename] = ...
    EDmakefirstordertfs(firstorderpathdata,planedata,edgedata,Sdata,...
    Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters)
% function [tfdirect,tfgeom,tfdiff,timingdata,outpar,existingfilename] = EDmakefirstordertfs(firstorderpathdata,...
%     frequencies,Rstart,difforder,envdata,Sdata,receivers,edgedata,EDversionnumber,inpar1,inpar2,inpar3)
% EDmakefirstordertfs calculates the direct sound, specular reflection, and
% first-order diffraction.
% 
% Input parameters:
%   firstorderpathdata      Struct generated by EDfindconvexGApaths
%   planedata,edgedata,Sdata,Rdata,envdata,controlparameters   structs
%   EDversionnumber 
%   filehandlingparameters  filehandlingparameters is a struct which
%                           contains the field showtext.
% 
% Output parameters:
%   tfdirect,tfgeom,tfdiff  Matrices, size [nfrequencies,nreceivers,nsources]
%                           (if doaddsources = 0) or [nfrequencies,nreceivers]
%                           (if doaddsources = 1)
%   timingdata              Vector, [1,3], containing times for the direct
%                           sound, spec. reflections, and first-order diffraction
%                           component generations
%   elapsedtimemaketfs
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
% Uses functions EDcoordtrans2, EDwedge1st_fd, EDrecycleresultfiles from EDtoolbox
% Uses function DataHash form Matlab Central
% 
% Peter Svensson 29 Nov. 2023 (peter.svensson@ntnu.no)
%
% [tfdirect,tfgeom,tfdiff,timingdata,elapsedtimemaketfs,existingfilename] = ...
% EDmakefirstordertfs(firstorderpathdata,planedata,edgedata,Snewdata,...
% Rnewdata,envdata,controlparameters,EDversionnumber,filehandlingparameters);

% 12 Jan. 2018 First complete version. Much simplified version of the
%                           previous ESIE2maketfs. Edgehits not handled
%                           yet.
% 15 Jan. 2018 Took the direct sound and spec refl amplitude 
%              (1, 0.5, 0.25) into account)
% 17 Jan. 2018 Had forgotten the Rstart factor for the direct sound and
% specular reflection.
% 17 Jan. 2018 Took the input parameter difforder into account
% 17 Jan 2018 Added showtext as input parameter. Fixed a bug which gave the
% wrong specular reflection amplitude.
% 18 Jan 2018 Changed input parameter to the Sdata struct, to also get
% the .sourceamplitudes field. Implemented the sourceamplitudes scale
% factor.
% 31 Jan 2018 Corrected the handling of sourceamplitudes (the frequency
% dependence was not handled before).
% 8 Feb 2018 Introduced the EDinputdatastruct. Changed input parameter from
% controlparameters to three of its fields.
% 12 Feb 2018 Adapted to a change in the EDwedge1st_fd (new output
% parameter).
% 25 Aug 2021 Introduced the new hash parameter calcfirstorderdiff.
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% Added the input parameter planedata, for future handling of piston
% sources.
% 3 Oct. 2023 Added the EDsettingshash as input parameter
% 27 Oct 2023 Substantial change to input parameter list: the whole structs
% 30 Oct. 2023 Fine-tuned the EDinputdatahash
% 30 Oct. 2023 Introduced the piston source
% 1 Nov. 2023 Finished the piston source for the direct sound
% 9 Nov. 2023 Added the piston gauss order and the source amplitudes to the
% EDinputdatahash.
% 10 Nov. 2023 For the diffraction, an unnecessary recalculation of the
% edge-related coordinates was done. The values were already computed and
% stored in the rSsho, thetaSsho, zSsho etc with the reftoshortlistS etc.
% 20 Nov. 2023 Same as on 10 Nov., but now for piston sources
% 24 Nov. 2023 For piston sources, a check is made if the direct sound is
% identical for many components. Major time savings for mesh geometries.
% 29 Nov. 2023 Adapted to the name change of the field pistongausspoints to
% pistongaussorder. Fixed a small bug related to the size of a matrix which
% generates a short list for polygon pistons.

t00 = clock;

showtext = filehandlingparameters.showtext;

difforder = controlparameters.difforder;
frequencies = controlparameters.frequencies;
Rstart = controlparameters.Rstart;
directsound = controlparameters.directsound;

% New parameter for the hash 25 Aug. 2021: calcfirstorderdiff
calcfirstorderdiff = double(difforder > 0);
    
EDinputdatastruct = struct('corners',planedata.corners,'planecorners',...
    planedata.planecorners,'offedges',edgedata.offedges,...
    'sources',Sdata.coordinates,'Snedgesubs',Sdata.nedgesubs,...
    'doallSRcombinations',Sdata.doallSRcombinations,...
    'sourcetype',Sdata.sourcetype,...
    'pistoncornercoordinates',Sdata.pistoncornercoordinates,...
    'pistoncornernumbers',Sdata.pistoncornernumbers,...
    'pistongaussorder',Sdata.pistongaussorder,...
    'receivers',Rdata.coordinates,'Rnedgesubs',Rdata.nedgesubs,...
    'cair',envdata.cair,'frequencies',frequencies,'Rstart',Rstart,...
     'difforder',difforder,'directsound',directsound,'EDversionnumber',EDversionnumber);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. 

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else    
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_tf',EDinputdatahash);    
end

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'];

if foundmatch == 1
%        eval(['load ''',existingfilename,''''])
    eval(['load(''',existingfilename,''')'])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    elapsedtimemaketfs_new = etime(clock,t00);
    elapsedtimemaketfs = [elapsedtimemaketfs_new elapsedtimemaketfs];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No existing file can be reused

timingdata = zeros(1,3);

nfrequencies = length(frequencies);
nreceivers = size(Rdata.coordinates,1);
nsources = size(Sdata.coordinates,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct sound

if showtext >= 2
   disp(['      Generating direct sound components']) 
end
t00sub = clock;

if Sdata.doaddsources == 1
    tfdirect = zeros(nfrequencies,nreceivers);
else
    tfdirect = zeros(nfrequencies,nreceivers,nsources);
end

kvec = 2*pi*frequencies(:)/envdata.cair;

if firstorderpathdata.ncomponents(1) > 0
    ncomponents = size(firstorderpathdata.directsoundlist(:,1),1);

    distvecs = Rdata.coordinates(firstorderpathdata.directsoundlist(:,2),:) - ... 
            Sdata.coordinates(firstorderpathdata.directsoundlist(:,1),:);

    %     Moved to inside the monopole if-test
% % %     % alldists will be a matrix of size [nreceivers,nsources]
% % %     if ncomponents == 1
% % %         alldists = norm(distvecs);
% % %     else
% % %         alldists = sqrt( sum(distvecs.^2,2) ); 
% % %     end
    
    maxrecnumber = max( firstorderpathdata.directsoundlist(:,2) );   

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Monopole

     if strcmp(Sdata.sourcetype,'monopole') == 1

         % alldists will be a matrix of size [nreceivers,nsources]
        if ncomponents == 1
            alldists = norm(distvecs);
        else
            alldists = sqrt( sum(distvecs.^2,2) ); 
        end

         if ncomponents > nfrequencies && Sdata.doaddsources == 1
            for kk = 1:nfrequencies  
                alltfs = exp(-1i*kvec(kk)*(alldists-Rstart))./alldists...
                    .*firstorderpathdata.directsoundlist(:,3).*Sdata.sourceamplitudes( firstorderpathdata.directsoundlist(:,1),kk );
                %            if Sdata.doaddsources == 1
                tfdirect(kk,1:maxrecnumber) = accumarray(firstorderpathdata.directsoundlist(:,2),alltfs);
                %            else
                %               tfdirect(ii,firstorderpathdata.directsoundlist(:,2),firstorderpathdata.directsoundlist(:,1)) = alltfs;
                %            end
            end
        else   % if ncomponents > nfrequencies && Sdata.doaddsources == 1
            for ii = 1:ncomponents 
                if nfrequencies > 1    
                    alltfs = exp(-1i*kvec*(alldists(ii)-Rstart))./alldists(ii)...
                    .*firstorderpathdata.directsoundlist(ii,3).*(Sdata.sourceamplitudes( firstorderpathdata.directsoundlist(ii,1),: ).');
                else
                   alltfs = exp(-1i*kvec*(alldists(ii)-Rstart))./alldists(ii)...
                        .*firstorderpathdata.directsoundlist(ii,3).*(Sdata.sourceamplitudes( firstorderpathdata.directsoundlist(ii,1),: ));         
                end
                if Sdata.doaddsources == 1
                    tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) = ...
                        tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) + alltfs;
                else
                    tfdirect(:,firstorderpathdata.directsoundlist(ii,2),firstorderpathdata.directsoundlist(ii,1)) = alltfs;                
                end
            end
        end  % if ncomponents > nfrequencies && Sdata.doaddsources == 1

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Polygonpiston

     elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
        if Rstart ~= 0
            disp(['WARNING! Rstart is not zero. Are you sure that is what you want, for the polygon piston source?'])
            pause
        end

        pistonmethod = 'new';

        % We calculate all the vertical distances from the piston plane to
        % all the receivers = zreceiver.
        %           zreceiver = r*cos(theta)
        % where cos(theta) = dot( rvec,nvec)/(norm(rvec)*norm(nvec)).
        % To define rvec, the midpoint of the piston can be used. 
        % So, zreceiver = dot(rvec,nvec)
        %               = dot( distvecs(ii,:),planedata.planeeqs(pistonplane(ii),1:3) );
        %               = the expression below
        zreceiver = sum(distvecs.*planedata.planeeqs(Sdata.pistonplanes...
            (firstorderpathdata.directsoundlist(:,1)),1:3),2 );
        clear distvecs
         
        % We also calculate all the pistonplane receiver projection
        % positions = Rplane. This projection is (with receiver = Rpos)
        %             Rplane = Rpos + t*nvec
        % Then Rplane and Rpos and Spos (midpoint of piston)
        % form a right-angled triangle:
        %       dot(Rpos - Rplane,Spos - Rplane) = 0 ->
        %       dot(-t*nvec,Spos - Rpos -t*nvec) = 0 ->
        %       -t*dot(nvec,Spos) + t*dot(nvec,Rpos) +
        %       t^2*dot(nvec,nvec) = 0 ->
        %       t^2 = t*dot(nvec,Spos-Rpos) -> 
        %       t = dot(nvec,Spos-Rpos) ->
        %
        %   Rplane = Rpos + nvec*dot(nvec,Spos - Rpos)
        %          = Rpos - nvec*dot(nvec,Rpos - Spos)
        %          = Rpos - nvec*zreceiver

        longlistpistonnumber = firstorderpathdata.directsoundlist(:,1);
        pistonplanenumber = Sdata.pistonplanes(longlistpistonnumber);
        Rplane = Rdata.coordinates(firstorderpathdata.directsoundlist(:,2),:) - ...
        zreceiver.*planedata.planeeqs(pistonplanenumber,1:3);

        if pistonmethod == 'new'

            % Addition 24 Nov. 2023: identify identical piston-to-receiver
            % geometries and create a shortlist. Combinations where abs(dx),
            % abs(dy), abs(dz) and piston-dx, piston-dy, piston-dz are identical can use the same
            % tf. dx = distance between piston center point and receiver point.
    
            longlistrecnumber = firstorderpathdata.directsoundlist(:,2);
    
            dx = Sdata.coordinates(longlistpistonnumber,1) - Rdata.coordinates(longlistrecnumber,1);
            dy = Sdata.coordinates(longlistpistonnumber,2) - Rdata.coordinates(longlistrecnumber,2);
            dz = Sdata.coordinates(longlistpistonnumber,3) - Rdata.coordinates(longlistrecnumber,3);
            npistons = size(Sdata.pistoncornernumbers,1);
            pistonsizes = zeros(npistons,3);
            for jj = 1:npistons
                pistonsizes(jj,:) = range(Sdata.pistoncornercoordinates(Sdata.pistoncornernumbers(jj,:),:));
            end
            pistonsizes = round(pistonsizes*1e5)/1e5;
    
% %             A = [round( abs(dx(:))*1e5)/1e5 round(abs(dy(:))*1e5)/1e5 round(abs(dz(:))*1e5)/1e5 pistonsizes(longlistpistonnumber)];
% %             [shortlist_pistontorec,examplefromshortlist_pistontorec,reftoshortlist_pistontorec] = unique(A,'rows');
            [shortlist_pistontorec,examplefromshortlist_pistontorec,reftoshortlist_pistontorec] = ...
                unique([round( abs(dx(:))*1e5)/1e5 round(abs(dy(:))*1e5)/1e5 round(abs(dz(:))*1e5)/1e5 pistonsizes(longlistpistonnumber,:)],'rows');
  
            [sortedvec,isort] = sort(reftoshortlist_pistontorec);
    
%             dsort = diff(sortedvec);
%             istep = find(dsort);
            istep = find(diff(sortedvec));
            listofalloccurences = cell(length(istep)+1,1);
            for jj = 1:length(istep)+1
                if jj == 1
                    istart = 1;
                else
                    istart = istep(jj-1)+1;
                end
                if jj == length(istep)+1
                   iend = length(sortedvec);
                else
                   iend = istep(jj); 
                end
            
                listofalloccurences{jj} = uint32(isort(istart:iend));        
            end
    
            nshortlistcomponents = size(shortlist_pistontorec,1);
            if showtext >= 2
                disp(['      Calculating tfs for ',int2str(nshortlistcomponents),' unique piston-to-points tfs instead of all ',int2str(ncomponents)])
            end
    
            for nn = 1:nshortlistcomponents
                examplepistonnumber = longlistpistonnumber(examplefromshortlist_pistontorec(nn));
                examplerecnumber = longlistrecnumber(examplefromshortlist_pistontorec(nn));
                npistonedges = Sdata.ncornersperpiston(examplepistonnumber);
                pistonplane = Sdata.pistonplanes(examplepistonnumber);
                        
                checkinsidehit = 1;
                checkedgehit = 0;
                lineintsummation = zeros(length(kvec),1);
                for jj = 1:npistonedges
    %                        disp(['           Edge ',int2str(jj)])
                    startnumber = Sdata.pistoncornernumbers(examplepistonnumber,jj);
                    if jj < Sdata.ncornersperpiston(examplepistonnumber)
                        endnumber = Sdata.pistoncornernumbers(examplepistonnumber,jj+1);
                    else
                        endnumber = Sdata.pistoncornernumbers(examplepistonnumber,1);
                    end
                    c1 = Sdata.pistoncornercoordinates(startnumber,:);
                    c2 = Sdata.pistoncornercoordinates(endnumber,:);
                    edgelength = norm( c2-c1 );
                    Rplaneshift = Rplane(examplefromshortlist_pistontorec(nn),:) - c1;
                    Rshift = Rdata.coordinates(examplerecnumber,:) - c1;
                    rotmat = Sdata.pistonrotationmatrices{Sdata.reftorotationmatrices(examplepistonnumber,jj)};
                    Rplanemod = rotmat*Rplaneshift(:);
                    Rmod = rotmat*Rshift(:);
                    if Rplanemod(2) < 0
                        checkinsidehit = 0;
                    elseif Rplanemod(2) == 0
                        checkedgehit = 1;
                    end
                    % nvecmod will be [0 0 1]
                    % Rplanemod will always have z = 0
                    % For inside-piston receivers, Rplanemod(2) is
                    % positive for all edges.
                    % Rmod will have z = zreceiver which was calculated
                    % further up (maybe it wasn't needed up there and 
                    % can be skipped?).
                    % Integration: x will be the axis along the edge;
                    % range: 0,edgelength.
                    % dist to rec = sqrt( (x-Rmod(1)).^2 + Rmod(2)^2 + Rmod(3)^2 )
                    % psi = Rplanemod(2)
    
                    for kk = 1:length(kvec)
                        intval = integral(@(x) EDpistonedgeintegrand(x,Rplanemod(2)^2,Rmod(3)^2,kvec(kk)),0-Rplanemod(1),edgelength-Rplanemod(1),'RelTol',1e-8);
                        lineintsummation(kk) = lineintsummation(kk) + intval*Rplanemod(2);
                    end
                    
                end
                alltfs = -lineintsummation/2/pi;
     
                if checkinsidehit == 1
                     alltfs = alltfs + exp(-1i*kvec*(Rmod(3)-Rstart));
                end
    
                alltfs = alltfs.*(-2i)*envdata.cair/Sdata.pistonareas(examplepistonnumber)./frequencies(:);
             
                if Sdata.doaddsources == 1
                    % old pistonmethod:
%                     tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) = ...
%                       tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) + alltfs;
                    for mm = 1:length(listofalloccurences{nn})
                        rep_recno = longlistrecnumber(listofalloccurences{nn}(mm));
                        tfdirect(:,rep_recno) = tfdirect(:,rep_recno) + alltfs;
                    end
                else    
                    % old pistonmethod:
     %               tfdirect(:,firstorderpathdata.directsoundlist(ii,2),firstorderpathdata.directsoundlist(ii,1)) = alltfs;

                    rep_recno = longlistrecnumber(listofalloccurences{nn});
                    rep_pistonno = longlistpistonnumber(listofalloccurences{nn});
%                     savetemp
%                     pause
                    for nf = 1:nfrequencies
                        addmatrix = accumarray([rep_recno rep_pistonno],alltfs(nf));
                        [nadd1,nadd2] = size(addmatrix);
                        tfdirect(nf,1:nadd1,1:nadd2) = tfdirect(nf,1:nadd1,1:nadd2) + reshape(addmatrix,1,nadd1,nadd2);
                    end
    
                end   % if Sdata.doaddsources == 1
            
            end     % for nn = 1:nshortlistcomponents
        end   % if pistonmethod == 'new'

        if pistonmethod == 'old'
            for ii = 1:ncomponents 
                pistonnumber = firstorderpathdata.directsoundlist(ii,1);
                recnumber = firstorderpathdata.directsoundlist(ii,2);
                npistonedges = Sdata.ncornersperpiston(pistonnumber);
                pistonplane = Sdata.pistonplanes(pistonnumber);
                        
                checkinsidehit = 1;
                checkedgehit = 0;
                lineintsummation = zeros(length(kvec),1);
                for jj = 1:npistonedges
    %                        disp(['           Edge ',int2str(jj)])
                    startnumber = Sdata.pistoncornernumbers(pistonnumber,jj);
                    if jj < Sdata.ncornersperpiston(pistonnumber)
                        endnumber = Sdata.pistoncornernumbers(pistonnumber,jj+1);
                    else
                        endnumber = Sdata.pistoncornernumbers(pistonnumber,1);
                    end
                    c1 = Sdata.pistoncornercoordinates(startnumber,:);
                    c2 = Sdata.pistoncornercoordinates(endnumber,:);
                    edgelength = norm( c2-c1 );
                    Rplaneshift = Rplane(ii,:) - c1;
                    Rshift = Rdata.coordinates(firstorderpathdata.directsoundlist(ii,2),:) - c1;
                    rotmat = Sdata.pistonrotationmatrices{Sdata.reftorotationmatrices(pistonnumber,jj)};
                    Rplanemod = rotmat*Rplaneshift(:);
                    Rmod = rotmat*Rshift(:);
                    if Rplanemod(2) < 0
                        checkinsidehit = 0;
                    elseif Rplanemod(2) == 0
                        checkedgehit = 1;
                    end
                    % nvecmod will be [0 0 1]
                    % Rplanemod will always have z = 0
                    % For inside-piston receivers, Rplanemod(2) is
                    % positive for all edges.
                    % Rmod will have z = zreceiver which was calculated
                    % further up (maybe it wasn't needed up there and 
                    % can be skipped?).
                    % Integration: x will be the axis along the edge;
                    % range: 0,edgelength.
                    % dist to rec = sqrt( (x-Rmod(1)).^2 + Rmod(2)^2 + Rmod(3)^2 )
                    % psi = Rplanemod(2)
    
                    for kk = 1:length(kvec)
                        intval = integral(@(x) EDpistonedgeintegrand(x,Rplanemod(2)^2,Rmod(3)^2,kvec(kk)),0-Rplanemod(1),edgelength-Rplanemod(1),'RelTol',1e-8);
                        lineintsummation(kk) = lineintsummation(kk) + intval*Rplanemod(2);
                    end
                end
                alltfs = -lineintsummation/2/pi;
                              
                if checkinsidehit == 1
                     alltfs = alltfs + exp(-1i*kvec*(Rmod(3)-Rstart));
                end
    
                alltfs = alltfs.*(-2i)*envdata.cair/Sdata.pistonareas(pistonnumber)./frequencies(:);
                
               if Sdata.doaddsources == 1
                  tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) = ...
                      tfdirect(:,firstorderpathdata.directsoundlist(ii,2)) + alltfs;
               else             
                  tfdirect(:,firstorderpathdata.directsoundlist(ii,2),firstorderpathdata.directsoundlist(ii,1)) = alltfs;
               end
            end    % for ii = 1:ncomponents 

        end    % if pistonmethod == 'old'

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Only monopole and polygonpiston have been implemented

     else
        error('ERROR: Only monopole and polygonpiston sourcetypes have been implemented.')
     end

end   % if firstorderpathdata.ncomponents(1) > 0


timingdata(1) = etime(clock,t00sub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specular reflections

if showtext >= 2
   disp(['      Generating specular reflection components']) 
end

t00sub = clock;

if Sdata.doaddsources == 1
    tfgeom = zeros(nfrequencies,nreceivers);
else
    tfgeom = zeros(nfrequencies,nreceivers,nsources);
end

if firstorderpathdata.ncomponents(2) > 0 & strcmp(Sdata.sourcetype,'monopole')

    ncomponents = size(firstorderpathdata.specrefllist(:,1),1);

    distvecs = firstorderpathdata.specreflIScoords - ...
        Rdata.coordinates(firstorderpathdata.specrefllist(:,2),:);

    if ncomponents == 1
       alldists = norm(distvecs);
    else
        alldists = sqrt( sum(distvecs.^2,2) );
    end

    maxrecnumber = max( firstorderpathdata.specrefllist(:,2) );

    if ncomponents > nfrequencies && Sdata.doaddsources == 1
        for ii = 1:nfrequencies    
                alltfs = exp(-1i*kvec(ii)*(alldists-Rstart))./alldists...
                    .*firstorderpathdata.specrefllist(:,3).*(Sdata.sourceamplitudes( firstorderpathdata.specrefllist(:,1),ii ));                
                tfgeom(ii,1:maxrecnumber) = accumarray(firstorderpathdata.specrefllist(:,2),alltfs);
        end
    else
        for ii = 1:ncomponents 
            alltfs = exp(-1i*kvec*(alldists(ii)-Rstart))./alldists(ii)...
                .*firstorderpathdata.specrefllist(ii,3).*(Sdata.sourceamplitudes( firstorderpathdata.specrefllist(ii,1),: ).');
            
           if Sdata.doaddsources == 1
              tfgeom(:,firstorderpathdata.specrefllist(ii,2)) = ...
                  tfgeom(:,firstorderpathdata.specrefllist(ii,2)) + alltfs;
           else
              tfgeom(:,firstorderpathdata.specrefllist(ii,2),firstorderpathdata.specrefllist(ii,1)) = alltfs;
           end
       end
       
        
    end
end

timingdata(2) = etime(clock,t00sub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffraction

t00sub = clock;

if difforder > 0

    if showtext >= 2
        disp(['      Generating first-order diffraction components']) 
    end
    
    if Sdata.doaddsources == 1
        tfdiff = zeros(nfrequencies,nreceivers);
    else
        tfdiff = zeros(nfrequencies,nreceivers,nsources);
    end

    iv = find(firstorderpathdata.edgeisactive);

    sou_vs_edges = sign(squeeze(sum(firstorderpathdata.diffpaths,1)));
    if size(sou_vs_edges,2) == 1
        sou_vs_edges = sou_vs_edges.';
    end
    rec_vs_edges = sign(squeeze(sum(firstorderpathdata.diffpaths,2)));
    if size(rec_vs_edges,2) == 1
        rec_vs_edges = rec_vs_edges.';
    end

    for ii = 1:length(iv)
        cylcoordS = zeros(nsources,3);
        cylcoordR = zeros(nreceivers,3);

        edgenumber = iv(ii);
        if showtext >= 2
            disp(['         Edge no. ',int2str(edgenumber)])
        end
        edgecoords = [edgedata.edgestartcoords(edgenumber,:);edgedata.edgeendcoords(edgenumber,:)];

        sourceandreceivercombos = squeeze(firstorderpathdata.diffpaths(:,:,edgenumber));
        iv2 = find(sourceandreceivercombos);
        [Rnumber,Snumber] = ind2sub([nreceivers,nsources],iv2);

        ivS = find(sou_vs_edges(:,edgenumber));
        ivR = find(rec_vs_edges(:,edgenumber));
 
        if strcmp(Sdata.sourcetype,'monopole')
 
            % 10 Nov. 2023 Instead of recalculating the edge-related
            % coordinates, they are read from the available shortlists
%            [rs,thetas,zs,rr,thetar,zr] = EDcoordtrans2(Sdata.coordinates(ivS,:),Rdata.coordinates(ivR,:),edgecoords,edgedata.edgenvecs(edgenumber,:));
            rs     = Sdata.rSsho(Sdata.reftoshortlistS(edgenumber,ivS));
            thetas = Sdata.thetaSsho(Sdata.reftoshortlistS(edgenumber,ivS));
            zs     = Sdata.zSsho(Sdata.reftoshortlistS(edgenumber,ivS));

            rr     = Rdata.rRsho(Rdata.reftoshortlistR(edgenumber,ivR));
            thetar = Rdata.thetaRsho(Rdata.reftoshortlistR(edgenumber,ivR));
            zr     = Rdata.zRsho(Rdata.reftoshortlistR(edgenumber,ivR));


            cylcoordS(ivS,:) = [rs thetas zs];
            cylcoordR(ivR,:) = [rr thetar zr];
    
            for jj = 1:length(iv2)
               if showtext >= 2
                    disp(['   Source no. ',int2str(Snumber(jj)),', Rec. no. ',int2str(Rnumber(jj))])
               end
                [tfnew,singularterm,zfirst] = EDwedge1st_fd(envdata.cair,frequencies,edgedata.closwedangvec(edgenumber),...
                    cylcoordS(Snumber(jj),1),cylcoordS(Snumber(jj),2),cylcoordS(Snumber(jj),3),...
                    cylcoordR(Rnumber(jj),1),cylcoordR(Rnumber(jj),2),cylcoordR(Rnumber(jj),3),...
                    edgedata.edgelengthvec(edgenumber)*[0 1],'n',Rstart,[1 1]);  
    
                tfnew = tfnew.*(Sdata.sourceamplitudes( Snumber(jj),: ).');

                if Sdata.doaddsources == 1
                    tfdiff(:,Rnumber(jj)) =  tfdiff(:,Rnumber(jj)) + tfnew;                       
                else
                    tfdiff(:,Rnumber(jj),Snumber(jj)) =  tfdiff(:,Rnumber(jj),Snumber(jj)) + tfnew;           
                end
            end
        else % polygonpiston
            npistonpoints = size(Sdata.pistongausscoordinates{ivS},1);

            % 20 Nov. 2023 Use the already existing short lists of
            % edge-related coordinates instead of recalculating them
%            [rr,thetar,zr] = EDcoordtrans1(Rdata.coordinates(ivR,:),edgecoords,edgedata.edgenvecs(edgenumber,:));
            rr = Rdata.rRsho(Rdata.reftoshortlistR(edgenumber,ivR));
            thetar = Rdata.thetaRsho(Rdata.reftoshortlistR(edgenumber,ivR));
            zr = Rdata.zRsho(Rdata.reftoshortlistR(edgenumber,ivR));

            cylcoordR(ivR,:) = [rr(:) thetar(:) zr(:)];
                
            for jj = 1:length(iv2)

                tfsumpiston = zeros(length(frequencies),1); 
                for np = 1:npistonpoints

                    % 20 Nov. 2023 Use the already existing short lists of
                    % edge-related coordinates instead of recalculating them
%                    [rs,thetas,zs] = EDcoordtrans1(Sdata.pistongausscoordinates{ivS}(np,:),edgecoords,edgedata.edgenvecs(edgenumber,:));
                    rs     = Sdata.rSsho(Sdata.reftoshortlistS(edgenumber,ivS,np));
                    thetas = Sdata.thetaSsho(Sdata.reftoshortlistS(edgenumber,ivS,np));
                    zs     = Sdata.zSsho(Sdata.reftoshortlistS(edgenumber,ivS,np));

                    [tfonepistonpoint,singularterm,zfirst] = EDwedge1st_fd(envdata.cair,frequencies,edgedata.closwedangvec(edgenumber),...
                        rs(:),thetas(:),zs(:    ),cylcoordR(Rnumber(jj),1),cylcoordR(Rnumber(jj),2),cylcoordR(Rnumber(jj),3),...
                        edgedata.edgelengthvec(edgenumber)*[0 1],'n',Rstart,[1 1]);  
                    tfsumpiston = tfsumpiston + tfonepistonpoint*Sdata.pistongaussweights{ivS}(np);
                end

                tfnew = tfsumpiston.*(Sdata.sourceamplitudes( Snumber(jj),: ).');
                if Sdata.doaddsources == 1
                    tfdiff(:,Rnumber(jj)) =  tfdiff(:,Rnumber(jj)) + tfnew;                       
                else
                    tfdiff(:,Rnumber(jj),Snumber(jj)) =  tfdiff(:,Rnumber(jj),Snumber(jj)) + tfnew;           
                end
            end


%            pause

        end        

    end
else
   tfdiff = zeros(size(tfgeom));
end

timingdata(3) = etime(clock,t00sub);

elapsedtimemaketfs = etime(clock,t00);

eval(['save(''',desiredname,''',''tfdirect'',''tfgeom'',''tfdiff'',''timingdata'',''EDinputdatahash'',''elapsedtimemaketfs'');'])


