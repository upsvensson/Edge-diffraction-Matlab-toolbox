function [irdirect,irgeom,irdiff,timingdata,elapsedtimemakeirs,existingfilename] ...
= EDmakefirstorderirs(firstorderpathdata,planedata,edgedata,Sdata,...
  Rdata,envdata,controlparameters,EDversionnumber,...
  filehandlingparameters)
% EDmakefirstorderirs - Computes the direct sound, specular reflection and first-order diffraction IRs.
% The impulse responses are computed for the paths in firstorderpathdata,
% with the sampling frequency controlparameters.fs. First-order
% diffraction is computed by EDwedge1st_ir. The results are saved in the
% stage file <filestem>_ir.mat, which is reused when the input data are
% the same.
%
% Input parameters:
%   firstorderpathdata, planedata, edgedata, Sdata, Rdata
%                   structs, see help EDdatastructs (planedata is not
%                   used)
%   envdata         struct; the field .cair is used.
%   controlparameters  struct; the fields .fs, .Rstart, .difforder
%                   (whether it is larger than 0), .directsound and
%                   .saveindividualfirstdiff are used.
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem, .showtext and
%                   .suppressresultrecycling are used.
%
% Output parameters:
%   irdirect, irgeom, irdiff  Matrices with the direct sound, the specular
%                   reflections and the first-order diffraction, of size
%                   [nsamples,nreceivers] if the sources are added
%                   (Sdata.doaddsources = 1), otherwise
%                   [nsamples,nreceivers,nsources].
%                   If controlparameters.saveindividualfirstdiff = 1, irdiff
%                   is instead a cell array, irdiff{ireceiver,isource}, where
%                   each cell has the fields .irvectors [nsamples,ndiff]
%                   (one IR per edge), .edgenumbers [1,ndiff] and
%                   .locations [3,ndiff] (the point on the edge that gives
%                   the first arrival).
%   timingdata      [1,3] the times for the direct sound, the specular
%                   reflections and the first-order diffraction.
%   elapsedtimemakeirs
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [irdirect,irgeom,irdiff,timingdata,elapsedtimemakeirs,existingfilename] = ...
%    EDmakefirstorderirs(firstorderpathdata,planedata,edgedata,Sdata,...
%    Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters);

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
% 18 Jan 2018 Changed input parameter to the Sindata struct, to also get
% the .sourceamplitudes field. Implemented the sourceamplitudes scale
% factor.
% 28 Jan 2018 First version; modified from EDmakefirstordertfs
% 28 Feb 2018 Introducing into the EDtoolbox master
% 7 Apr 2018 Introduced the use of saveindividualfirstdiff
% 9 Apr 2018 Little bug fixed, on line 358, thanks to Ville Pulkki
% 12 Apr 2018 Little bug fixed, around line 334, thanks to Ville Pulkki
% 6 May 2018 Little bug fixed, around line 331, thanks to Ville Pulkki
% 25 Aug 2021 Introduced the new hash parameter calcfirstorderdiff.
% 29 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 3 Oct. 2023 Added the EDsettingshash as input parameter.
% 27 Oct. 2023 Changed the input parameter list substantially: whole
% structs.
% 30 Oct. 2023 Fine-tuned the EDinputdatahash
% 12 Oct. 2026 New EDinputdatahash: chained to the hashes of
% firstorderpathdata, edgedata, Sdata and Rdata, with a local stageversion
% instead of the overall EDversionnumber. sourceamplitudes and
% doaddsources were missing in the hash and have been added. The input
% parameter planedata is no longer used. Fixed: elapsedtimemakeirs was
% only the time for the diffraction part (t00 was reused for the section
% timers). load/save without eval.
% 14 Oct. 2026 The stageversion is taken from EDstageversions, and it
% also contains the versions of the subroutines that this function uses.
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).
% 15 Oct. 2026 The list of called functions was replaced by a See also line.

t00 = clock;
t00_total = t00;

showtext = filehandlingparameters.showtext;

fs = controlparameters.fs;
Rstart = controlparameters.Rstart;
difforder = controlparameters.difforder;
directsound = controlparameters.directsound;
saveindividualfirstdiff = controlparameters.saveindividualfirstdiff;

% New parameter for the hash 25 Aug. 2021: calcfirstorderdiff
calcfirstorderdiff = double(difforder > 0);

% The version of this function, and of the subroutines it uses, are
% defined in EDstageversions. They are part of the hash, so that saved
% stage files are not reused after a change in the computations.
stageversion = EDstageversions('EDmakefirstorderirs');

EDinputdatastruct = struct('stage','EDmakefirstorderirs','stageversion',stageversion,...
    'pathshash',EDparenthash(firstorderpathdata),'edgedatahash',EDparenthash(edgedata),...
    'Sdatahash',EDparenthash(Sdata),'Rdatahash',EDparenthash(Rdata),...
    'cair',envdata.cair,'fs',fs,'Rstart',Rstart,'calcfirstorderdiff',calcfirstorderdiff,...
    'directsound',directsound,'saveindividualfirstdiff',saveindividualfirstdiff,...
    'sourceamplitudes',Sdata.sourceamplitudes,'doaddsources',Sdata.doaddsources);
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
		'_ir',EDinputdatahash);
end

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_ir.mat'];

if foundmatch == 1
	filedata = load(existingfilename,'irdirect','irgeom','irdiff','timingdata','elapsedtimemakeirs');
	irdirect = filedata.irdirect;
	irgeom = filedata.irgeom;
	irdiff = filedata.irdiff;
	timingdata = filedata.timingdata;
	elapsedtimemakeirs = filedata.elapsedtimemakeirs;
	if ~strcmp(existingfilename,desiredname)
		copyfile(existingfilename,desiredname);
	end
	elapsedtimemakeirs_new = etime(clock,t00);
	elapsedtimemakeirs = [elapsedtimemakeirs_new elapsedtimemakeirs];
	return
end

%----------------------------------------------------------------
% No existing file can be used

timingdata = zeros(1,3);

nreceivers = size(Rdata.coordinates,1);
nsources = size(Sdata.coordinates,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specular reflections

if showtext >= 2
   disp(['      Generating specular reflection components']) 
end

t00 = clock;

if firstorderpathdata.ncomponents(2) > 0

    ncomponents = size(firstorderpathdata.specrefllist(:,1),1);

    distvecs = firstorderpathdata.specreflIScoords - ...
        Rdata.coordinates(firstorderpathdata.specrefllist(:,2),:);

    if ncomponents == 1
       alldists = norm(distvecs);
    else
        alldists = sqrt( sum(distvecs.^2,2) );
    end

    allslotnumberfracs = (alldists-Rstart)/envdata.cair*fs+1;
    allamps = 1./alldists.*firstorderpathdata.specrefllist(:,3).*...
        Sdata.sourceamplitudes( firstorderpathdata.specrefllist(:,1) );
    allslotnumbers = floor(allslotnumberfracs);
    if any(allslotnumbers<1)
        error('ERROR: Rstart is set to too low a value!')
    end
    allslotnumberfracs = allslotnumberfracs - allslotnumbers;
    nirlength = max(allslotnumbers + 2);
    
    if Sdata.doaddsources == 1
        irgeom = zeros(nirlength,nreceivers);
    else
        irgeom = zeros(nirlength,nreceivers,nsources);
    end
    
    for ii = 1:ncomponents
        if Sdata.doaddsources == 0
            irgeom(allslotnumbers(ii),firstorderpathdata.specrefllist(ii,2),...
                firstorderpathdata.specrefllist(ii,1)) = ...
                  irgeom(allslotnumbers(ii),firstorderpathdata.specrefllist(ii,2),...
                  firstorderpathdata.specrefllist(ii,1)) + allamps(ii)*(1-allslotnumberfracs(ii));
            irgeom(allslotnumbers(ii)+1,firstorderpathdata.specrefllist(ii,2),...
                firstorderpathdata.specrefllist(ii,1)) = ...
                  irgeom(allslotnumbers(ii)+1,firstorderpathdata.specrefllist(ii,2),...
                  firstorderpathdata.specrefllist(ii,1)) + allamps(ii)*(allslotnumberfracs(ii));            
        else
           irgeom(allslotnumbers(ii),firstorderpathdata.specrefllist(ii,2)) = ...
                  irgeom(allslotnumbers(ii),firstorderpathdata.specrefllist(ii,2))...
                   + allamps(ii)*(1-allslotnumberfracs(ii));
            irgeom(allslotnumbers(ii)+1,firstorderpathdata.specrefllist(ii,2)) = ...
                  irgeom(allslotnumbers(ii)+1,firstorderpathdata.specrefllist(ii,2))...
                   + allamps(ii)*(allslotnumberfracs(ii));            
        end
        
    end
  
else
    nirlength = 2;
    if Sdata.doaddsources == 1
        irgeom = zeros(nirlength,nreceivers);
    else
        irgeom = zeros(nirlength,nreceivers,nsources);
    end    
end

timingdata(2) = etime(clock,t00);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Direct sound

if showtext >= 2
   disp(['      Generating direct sound components']) 
end
t00 = clock;

% alldists, allslotnumberfracs, allamps. allslotnumbers will be matrices 
% of size [nreceivers,nsources]

if firstorderpathdata.ncomponents(1) > 0
    ncomponents = size(firstorderpathdata.directsoundlist(:,1),1);

    distvecs = Sdata.coordinates(firstorderpathdata.directsoundlist(:,1),:) - ...
        Rdata.coordinates(firstorderpathdata.directsoundlist(:,2),:);

    if ncomponents == 1
       alldists = norm(distvecs);
    else
        alldists = sqrt( sum(distvecs.^2,2) ); 
    end
    
    allslotnumberfracs = (alldists-Rstart)/envdata.cair*fs+1;
    allamps = 1./alldists.*firstorderpathdata.directsoundlist(:,3).*...
        Sdata.sourceamplitudes( firstorderpathdata.directsoundlist(:,1) );
    allslotnumbers = floor(allslotnumberfracs);
    if any(allslotnumbers<1)
        error('ERROR: Rstart is set to too low a value!')
    end
    allslotnumberfracs = allslotnumberfracs - allslotnumbers;
    nirlengthdirect = max(allslotnumbers + 2);
    
    if nirlengthdirect > nirlength
        if Sdata.doaddsources == 1
            irgeom = [irgeom;zeros(nirlengthdirect-nirlength,nreceivers)];
        else
            irgeom = [irgeom;zeros(nirlengthdirect-nirlength,nreceivers,nsources)];
        end    
        nirlength = nirlengthdirect;
    end
    
    if Sdata.doaddsources == 1
        irdirect = zeros(nirlength,nreceivers);
    else
        irdirect = zeros(nirlength,nreceivers,nsources);
    end
	    
    for ii = 1:ncomponents
        if Sdata.doaddsources == 0
            irdirect(allslotnumbers(ii),firstorderpathdata.directsoundlist(ii,2),...
                firstorderpathdata.directsoundlist(ii,1)) = ...
                  irdirect(allslotnumbers(ii),firstorderpathdata.directsoundlist(ii,2),...
                  firstorderpathdata.directsoundlist(ii,1)) + allamps(ii)*(1-allslotnumberfracs(ii));
            irdirect(allslotnumbers(ii)+1,firstorderpathdata.directsoundlist(ii,2),...
                firstorderpathdata.directsoundlist(ii,1)) = ...
                  irdirect(allslotnumbers(ii)+1,firstorderpathdata.directsoundlist(ii,2),...
                  firstorderpathdata.directsoundlist(ii,1)) + allamps(ii)*(allslotnumberfracs(ii));            
        else
           irdirect(allslotnumbers(ii),firstorderpathdata.directsoundlist(ii,2)) = ...
                  irdirect(allslotnumbers(ii),firstorderpathdata.directsoundlist(ii,2))...
                   + allamps(ii)*(1-allslotnumberfracs(ii));
            irdirect(allslotnumbers(ii)+1,firstorderpathdata.directsoundlist(ii,2)) = ...
                  irdirect(allslotnumbers(ii)+1,firstorderpathdata.directsoundlist(ii,2))...
                   + allamps(ii)*(allslotnumberfracs(ii));            
        end
        
    end

else
    if Sdata.doaddsources == 1
        irdirect = zeros(nirlength,nreceivers);
    else
        irdirect = zeros(nirlength,nreceivers,nsources);
    end  
end

timingdata(1) = etime(clock,t00);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diffraction

t00 = clock;

if difforder > 0

    if showtext >= 2
        disp(['      Generating first-order diffraction components']) 
    end
    
    if saveindividualfirstdiff == 0
        if Sdata.doaddsources == 1
            irdiff = zeros(nirlength,nreceivers);
        else
            irdiff = zeros(nirlength,nreceivers,nsources);
        end
    else
        if Sdata.doaddsources == 1
            irdiff = cell(nreceivers,1);
            for ii = 1:nreceivers
               irdiff{ii,1} = struct; 
               irdiff{ii,1}.irvectors = zeros(nirlength,1); 
               irdiff{ii,1}.locations = zeros(3,1);
               irdiff{ii,1}.edgenumbers = zeros(1,1);
            end
        else
            irdiff = cell(nreceivers,nsources);
            for ii = 1:nreceivers
                for jj = 1:nsources
                    irdiff{ii,jj} = struct;
                    irdiff{ii,jj}.irvectors = zeros(nirlength,1); 
                    irdiff{ii,jj}.locations = zeros(3,1);
                    irdiff{ii,jj}.edgenumbers = zeros(1,1);
                end                
            end
        end 
        componentcounter = zeros(nreceivers,nsources);
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
        edgecoords = [edgedata.edgestartcoords(edgenumber,:);edgedata.edgeendcoords(edgenumber,:)];

        sourceandreceivercombos = squeeze(firstorderpathdata.diffpaths(:,:,edgenumber));
        iv2 = find(sourceandreceivercombos);
        [Rnumber,Snumber] = ind2sub([nreceivers,nsources],iv2);

        ivS = find(sou_vs_edges(:,edgenumber));
        ivR = find(rec_vs_edges(:,edgenumber));
        [rs,thetas,zs,rr,thetar,zr] = EDcoordtrans2(Sdata.coordinates(ivS,:),Rdata.coordinates(ivR,:),edgecoords,edgedata.edgenvecs(edgenumber,:));

        cylcoordS(ivS,:) = [rs thetas zs];
        cylcoordR(ivR,:) = [rr thetar zr];

        for jj = 1:length(iv2)
            [irnew,singularterm,singularterm,zfirst] = EDwedge1st_ir(fs,envdata.cair,edgedata.closwedangvec(edgenumber),...
                cylcoordS(Snumber(jj),1),cylcoordS(Snumber(jj),2),cylcoordS(Snumber(jj),3),...
                cylcoordR(Rnumber(jj),1),cylcoordR(Rnumber(jj),2),cylcoordR(Rnumber(jj),3),...
                edgedata.edgelengthvec(edgenumber)*[0 1],'n',Rstart,[1 1]);  

            irnew = irnew.*Sdata.sourceamplitudes( Snumber(jj) );
            nirlengthnew = length(irnew);

            if nirlengthnew > nirlength
                if saveindividualfirstdiff == 0
                    if Sdata.doaddsources == 1
                        irgeom   = [irgeom;  zeros(nirlengthnew-nirlength,nreceivers)];
                        irdirect = [irdirect;zeros(nirlengthnew-nirlength,nreceivers)];
                        irdiff = [irdiff;zeros(nirlengthnew-nirlength,nreceivers)];                    
                    else
                        irgeom   = [irgeom;  zeros(nirlengthnew-nirlength,nreceivers,nsources)];
                        irdirect = [irdirect;zeros(nirlengthnew-nirlength,nreceivers,nsources)];
                        irdiff = [irdiff;zeros(nirlengthnew-nirlength,nreceivers,nsources)];
                    end    
                else
                    if Sdata.doaddsources == 1
                        irgeom   = [irgeom;  zeros(nirlengthnew-nirlength,nreceivers)];
                        irdirect = [irdirect;zeros(nirlengthnew-nirlength,nreceivers)];
                        % Bug fixed 6 May 2018: before, the ncolumns in the
                        % loop was a fixed value 1, which gave errors when
                        % there was only one source.
                        for kk = 1:nreceivers
                            ncolumns = size(irdiff{kk,1}.irvectors,2);
                            irdiff{kk,1}.irvectors = [irdiff{kk,1}.irvectors;zeros(nirlengthnew-nirlength,ncolumns)]; 
                        end
                    else
                        irgeom   = [irgeom;  zeros(nirlengthnew-nirlength,nreceivers,nsources)];
                        irdirect = [irdirect;zeros(nirlengthnew-nirlength,nreceivers,nsources)];
                        for kk = 1:nreceivers
                            for mm = 1:nsources
                                ncolumns = size(irdiff{kk,mm}.irvectors,2);
                                irdiff{kk,mm}.irvectors = [irdiff{kk,mm}.irvectors;zeros(nirlengthnew-nirlength,ncolumns)];                
                            end                            
                        end
                    end                        
                end
                nirlength = nirlengthnew;                
            end    
            if saveindividualfirstdiff == 0    
                if Sdata.doaddsources == 1
                    irdiff(1:nirlengthnew,Rnumber(jj)) =  irdiff(1:nirlengthnew,Rnumber(jj)) + irnew;                       
                else
                    irdiff(1:nirlengthnew,Rnumber(jj),Snumber(jj)) =  irdiff(1:nirlengthnew,Rnumber(jj),Snumber(jj)) + irnew;           
                end
            else
                componentcounter(Rnumber(jj),Snumber(jj)) = componentcounter(Rnumber(jj),Snumber(jj)) + 1;
                currentcounter = componentcounter(Rnumber(jj),Snumber(jj));
                edgestart = edgedata.edgestartcoords(edgenumber,:);
                edgeend   = edgedata.edgeendcoords(edgenumber,:);                
                edgevector = edgeend - edgestart;
                xyzfirst = edgestart + zfirst*edgevector/norm(edgevector);
                
                if Sdata.doaddsources == 1
                    irdiff{Rnumber(jj),1}.irvectors(1:nirlengthnew,currentcounter) =  irnew;  
                    irdiff{Rnumber(jj),1}.edgenumbers(currentcounter) = edgenumber;
                    irdiff{Rnumber(jj),1}.locations(:,currentcounter) = xyzfirst(:);
                else
                    irdiff{Rnumber(jj),Snumber(jj)}.irvectors(1:nirlengthnew,currentcounter) =  irnew;           
                    irdiff{Rnumber(jj),Snumber(jj)}.edgenumbers(currentcounter) = edgenumber;
                    irdiff{Rnumber(jj),Snumber(jj)}.locations(:,currentcounter) = xyzfirst(:);
                end                
            end
        end

    end
else
   irdiff = zeros(size(irgeom));
end

timingdata(3) = etime(clock,t00);
    
elapsedtimemakeirs = etime(clock,t00_total);

save(desiredname,'irdirect','irgeom','irdiff','timingdata','EDinputdatahash','elapsedtimemakeirs','stageversion','EDversionnumber');



