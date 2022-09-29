function EDmain_convex_time(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDmain_convex_time - Calculates the specular and edge diffraction IRs 
% for convex, rigid scattering objects, using the image source method for
% the specular reflection, and the explicit integral formulations for first-
% and higher-order diffraction.
%
% Input parameters are six structs with fields as specified below:
%   geoinputdata        .geoinputfile        (obligatory)
%                       As an alternative to specifying geoinputfile, it is
%                       possible to specify .corners, .planecorners.
%                       .firstcornertoskip   (default: 1e6)
%   Sinputdata          .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%                       .sourceamplitudes    (default:
%                        ones(nsources,1)
%   Rinputdata          .coordinates         (obligatory)
%   envdata             .cair                (default: 344)
%                       .rhoair              (default: 1.21)
%   controlparameters   .fs                  (default: 44100)
%                       .directsound         (default: 1 = yes)
%                       .difforder           (default: 2)
%                       .savealldifforders   (default: 0)
%                       .docalcir            (default: 1)
%                       .skipfirstorder      (default: 0)
%                       .Rstart              (default: 0)
%                       .discretizationtype  (default: 2 = G-L)
%                       .ngauss              (default: 16)
%   filehandlingparameters (optional)    
%                       .outputdirectory  (default: /result, in the folder of the geoinputfile)  
%                       .filestem        (default: name of the cad-file, with an underscore + running integer)
%                       .showtext            (default: 1)
%                       .suppressresultrecycling   (default: 0)
%                       .savecadgeofile      (default: 0)
%                       .saveSRdatafiles     (default: 1)
%                       .saveeddatafile      (default: 1)
%                       .saveed2datafile     (default: 1)
%                       .savepathsfile       (default: 1)
%                       .savelogfile         (default: 1)
%                       .savediff2result      (default: 0)
% 
% Uses the functions EDgetversion,EDcheckinputstructs, EDreadcad, EDreadgeomatrices,
% EDedgeo, EDSorRgeo, EDfindconvexGApaths, EDmakefirstorderirs, EDed2geo,
% EDinteg_submatrixstructure, EDintegralequation_convex_ir from EDtoolbox
% Uses the functions DataHash from Matlab Central
% 
% Peter Svensson 24 Sep 2022 (peter.svensson@ntnu.no)
%
% EDmain_convex_time(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 A first version = a trimmed version of ESIE2main,
% specialized for the scattering from a convex object.
% 28 Nov. 2017 Removed the JJ. Added the handling of timing data and
% writing to a logfile.
% 29 Nov. 2017 Added the semicolon after fclose(fid). Made some adjustments to
% handle input matrices instead of a cad file.
% 30 Nov. 2017 Fixed a missing part with calctfs and doaddsources = 1.
% Changed one field to savelogfile instead of logfilename. Removed the
% lineending field.
% 4 Dec. 2017 Added more info in the printouts and logfiles.
% 5 Dec. 2017 Added more info in the printouts and logfiles.
% 10 Dec. 2017 Added the timingstruct to the output file
% 13 Dec. 2017 Added the sourceamplitudes field to the Sinputdata struct
% 15 Dec. 2017 Added the number of non-zero matrix elements to the log
% file, and on-screen.
% 12 Jan. 2018 Made substantial changes to the GA part; calling new
% functions EDfindconvexGApaths and EDmakefirstordertfs
% 15 Jan. 2018 Started numbering the versions. Saved all the setup structs
% in the final result file.
% 15 Jan. 2018 Changed the parametername from EDversion to EDversionnumber
% because EDversion is a function.
% 16 Jan. 2018 Adapted to the name-change for the function EDgetversion
% 17 Jan 2018 Added the handling of the planecornertype input field.
% 17 Jan. 2018 Added the showtext input parameter to the "findGA..." and
% "makefirstordertfs".
% 18 Jan 2018 Changed input parameters to EDmakefirstordertfs. 
% 18 Jan 2018 Changed the order of the blocks, so that GA and diff1 comes
% before ESIE/HOD.
% 22 Jan 2018 Adapted to the removal of parameter planecornertype in
% EDreadcad and EDreadgeomatrices.
% 22 Jan 2018 Removed the input parameter nedgepoints_visibility. It is
% forced to the value 2 at the input to EDSorRgeo (but it should preferrably
% be avoided altogether for convex geometries.
% 22 Jan 2018 Corrected the handling of docalctf
% 23 Jan 2018 Split up the compstr == 'sun' || .... into two if-s. Hint
% from Jan Slechta.
% 23 Jan 2018 Added the savepathsfile handling; the parameter was already
% in place, but its handling was not implemented until now.
% 23 Jan 2018 Version 0.101 because of a bug in EDfindconvexGApaths
% 25 Jan 2018 Version 0.103: new version of EDwedge1st_fd.m
% 25 Jan 2018 Version 0.104: new version of ED_submatrixstructure, to
% enforce at least two edge points per edge (it could be one for very short
% edges before).
% 25 Jan 2018 Verson 0.105: new version of EDfindconvexGApaths. Previously,
% the direct sound was computed even if .directsound = 0.
% 26 Jan 2018 Version 0.106: fixed bug: when HOD was not calculated, an
% empty tfinteqdiff should have been saved, but that was not the case.
% Fixed now. Also, in EDwedge1st_fd, the case useserialexp2 had not been
% implemented. Fixed now.
% 26 Jan 2018 Change: "results" is not added to the outputdirectory in this
% function; it was added to the outputdirectory path as default.
% 26 Jan 2018 V 0.107: introduced the parameter doallSRcombinations, for
% backscatter cases.
% 31 Jan 2018 v0.108: Corrected the handling of sourceamplitudes
% (frequency dependence had not been implemented)
% 2 Feb 2018 v0.108: The timingstruct wasn't saved anywhere. Now it
% is saved in the tfinteq output file.
% 5 Feb. 2018 Moved the saving for the tfinteq file so that it could
% include the last timingdata.
% 8 Feb 2018 v0.109: Introduced a new parameter:
% controlparameters.skipfirstorder (default = 0).
% 8 Feb 2018 v0.110 First experiments w DataHash for recycling intermediate files
% 9 Feb 2018 v0.110 Added file duplication for the eddata,Sdata,Rdata files
% 13 Feb 2018 v0.111 Replaced quad_gauss with lgwt from Matlab Central.
% 14 Feb 2018 v0.112 Changed the file saving procedure again: now the tf
% and tfinteq files will always contain output data, instead of reference
% to another file. Also changed name for the input structs to geoinputdata,
% instead of geofiledata; Sinputdata instead of Sindata; Rinputdata instead
% of Rindata.
% 14 Feb 2018 Small change to the log file text, indicating that three
% edge points per wavelength is needed.
% 14 Feb 2018 Added an error message: if no S/R could see any planes, then
% the program stopped, with a warning sign.
% 15 Feb 2018 Removed the setupfile. Instead, introduced a master setting
% variable: EDsettings, a cell variable with all the 6 input structs +
% EDtoolbox version number. EDsettings will be saved in the tf and tfinteq 
% files. Also, stopped the automatic saving, with
% inputdatahash, if the corresponding savexxxxfile = 0. Introduced the
% parameter .suppressresultrecycling with default = 0, and implemented th
% corresponding suppressing of the result recycling.
% 28 Feb 2018 First version of EDmain_convexESIEtime
% 15 Mar 2018 First version of EDmain_convex_time
% 21 Mar 2018 Introduced the new parameter .savealldifforders
% 21 Mar 2018 Removed the parameter .hodtype
% 21 Mar 2018 Completed the section around EDmakeHODirs, with recycling
% result files.
% 14 May 2018 CLeaned up the lineending
% 21 May 2019 Adjusted to changes of the function EDfindHODpaths, with a
% new output parameter.
% 21 May 2019 Cleaned up some file recycling mistakes.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 20 Jan 2021 Fixed: the field .planerefltypes wasn't forwarded to
%             EDreadgeomatrices
% 25 Aug. 2021 Small improvement: introduced the parameter
% calcfirstorderdiff in the hash for EDmakefirstordertfs. Makes recycling
% possible a bit more often.
% 24 Sep. 2022 Changed the definition of the input data hash for the HODir.
% The previous version created a huge struct, and failed to identify
% existing results that could be reused.

[EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input data, assign default values if needed

if nargin < 6
    filehandlingparameters = struct('showtext',0);
    if nargin < 5
        error('ERROR: The first five input parameters to EDmain_convex_time must be specified')  
    end
end

% compstr = computer;
% compstr = lower(compstr(1:3));
% if compstr == 'mac'  
% 	lineending = 13;
% elseif compstr == 'sun' 
% 	lineending = 10;    
% elseif compstr == 'sol'            
% 	lineending = 10;
% elseif compstr == 'pcw'
% 	lineending = [13,10];
% else
%     error('ERROR: Not implemented for this computer type yet')	
% end
lineending = 10;
if ispc == 1
   lineending = [13,10];
end

[geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,4);

EDsettings = cell(7,1);
EDsettings{1} = geoinputdata;
EDsettings{2} = Sinputdata;
EDsettings{3} = Rinputdata;
EDsettings{4} = envdata;
EDsettings{5} = controlparameters;
EDsettings{6} = filehandlingparameters;
EDsettings{7} = EDversionnumber;

if filehandlingparameters.savelogfile == 1
    logfilename = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_log.txt'];
end

if filehandlingparameters.showtext >= 1
	disp('    ');disp('####################################################################')
              disp(['#  EDmain_convex_time, v. ',num2str(EDversionnumber),' (',lastsavedate,')'])
              disp(['#  filestem for results: ',filehandlingparameters.filestem])
              disp(' ')
end
if filehandlingparameters.savelogfile == 1
    fid = fopen(logfilename,'w');
    if fid == -1
        disp('This file is not possible to open - check that it isn''t opened by any program!')
    	return
    end
    fwrite(fid,['####################################################################',lineending],'char');
    fwrite(fid,['#  EDmain_convex_time, v. ',num2str(EDversionnumber),' (',lastsavedate,')',lineending],'char');
    fwrite(fid,['#  filestem for results: ',filehandlingparameters.filestem,lineending],'char');
    fwrite(fid,[' ',lineending],'char');
end

% nfrequencies = length(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

if isfield(geoinputdata,'geoinputfile')
    if filehandlingparameters.showtext >= 1
        disp('   Creating the planedata struct from the CAD file')
    end
    t00 = clock;
    [planedata,extraCATTdata] = EDreadcad(geoinputdata.geoinputfile,0);
    if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
        error('ERROR: EDmain_convex_time can only be used for convex scatterers, including a single thin plate')
    end
    if filehandlingparameters.savecadgeofile == 1
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_cadgeo.mat'];
        eval(['save ',desiredname,' planedata extraCATTdata'])
    end
    t01 = etime(clock,t00);
    ncorners = size(planedata.corners,1);
    nplanes = (size(planedata.planecorners,1));
    if filehandlingparameters.showtext >= 1
        disp(['      ',int2str(ncorners),' corners and ',int2str(nplanes),' planes'])
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDreadcad (',int2str(ncorners),' corners and ',int2str(nplanes),' planes), time: ',num2str(t01),' s',lineending],'char');
    end
else
    if filehandlingparameters.showtext >= 1
        disp('   Creating the planedata struct from the input geometry matrices')
    end
    t00 = clock;
    % Fix on 20 Jan 2021: the field .planerefltypes wasn't forwarded to
    % EDreadgeomatrices
    planedata = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners,geoinputdata.planerefltypes);    
    ncorners = size(planedata.corners,1);
    nplanes = size(planedata.planecorners,1);
    t01 = etime(clock,t00);
    if filehandlingparameters.showtext >= 1
        disp(['      ',int2str(ncorners),' corners and ',int2str(nplanes),' planes'])
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDreadgeomatrices (',int2str(ncorners),' corners and ',int2str(nplanes),' planes), time: ',num2str(t01),' s',lineending],'char');
    end    
end
timingstruct = struct('geoinput',t01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the planedata struct and create an edgedata struct

if filehandlingparameters.showtext >= 1	
	disp('   Creating the edgedata struct ')
end

t00 = clock;
if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
else
    EDedgeoinputstruct = struct('planedata',planedata,'firstcornertoskip',geoinputdata.firstcornertoskip,...
        'listofcornerstoskip',[],'planeseesplanestrategy',0,'EDversionnumber',EDversionnumber);
    EDedgeoinputhash = DataHash(EDedgeoinputstruct);
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_eddata',EDedgeoinputhash);
end
desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_eddata.mat'];

if foundmatch == 1
    eval(['load ''',existingfilename,''''])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
else
    [edgedata,planedata,EDinputdatahash] = EDedgeo(planedata,EDversionnumber,geoinputdata.firstcornertoskip,[],0,filehandlingparameters.showtext);
    if filehandlingparameters.saveeddatafile == 1
        eval(['save ''',desiredname,''' planedata edgedata EDinputdatahash'])
    end
end
nedges = size(edgedata.edgecorners,1);
t01 = etime(clock,t00);
timingstruct.edgedata = t01;
if filehandlingparameters.showtext >= 1
    if foundmatch == 1
        disp(['      Recycled and duplicated ',existingfilename])
    end
     disp(['      ',int2str(nedges),' edges'])
end
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDedgeo, (',int2str(nedges),' edges), time: ',num2str(t01),' s',lineending],'char');
    if foundmatch == 1
        fwrite(fid,['      by recycling and duplicating ',existingfilename,lineending],'char');        
    end
end
if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
    error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Sdata and Rdata structs

if filehandlingparameters.showtext >= 1	
	disp('   Creating the Sdata struct ')
end
t00 = clock;
if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
else
    EDSorRgeoinputstruct = struct('planedata',planedata,'edgedata',edgedata, 'pointcoords',Sinputdata.coordinates,...
        'nedgesubs',2,'EDversionnumber',EDversionnumber);
    EDSorRgeoinputhash = DataHash(EDSorRgeoinputstruct);
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_Sdata',EDSorRgeoinputhash);
end
desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_Sdata.mat'];
if foundmatch == 1
    eval(['load ''',existingfilename,''''])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
else
    [Sdata,EDinputdatahash] = EDSorRgeo(planedata,edgedata,Sinputdata.coordinates,'S',EDversionnumber,2,filehandlingparameters.showtext);
    if filehandlingparameters.saveSRdatafiles == 1
        eval(['save ''',desiredname,''' Sdata EDinputdatahash'])
    end    
end
nsources = size(Sdata.sources,1);
t01 = etime(clock,t00);
timingstruct.Sdata = t01;
if any(any(Sdata.visplanesfroms)) == 0
   error('ERROR: No source can see any of the planes of the scattering object. Please check your geometry.') 
end

if filehandlingparameters.showtext >= 1
    if foundmatch == 1
        disp(['      Recycled and duplicated ',existingfilename])
    end
     disp(['      ',int2str(nsources),' source(s)'])
end
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDSorRgeo(S), (',int2str(nsources),' source(s)), time: ',num2str(t01),' s',lineending],'char');
    if foundmatch == 1
        fwrite(fid,['      by recycling and duplicating ',existingfilename,lineending],'char');        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if filehandlingparameters.showtext >= 1
	disp('   Creating the Rdata struct ')
end
t00 = clock;
if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
else
    EDSorRgeoinputstruct = struct('planedata',planedata,'edgedata',edgedata, 'pointcoords',Rinputdata.coordinates,...
        'nedgesubs',2,'EDversionnumber',EDversionnumber);
    EDSorRgeoinputhash = DataHash(EDSorRgeoinputstruct);

    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_Rdata',EDSorRgeoinputhash);
end
desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_Rdata.mat'];
if foundmatch == 1
    eval(['load ''',existingfilename,''''])
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
else
    [Rdata,EDinputdatahash] = EDSorRgeo(planedata,edgedata,Rinputdata.coordinates,'R',EDversionnumber,2,filehandlingparameters.showtext);
    if filehandlingparameters.saveSRdatafiles == 1
        eval(['save ''',desiredname,''' Rdata EDinputdatahash'])
    end
end
nreceivers = size(Rdata.receivers,1);
t01 = etime(clock,t00);
timingstruct.Rdata = t01;
if any(any(Rdata.visplanesfromr)) == 0
   error('ERROR: No receiver can see any of the planes of the scattering object. Please check your geometry.') 
end

if filehandlingparameters.showtext >= 1
    if foundmatch == 1
        disp(['      Recycled and duplicated ',existingfilename])
    end
     disp(['      ',int2str(nreceivers),' receiver(s)'])
end
if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   EDSorRgeo(R), (',int2str(nreceivers),' receiver(s)), time: ',num2str(t01),' s',lineending],'char');
    if foundmatch == 1
        fwrite(fid,['      by recycling and duplicating ',existingfilename,lineending],'char');        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the paths for direct sound, first-order specular, and first-order
% diffraction paths.

if controlparameters.skipfirstorder == 0
    if filehandlingparameters.showtext >= 1	
        disp('   Find the (first-order) GA paths')
    end

    t00 = clock;
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else    
        EDfindconvGApathsinputstruct = struct('planedata',planedata,'edgedata',edgedata,...
            'sources',Sdata.sources,'visplanesfromS',Sdata.visplanesfroms,'vispartedgesfromS',Sdata.vispartedgesfroms,...
            'receivers',Rdata.receivers,'visplanesfromR',Rdata.visplanesfromr,'vispartedgesfromR',Rdata.vispartedgesfromr,...
            'difforder',controlparameters.difforder,'directsound',controlparameters.directsound,...
            'doallSRcombinations',Sinputdata.doallSRcombinations,'EDversionnumber',EDversionnumber);
        EDfindconvGApathsinputhash = DataHash(EDfindconvGApathsinputstruct);
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_paths',EDfindconvGApathsinputhash);
    end
    if foundmatch == 1
        eval(['load ''',existingfilename,''''])
    else
        [firstorderpathdata,EDinputdatahash] = EDfindconvexGApaths(planedata,edgedata,...
            Sdata.sources,Sdata.visplanesfroms,Sdata.vispartedgesfroms,...
            Rdata.receivers,Rdata.visplanesfromr,Rdata.vispartedgesfromr,...
            controlparameters.difforder,controlparameters.directsound,Sinputdata.doallSRcombinations,...
            EDversionnumber,filehandlingparameters.showtext);
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_paths.mat'];
        if filehandlingparameters.savepathsfile == 1
            eval(['save ''',desiredname,''' firstorderpathdata EDinputdatahash'])   
        end
    end
    t01 = etime(clock,t00);
    timingstruct.findpaths = t01;
    if filehandlingparameters.showtext >= 1	&& foundmatch == 1
        disp(['      Recycled ',existingfilename])
    end    
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDfindconvexGApaths',lineending],'char');
        fwrite(fid,['                                Total time: ',num2str(t01),' s',lineending],'char');
        if foundmatch == 1
            fwrite(fid,['      by recycling ',existingfilename,lineending],'char');        
        end
    end
else
    if filehandlingparameters.showtext >= 1	
        disp('   The first-order GA and diff paths are not identified because skipfirstorder was set to 1')
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDfindconvexGApaths was not run because skipfirstorder was set to 1',lineending],'char');
    end
    timingstruct.findpaths = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the first-order specular, and first-order
% diffraction irs.

if controlparameters.docalcir == 1 && controlparameters.skipfirstorder == 0
    if filehandlingparameters.showtext >= 1	
        disp('   Generate the (first-order) GA and diff irs.')
    end

    t00 = clock;

    % New parameter for the hash 25 Aug. 2021: calcfirstorderdiff
    calcfirstorderdiff = double(controlparameters.difforder > 0);

    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
        EDfirstorderirsinputstruct = struct('firstorderpathdata',firstorderpathdata,'edgedata',edgedata,...
            'fs',controlparameters.fs,'Rstart',controlparameters.Rstart,...
            'calcfirstorderdiff',calcfirstorderdiff,'envdata',envdata,'Sinputdata',Sinputdata,...
            'receivers',Rdata.receivers,'saveindividualfirstdiff',controlparameters.saveindividualfirstdiff,'EDversionnumber',EDversionnumber);
        EDfirstorderirsinputhash = DataHash(EDfirstorderirsinputstruct);
        [foundmatch,recycledresultsfile] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_ir',EDfirstorderirsinputhash);
    end
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_ir.mat'];    
    if foundmatch == 1
        currenttimingstruct = timingstruct;
         eval(['load ''',recycledresultsfile,''' irdirect irgeom irdiff timingstruct'])
         EDinputdatahash = EDfirstorderirsinputhash;
        t01 = etime(clock,t00);
        currenttimingstruct.makeirs = [t01 timingstruct.makeirs(2:4)];
        timingstruct = currenttimingstruct;
        eval(['save ''',desiredname,''' irdirect irgeom irdiff timingstruct recycledresultsfile EDsettings EDinputdatahash'])
    else
        [irdirect,irgeom,irdiff,timingdata,EDinputdatahash] = EDmakefirstorderirs(firstorderpathdata,...
            controlparameters.fs,controlparameters.Rstart,controlparameters.difforder,envdata,Sinputdata,Rdata.receivers,...
            edgedata,controlparameters.saveindividualfirstdiff,EDversionnumber,filehandlingparameters.showtext);        
%         [tfdirect,tfgeom,tfdiff,timingdata,EDinputdatahash] = EDmakefirstordertfs(firstorderpathdata,...
%             controlparameters.frequencies,controlparameters.Rstart,controlparameters.difforder,envdata,Sinputdata,Rdata.receivers,...
%             edgedata,EDversionnumber,filehandlingparameters.showtext);
        recycledresultsfile = '';
        t01 = etime(clock,t00);
        timingstruct.makeirs = [t01 timingdata];
        eval(['save ''',desiredname,''' irdirect irgeom irdiff timingstruct recycledresultsfile EDsettings EDinputdatahash'])
    end
    if filehandlingparameters.showtext >= 1 && foundmatch == 1
        disp(['      Recycled ',recycledresultsfile])
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDmakefirstorderirs',lineending],'char');
        if foundmatch == 1
            fwrite(fid,['      by recycling ',recycledresultsfile,lineending],'char');     
            fwrite(fid,['                                Total time: ',num2str(t01),' s',lineending],'char');          
        else
            fwrite(fid,['                                Total time: ',num2str(t01),' s. Parts, as below',lineending],'char');
            fwrite(fid,['                                Generate the direct sound: ',num2str(timingdata(1)),' s',lineending],'char');
            fwrite(fid,['                                Generate the specular reflections: ',num2str(timingdata(2)),' s',lineending],'char');
            fwrite(fid,['                                Generate the first-order diffraction: ',num2str(timingdata(3)),' s',lineending],'char');
        end
    end
else
    if filehandlingparameters.showtext >= 1	
        disp('   First-order GA and diff irs are not generated because docalcir was set to 0, or skipfirstorder was set to 1')
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDmakefirstorderirs was not run because docalcir was set to 0, or skipfirstorder was set to 1',lineending],'char');
    end
    timingstruct.makeirs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add edge-to-edge visibility data

if controlparameters.difforder > 1
    if filehandlingparameters.showtext >= 1	
        disp('   Creating the edgetoedgedata struct ')
    end
    t00 = clock;

    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
        EDed2geoinputstruct = struct('planedata',planedata,'edgedata',edgedata,...
            'Sdata',Sdata,'Rdata',Rdata,'specorder',1,'difforder',2,...
            'nedgesubs',2,'ndiff2batches',1,'EDversionnumber',EDversionnumber);
        EDed2geoinputhash = DataHash(EDed2geoinputstruct);
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_ed2data',EDed2geoinputhash);
    end    
    if foundmatch == 1
        eval(['load ''',existingfilename,''''])
    else 
        [edgetoedgedata,EDinputdatahash] = EDed2geo(edgedata,planedata,Sdata,Rdata,1,2,EDversionnumber,2,1,filehandlingparameters.showtext);    
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_ed2data.mat'];
        if filehandlingparameters.saveed2datafile == 1
            eval(['save ''',desiredname,''' planedata edgedata edgetoedgedata EDinputdatahash'])
        end
    end
    t01 = etime(clock,t00);
    timingstruct.edgetoedgedata = t01;
    if filehandlingparameters.showtext >= 1	&& foundmatch == 1
        disp(['      Recycled ',existingfilename])
    end    
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDed2geo, time: ',num2str(t01),' s',lineending],'char');
        if foundmatch == 1
                fwrite(fid,['      by recycling ',existingfilename,lineending],'char');        
        end
    end
else
    if filehandlingparameters.showtext >= 1	
        disp(['   Skipping the edgetoedgedata struct, since difforder = ',int2str(controlparameters.difforder)])
    end
    timingstruct.edgetoedgedata = 0;
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDed2geo was not run',lineending],'char');
    end    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the higher-order diffraction paths

if controlparameters.difforder > 1 && controlparameters.docalcir == 1
     
    if filehandlingparameters.showtext >= 1	
        disp('   Finding the higher-order diffraction paths ')
    end
    t00 = clock;
    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
%         EDhodpathsinputstruct = struct();
%         EDhodpathsinputhash = DataHash(EDhodpathsinputstruct);
        EDhodpathsinputstruct = struct('difforder',controlparameters.difforder,...
        'edgeseespartialedge',int8(sign(edgetoedgedata.edgeseespartialedge)),...
        'vispartedgesfroms',sign(Sdata.vispartedgesfroms),'vispartedgesfromr',sign(Rdata.vispartedgesfromr),'EDversionnumber',EDversionnumber);
        EDhodpathsinputhash = DataHash(EDhodpathsinputstruct);
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_hodpaths',EDhodpathsinputhash);
    end
    if foundmatch == 1
        eval(['load ''',existingfilename,''''])
    else
        [hodpaths,hodpathsalongplane,EDinputdatahash] = EDfindHODpaths(int8(sign(edgetoedgedata.edgeseespartialedge)),sign(Sdata.vispartedgesfroms),...
            sign(Rdata.vispartedgesfromr),controlparameters.difforder,EDversionnumber);
%         [hodpaths,EDinputdatahash] = EDfindHODpaths(uint8(sign(edgetoedgedata.edgeseespartialedge)),sign(Sdata.vispartedgesfroms),...
%             sign(Rdata.vispartedgesfromr),controlparameters.difforder,EDversionnumber);
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_hodpaths.mat'];
        if filehandlingparameters.savehodpaths == 1
            eval(['save ''',desiredname,''' hodpaths hodpathsalongplane EDinputdatahash'])
        end
    end
    t01 = etime(clock,t00);
    timingstruct.hodpaths = t01;
    if filehandlingparameters.showtext >= 1    
         if foundmatch == 1
            disp(['      Recycled ',existingfilename])
         end
         disp(['      Found higher-order diffraction paths up to order ',int2str(controlparameters.difforder)]) 
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDfindHODpaths, diffraction order ',int2str(controlparameters.difforder),', time: ',num2str(t01),' s',lineending],'char');
        if foundmatch == 1
                fwrite(fid,['      by recycling ',existingfilename,lineending],'char');        
        end
    end
else
    if filehandlingparameters.showtext >= 1	
        disp(['   Skipping EDfindHODpaths, since difforder = ',int2str(controlparameters.difforder),' (or docalcir was set to 0)'])
        timingstruct.hodpaths = 0;
    end    
    if filehandlingparameters.savelogfile >= 1	        
        fwrite(fid,['   EDfindHODpaths was not run',lineending],'char');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contributions one order at a time

if controlparameters.difforder > 1 && controlparameters.docalcir == 1
     
    if filehandlingparameters.showtext >= 1	
        disp('   Generating the higher-order diffraction irs ')
    end
    t00 = clock;
    elemsize = 2.^(-[0:controlparameters.difforder-1]);
    elemsize = elemsize*2;
    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
%         EDhodirsinputstruct = struct();
%         EDhodirinputhash = DataHash(EDhodirsinputstruct);
        EDhodirinputstruct = struct('difforder',controlparameters.difforder);       
        EDhodirinputstruct.hodpaths = hodpaths;
        EDhodirinputstruct.hodpathsalongplane = hodpathsalongplane;
        EDhodirinputstruct.elemsize = elemsize;
        EDhodirinputstruct.edgedata = edgedata;
        EDhodirinputstruct.Sdata = Sdata;
        EDhodirinputstruct.doaddsources = Sinputdata.doaddsources;
        EDhodirinputstruct.sourceamplitudes = Sinputdata.sourceamplitudes;
        EDhodirinputstruct.Rdata = Rdata;
        EDhodirinputstruct.cair = envdata.cair;
        EDhodirinputstruct.fs = controlparameters.fs;
        EDhodirinputstruct.Rstart = controlparameters.Rstart;
        EDhodirinputstruct.savealldifforders = controlparameters.savealldifforders;
        EDhodirinputstruct.EDversionnumber = EDversionnumber;
        EDhodirinputhash = DataHash(EDhodirinputstruct);
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_irhod',EDhodirinputhash);
    end
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_irhod.mat'];
    if foundmatch == 1
        currenttimingstruct = timingstruct;
         eval(['load ''',existingfilename,''' irhod'])
         EDinputdatahash = EDhodirinputhash;
        t01 = etime(clock,t00);
        currenttimingstruct.hodir = t01;
        timingstruct = currenttimingstruct;
        eval(['save ''',desiredname,'''  irhod timingstruct EDsettings EDinputdatahash'])        
    else        
        [irhod,EDinputdatahash] = EDmakeHODirs(hodpaths,hodpathsalongplane,controlparameters.difforder,elemsize,edgedata,...
        edgetoedgedata,Sdata,Sinputdata.doaddsources,Sinputdata.sourceamplitudes,Rdata,envdata.cair,...
        controlparameters.fs,controlparameters.Rstart,controlparameters.savealldifforders,...
        filehandlingparameters.showtext,EDversionnumber);

        existingfilename =  '';
        t01 = etime(clock,t00);
        timingstruct.hodir = t01;
        eval(['save ''',desiredname,''' irhod timingstruct EDsettings EDinputdatahash'])
    end
    if filehandlingparameters.showtext >= 1    
         if foundmatch == 1
            disp(['      Recycled ',existingfilename])
         end
         disp(['      Generated higher-order diffraction irs up to order ',int2str(controlparameters.difforder)]) 
    end
    if filehandlingparameters.savelogfile == 1
        fwrite(fid,['   EDmakeHODirs, diffraction order ',int2str(controlparameters.difforder),', time: ',num2str(t01),' s',lineending],'char');
        if foundmatch == 1
                fwrite(fid,['      by recycling ',existingfilename,lineending],'char');        
        end
    end
else
    if filehandlingparameters.showtext >= 1	
        disp(['   Skipping EDmakeHODirs, since difforder = ',int2str(controlparameters.difforder),' (or docalcir was set to 0)'])
        timingstruct.hodpaths = 0;
    end    
    if filehandlingparameters.savelogfile >= 1	        
        fwrite(fid,['   EDmakeHODirs was not run',lineending],'char');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if filehandlingparameters.savelogfile == 1
    fclose(fid);
end


