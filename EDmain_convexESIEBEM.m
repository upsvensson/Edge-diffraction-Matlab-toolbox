function EDmain_convexESIEBEM(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDmain_convexESIEBEM - Calculates the total scattered field, in the frequency domain,
% for convex, rigid scattering objects, using the ESIEBEM approach. As a first
% step, intermediate receiver positions are distrbuted along the surfaces of
% the object, and the ordinary ESIE is used to calculate the sound field in
%  those positions. In a second step, this surface sound field is
%  propagated to the original external receiver positions using the
%  Helmholtz integral.
%
% Input parameters are six structs with fields as specified below:
%   geoinputdata        .geoinputfile        (obligatory)
%                       As an alternative to specifying geoinputfile, it is
%                       possible to specify .corners, .planecorners.
%                       .firstcornertoskip   (default: 1e6)
%   Sinputdata          .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%                       .sourceamplitudes    (default:
%                        ones(nsources,nfrequencies)
%   Rinputdata          .coordinates         (obligatory)
%   envdata             .cair                (default: 344)
%                       .rhoair              (default: 1.21)
%   controlparameters   .fs                  (default: 44100)
%                       .directsound         (default: 1 = yes)
%                       .difforder           (default: 15)
%                       .docalctf            (default: 1)
%                       .docalcir            (default: 0)
%                       .skipfirstorder      (default: 0)
%                       .Rstart              (default: 0)
%                       .frequencies         (obligatory, if docalcftf = 1)
%                       .discretizationtype  (default: 2 = G-L)
%                       .ngauss              (default: 16)
%                       .surfacegaussorder   (default: 5)
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
% EDmessage, EDpostfunctext,EDedgeo, EDSorRgeo, EDfindconvexGApaths, EDmakefirstordertfs, EDed2geo,
% EDinteg_submatrixstructure, EDintegralequation_convex_tf from EDtoolbox
% Uses the functions DataHash from Matlab Central
% 
% Peter Svensson 27 Sep. 2023 (peter.svensson@ntnu.no)
%
% EDmain_convexESIEBEM(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

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
% 1 Mar 2018 First attempt at the ESIEBEM version
% 2 Mar 2018 Fixed two things, which gave expected results. First: the
% free-field component must be added to the surface-propagated pressure.
% Second: if Rstart must be taken into account in the free-field component.
% 2 Mar 2018 Completed the info printing, for the ESIEBEM part.
% 7 Mar 2018 Completed the handling of different numbers of sources,
% receivers, frequencies in the last ESIEBEM section.
% 4 Apr 2018 Removed the showtext parameter from the EDgensurfreceivers
% 14 May 2018 CLeaned up the lineending business.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 20 Jan 2021 Fixed: the field .planerefltypes wasn't forwarded to
%             EDreadgeomatrices
% 25 Aug. 2021 Small improvement: introduced the parameter
% calcfirstorderdiff in the hash for EDmakefirstordertfs. Makes recycling
% possible a bit more often.
% 27 Sep. 2023 Introduced the functions EDmessage and EDpostfunctext to
% make this main function more readable. 

[EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input data, assign default values if needed

if nargin < 6
    filehandlingparameters = struct('showtext',0);
    if nargin < 5
        error('ERROR: The first five input parameters to EDmain_convexESIEBEM must be specified')  
    end
end

if ispc == 1
   lineending = [13,10];
else
    lineending = 10;
end

[geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,3);

textline1 = ['EDmain_convexESIEBEM, v.',num2str(EDversionnumber),' (',lastsavedate,')'];
textline2 = ['filestem for results: ',filehandlingparameters.filestem];
textline3 = ['directory for results: ',filehandlingparameters.outputdirectory];
fid = EDmessage(filehandlingparameters,'sf',0,1,'',textline1,textline2,textline3);

nfrequencies = length(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

t00 = clock;
if isfield(geoinputdata,'geoinputfile')    
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the CAD file']);
    [planedata,extraCATTdata] = EDreadcad(geoinputdata.geoinputfile,0);
    textline1start = ['   EDreadcad: '];
    t01 = etime(clock,t00);
else
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the input geometry matrices']);
    planedata = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners,geoinputdata.planerefltypes); 
    textline1start = ['   EDreadgeomatrices: '];
    t01 = etime(clock,t00);
end
if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
     error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
end
ncorners = size(planedata.corners,1);
nplanes = (size(planedata.planecorners,1));
Timestring = [' Time: ',num2str(t01(1)),' s'];
if length(t01) > 1, Timestring = [Timestring,' (Orig.: ',num2str(t01(2)), 's)'];  end
textline1 = [textline1start,int2str(ncorners),' corners and ',int2str(nplanes),' planes.',Timestring];
    
EDmessage(filehandlingparameters,'sf',fid,1,'',textline1);
timingstruct = struct('geoinput',t01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store the original receiver points as "fieldpoints", and generate new
% receiver points at the surfaces of the polyhedron.

EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the new surface receiver positions']);
t00 = clock;
fieldpoints = Rinputdata.coordinates;
doallSRcombinations_fieldpoints = Sinputdata.doallSRcombinations;
directsound_fieldpoints = controlparameters.directsound;

% [surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
%     controlparameters.surfacegaussorder,EDversionnumber,filehandlingparameters.showtext);
[surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
    controlparameters.surfacegaussorder,EDversionnumber);
Rinputdata.coordinates = surfacerecs;
Rinputdata.nvecs = surfacerecnvecs;
Rinputdata.weights = surfacerecweights;
Sinputdata.doallSRcombinations = 1;
controlparameters.directsound = 1;
t01 = etime(clock,t00);
timingstruct.surfrecs = t01;

EDsettings = cell(7,1);
EDsettings{1} = geoinputdata;
EDsettings{2} = Sinputdata;
EDsettings{3} = Rinputdata;
EDsettings{4} = envdata;
EDsettings{5} = controlparameters;
EDsettings{6} = filehandlingparameters;
EDsettings{7} = EDversionnumber;

textline1 = ['       Created ',int2str(size(Rinputdata.coordinates,1)),' surface receiver positions'];
EDpostfunctext('EDgensurfreceivers',t01,existingfilename,...
    filehandlingparameters,fid,'',textline1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the planedata struct and create an edgedata struct

EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the edgedata struct']);
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
EDpostfunctext('EDedgeo',t01,existingfilename,...
    filehandlingparameters,fid,'',[int2str(nedges),' edges.'])

if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
    error('ERROR: EDmain_convexESIE can only be used for convex scatterers, including a single thin plate')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Sdata and Rdata structs

EDmessage(filehandlingparameters,'s',fid,1,'',[TTT,'Creating the Sdata struct']);
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
EDpostfunctext('EDSorRgeo',t01,existingfilename,...
    filehandlingparameters,fid,'',[int2str(nsources),' sources.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Rdata struct

EDmessage(filehandlingparameters,'s',fid,1,'',[TTT,'Creating the Rdata struct']);
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
 EDpostfunctext('EDSorRgeo',t01,existingfilename,...
    filehandlingparameters,fid,'',[int2str(nreceivers),' receivers.'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the paths for direct sound, first-order specular, and first-order
% diffraction paths.

if controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Find the (first-order) GA paths']);
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
    nspecs = size(firstorderpathdata.specrefllist,1);
    ndiffs = prod(size(firstorderpathdata.diffpaths));
    t01 = etime(clock,t00);
    timingstruct.findpaths = t01;
    EDpostfunctext('EDfindconvexGApaths',t01,existingfilename,...
        filehandlingparameters,fid,'',[int2str(nspecs),' spec. paths and ',...
        int2str(ndiffs),' diff. paths.'])

else
    textline1 = ['   EDfindconvexGApaths was not run because skipfirstorder was set to 1'];
    EDmessage(filehandlingparameters,'fs',fid,1,'',textline1);
    timingstruct.findpaths = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the first-order specular, and first-order
% diffraction tfs.

if controlparameters.docalctf == 1 && controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Generate the (first-order) GA and diff tfs']);
    t00 = clock;

    % New parameter for the hash 25 Aug. 2021: calcfirstorderdiff
    calcfirstorderdiff = double(controlparameters.difforder > 0);
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
        EDfirstordertfsinputstruct = struct('firstorderpathdata',firstorderpathdata,'edgedata',edgedata,...
            'frequencies',controlparameters.frequencies,'Rstart',controlparameters.Rstart,...
            'calcfirstorderdiff',calcfirstorderdiff,'envdata',envdata,'Sinputdata',Sinputdata,...
            'receivers',Rdata.receivers,'EDversionnumber',EDversionnumber);
        EDfirstordertfsinputhash = DataHash(EDfirstordertfsinputstruct);
        [foundmatch,recycledresultsfile] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_tf',EDfirstordertfsinputhash);
    end
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'];    
    if foundmatch == 1
        currenttimingstruct = timingstruct;
         eval(['load ''',recycledresultsfile,''' tfdirect tfgeom tfdiff timingstruct'])
         EDinputdatahash = EDfirstordertfsinputhash;
        t01 = etime(clock,t00);
        currenttimingstruct.maketfs = [t01 timingstruct.maketfs(2:4)];
        timingdata = timingstruct.maketfs(2:4);
        eval(['save ''',desiredname,''' tfdirect tfgeom tfdiff timingstruct recycledresultsfile EDsettings EDinputdatahash'])
    else
        [tfdirect,tfgeom,tfdiff,timingdata,EDinputdatahash] = EDmakefirstordertfs(firstorderpathdata,...
            controlparameters.frequencies,controlparameters.Rstart,controlparameters.difforder,envdata,Sinputdata,Rdata.receivers,...
        edgedata,EDversionnumber,filehandlingparameters.showtext);
        recycledresultsfile = '';
        t01 = etime(clock,t00);
        timingstruct.maketfs = [t01 timingdata];
        eval(['save ''',desiredname,''' tfdirect tfgeom tfdiff timingstruct recycledresultsfile EDsettings EDinputdatahash'])
    end
    textline1 = [int2str(nfrequencies),' frequencies.'];
    textline2 = ['      Direct sound: ',num2str(timingdata(1)),'s, spec.refl.: ',...
        num2str(timingdata(2)),'s, first-order diff.: ',num2str(timingdata(3)),'s'];
    EDpostfunctext('EDmakefirstordertfs',t01,existingfilename,filehandlingparameters,fid,'',textline1,textline2);    
else
    textline1 = ['   EDmakefirstordertfs was not run because docalctf was set to 0, or skipfirstorder was set to 1'];
    EDmessage(filehandlingparameters,'fs',fid,1,'',textline1);
    timingstruct.maketfs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add edge-to-edge visibility data

DIF = int2str(controlparameters.difforder);
if controlparameters.difforder > 1
   EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the edgetoedgedata struct']);
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
    EDpostfunctext('EDed2geo',t01,existingfilename,...
        filehandlingparameters,fid,'');
else
    textline1 = ['   EDed2geo was not run since difforder = ',DIF];
    EDmessage(filehandlingparameters,'sf',fid,1,'',textline1);
    timingstruct.edgetoedgedata = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for the integral equation: set up the submatrix structure

if controlparameters.difforder > 1 && controlparameters.docalctf == 1
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the integral equation Hsubmatrixdata struct']);
    t00 = clock;    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
        EDsubmatrixinputstruct = struct();
        EDsubmatrixinputhash = DataHash(EDsubmatrixinputstruct);
        EDsubmatrixinputstruct = struct('edgelengthvec',edgedata.edgelengthvec,'closwedangvec',edgedata.closwedangvec,...
        'inteq_ngauss',controlparameters.ngauss,'inteq_discretizationtype',controlparameters.discretizationtype,...
        'edgetoedgedata',edgetoedgedata,'planesatedge',edgedata.planesatedge,'EDversionnumber',EDversionnumber);
        EDsubmatrixinputhash = DataHash(EDsubmatrixinputstruct);
        [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_submatrixdata',EDsubmatrixinputhash);
    end
    if foundmatch == 1
        eval(['load ''',existingfilename,''''])
    else
        [Hsubmatrixdata,EDinputdatahash] = EDinteg_submatrixstructure(edgedata.edgelengthvec,edgedata.closwedangvec,...
        controlparameters.ngauss,controlparameters.discretizationtype,edgetoedgedata,edgedata.planesatedge,EDversionnumber,filehandlingparameters.showtext);
    
        desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_submatrixdata.mat'];
        if filehandlingparameters.savesubmatrixdata == 1
            eval(['save ''',desiredname,''' Hsubmatrixdata EDinputdatahash'])
        end
    end
%    nsousigs = Hsubmatrixdata.bigmatrixendnums(end);
%    nsubmatrices = size(Hsubmatrixdata.edgetripletlist,1);
    edgeelemsizes = edgedata.edgelengthvec./Hsubmatrixdata.nedgeelems;
    meanelemsize = mean(edgeelemsizes);
%    maxfreq = envdata.cair/(3*meanelemsize);
%    minedgeelemnumber = min(Hsubmatrixdata.nedgeelems);
%    maxedgeelemnumber = max(Hsubmatrixdata.nedgeelems);
    nonzeroelements = sum(prod(Hsubmatrixdata.nedgeelems(Hsubmatrixdata.edgetripletlist),2));
    t01 = etime(clock,t00);
    timingstruct.submatrixdata = t01;
    EDpostfunctext('EDinteg_submatrixstructure',t01,existingfilename,...
        filehandlingparameters,fid,'',...
        [int2str(Hsubmatrixdata.nuniquesubmatrices),' submatrices, out of ',int2str(size(Hsubmatrixdata.edgetripletlist,1)),', to compute.'],...
        ['      Edges discretized with: ',int2str(min(Hsubmatrixdata.nedgeelems)),'-',...
        int2str(max(Hsubmatrixdata.nedgeelems)),' gauss points. ',int2str(Hsubmatrixdata.bigmatrixendnums(end)),' edge source signals to compute.'],...
        ['      The IE matrix has ',int2str(nonzeroelements),' non-zero elements. Some may be identical due to symmetries'])    
    
 else
   textline1 = ['   EDinteg_submatrixstructure was not run because difforder = ',DIF,' or docalctf was set to 0'];
    EDmessage(filehandlingparameters,'sf',fid,1,'',textline1);
    timingstruct.submatrixdata = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contribution with the integral equation

if controlparameters.difforder > 1 && controlparameters.docalctf == 1
    textline2 = ['      ',int2str(length(controlparameters.frequencies)),' frequencies. Diffraction order: ',DIF];
    EDmessage(filehandlingparameters,'s',fid,1,'',['Calculating the HOD contribution with the FD integral equation'],textline2);

    t00 = clock;    
    if filehandlingparameters.suppressresultrecycling == 1
        foundmatch = 0;
    else
        EDinteqinputstruct = struct('envdata',envdata,...
            'planedata',planedata,'edgedata',edgedata,'edgetoedgedata',edgetoedgedata,...
            'Hsubmatrixdata',Hsubmatrixdata,'Sdata',Sdata,'doaddsources',Sinputdata.doaddsources,...
            'sourceamplitudes',Sinputdata.sourceamplitudes,'doallSRcombinations',Sinputdata.doallSRcombinations,...
            'Rdata',Rdata,'controlparameters',controlparameters,'EDversionnumber',EDversionnumber);
        EDinteqinputhash = DataHash(EDinteqinputstruct);
        [foundmatch,recycledresultsfile] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_tfinteq',EDinteqinputhash);
    end    
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'];
    
    if foundmatch == 1
        currenttimingstruct = timingstruct;
         eval(['load ''',recycledresultsfile,''' tfinteqdiff extraoutputdata timingstruct'])
         EDinputdatahash = EDinteqinputhash;
        t01 = etime(clock,t00);
        timingdata = timingstruct.integralequation(2:5);
        currenttimingstruct.integralequation = [t01 timingdata];
        timingstruct = currenttimingstruct;
        eval(['save ''',desiredname,'''  tfinteqdiff extraoutputdata recycledresultsfile timingstruct EDsettings EDinputdatahash'])        
    else
        [tfinteqdiff,timingdata,extraoutputdata,EDinputdatahash] = EDintegralequation_convex_tf(filehandlingparameters,...
            envdata,planedata,edgedata,edgetoedgedata,Hsubmatrixdata,Sdata,Sinputdata.doaddsources,Sinputdata.sourceamplitudes,Sinputdata.doallSRcombinations,...
            Rdata,controlparameters,EDversionnumber);
        recycledresultsfile =  '';
        t01 = etime(clock,t00);
        timingstruct.integralequation = [t01 timingdata];
        eval(['save ''',desiredname,''' tfinteqdiff  extraoutputdata timingstruct recycledresultsfile EDsettings EDinputdatahash'])
    end

    EDpostfunctext('EDintegralequation_convex_tf',t01,existingfilename,...
        filehandlingparameters,fid,'',...
        ['      ',int2str(nfrequencies),' frequencies. Diffraction order: ',DIF,'.'],...
        ['      H-matrix: ',num2str(timingdata(1)),' s. Q_firstterm: ',num2str(timingdata(2)),' s (per freq.)'],...
        ['      Qfinal: ',num2str(timingdata(3)),' s. P at receiver: ',num2str(timingdata(4)),' s.'])    
    timingstruct.integralequation = timingdata;

else
    textline1 = ['   EDintegralequation_convex_tf was not run since difforder = ',DIF,' or docalctf was set to 0.'];
    EDmessage(filehandlingparameters,'sf',fid,1,'',textline1);    
    timingstruct.integralequation = [0 0 0 0 0];
    if Sinputdata.doaddsources == 1 || nsources == 1
        tfinteqdiff = zeros(nfrequencies,nreceivers);
    else
        tfinteqdiff = zeros(nfrequencies,nreceivers,nsources);        
    end
    extraoutputdata = [];
    desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'];
    eval(['save ''',desiredname,''' tfinteqdiff extraoutputdata timingstruct'])    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the results for the surface sound pressure, and propagate to the
% field points = the original receiver points.
% 
% If directsound_fieldpoints = 1, then we also have to add the
% **free-field** component (not direct sound!), from the original source
% to the fieldpoints. 

EDmessage(filehandlingparameters,'s',fid,1,'',['Propagating the surface pressure to the fieldpoints = original receivers']);
t00 = clock;    

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'];
eval(['load ''',desiredname,''' tfdirect tfgeom tfdiff']);
desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'];
eval(['load ''',desiredname,''' tfinteqdiff']);
tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;

[distances,cosfi] = EDcalccosfi(fieldpoints,Rinputdata.coordinates,Rinputdata.nvecs);
nfieldpoints = size(fieldpoints,1);
nsurfacepoints = size(Rinputdata.coordinates,1);

if directsound_fieldpoints == 1
   ffdist = EDcalcdist(fieldpoints,Sinputdata.coordinates);
   if any(Sinputdata.sourceamplitudes~=1)
      error(['ERROR: sourceamplitudes have not been implemented yet in EDmain_convexESIEBEM']) 
   end
end

if Sinputdata.doaddsources == 1 || nsources == 1
    tffieldpoints = zeros(nfrequencies,nfieldpoints);
else
    tffieldpoints = zeros(nfrequencies,nfieldpoints,nsources);    
end
kvec = 2*pi*controlparameters.frequencies(:)./envdata.cair;
kvec = kvec(:,ones(1,nsurfacepoints));
onesvec = ones(1,nfrequencies);

if Sinputdata.doaddsources == 1 || nsources == 1
    % tftot has size [nfreq,nsurfacerec]
    % tffieldpoints has size [nfreq,nfieldpoints]    
    %   becomes tftot
    for ii = 1:nfieldpoints    
        tfonefp = -tftot.*exp(-1i*kvec.*distances(:,ii*onesvec).')./distances(:,ii*onesvec).'.*(1i*kvec + 1./distances(:,ii*onesvec).').*cosfi(:,ii*onesvec).'.*Rinputdata.weights(:,onesvec).'/4/pi;
        tffieldpoints(:,ii) = sum(tfonefp.').';
    end
    tftot = tffieldpoints;
    
    if directsound_fieldpoints == 1
        kvec = 2*pi*controlparameters.frequencies(:)./envdata.cair;
        kvec = kvec(:,ones(1,nfieldpoints));
        for jj = 1:nsources
            tftot = tftot + exp(-1i*kvec.*(ffdist(:,onesvec*jj).'-controlparameters.Rstart))./ffdist(:,onesvec*jj).';
        end
    end
else
    % tftot has size [nfreq,nsurfacerec,nsources]
    % tffieldpoints has size [nfreq,nfieldpoints,nsources]    
    %   becomes tftot
    for ii = 1:nfieldpoints    
        transfermatrix = -exp(-1i*kvec.*distances(:,ii*onesvec).')./distances(:,ii*onesvec).'.*(1i*kvec + 1./distances(:,ii*onesvec).').*cosfi(:,ii*onesvec).'.*Rinputdata.weights(:,onesvec).'/4/pi;
        for jj = 1:nsources
            if Sinputdata.doallSRcombinations == 1 || jj == ii
                tfonefp = squeeze(tftot(:,:,jj)).*transfermatrix;
                tffieldpoints(:,ii,jj) = sum(tfonefp.').';
            end
        end
    end
    tftot = tffieldpoints;

    checkvalue = (size(tftot,2)==1) || (size(tftot,3)==1);
    if directsound_fieldpoints == 1
        for kk = 1:nfrequencies
            kvec = 2*pi*controlparameters.frequencies(kk*ones(nfieldpoints,nsources))./envdata.cair;
            if checkvalue == 1 
                tftot(kk,:,:) = squeeze(tftot(kk,:,:)).' + exp(-1i*kvec.*(ffdist-controlparameters.Rstart))./ffdist;            
            else
                tftot(kk,:,:) = squeeze(tftot(kk,:,:)) + exp(-1i*kvec.*(ffdist-controlparameters.Rstart))./ffdist;                            
            end
        end
    end
    
end

t01 = etime(clock,t00);
timingstruct.esiebem = t01;
        
desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfesiebem.mat'];
eval(['save ''',desiredname,''' tftot timingstruct']);

if filehandlingparameters.savelogfile == 1
    fwrite(fid,['   Propagating the surface pressure to the fieldpoints = original receivers. Time: ',num2str(t01),' s',lineending],'char');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if filehandlingparameters.savelogfile == 1
    fclose(fid);
end
