function EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDmain_convex - Calculates the specular and edge diffraction IRs and/or TFs
% for convex, rigid scattering objects, using the image source method for
% the specular reflection and the explicit integral formulation for first-order
% diffraction. For TFs, the edge source integral equation is used for the
% second- and higher-order diffraction whereas for IRs, the "restarted-
% order" approach is used.
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
%                       .textshift           (default: 0)
%
% Output parameters:
%   EDres               (optional) If no output parameter is specified when
%   EDmain_convex is called, then the results are saved in files. If this
%   output parameter is specified, EDres will be a struct which contains
%   these fields (the irrelevant ones will have zeros):
%                       irdirect, irgeom, irdiff, irhod, irtot
%                       tfdirect, tfgeom, tfdiff, tfinteqdiff, tftot
%                       tftot_ESIEBEM
%                       
%                       xxdirect contains the direct sound
%                       xxgeom contains the specular reflection
%                       xxdiff contains the first-order diffraction
%                       irhod contains the second- and higher-order
%                           diffraction; either as a sum of all orders or with
%                           one cell for each order
%                       tfinteqdiff contains the sum of second- and
%                           higher-order diffraction
%                       xxtot is the sum of the individual sound field
%                           terms listed above
%                       tftot_ESIEBEM is the total sound field, when
%                           computed with the ESIEBEM method
% 
% Uses the functions EDgetversion,EDcheckinputstructs, EDreadcad, EDreadgeomatrices,
% EDedgeo, EDSorRgeo, EDfindconvexGApaths, EDmakefirstordertfs, EDed2geo,
% EDmessage, EDpostfunctext, EDinteg_submatrixstructure, EDintegralequation_convex_tf from EDtoolbox
% Uses the functions DataHash from Matlab Central
% 
% Peter Svensson 17 Apr 2024 (peter.svensson@ntnu.no)
%
% EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

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
% 14 May 2018 Cleaned up the lineending.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 25 Aug. 2021 Small improvement: introduced the parameter
% calcfirstorderdiff in the hash for EDmakefirstordertfs. Makes recycling
% possible a bit more often.
% 27 Sep. 2023 Introduced the functions EDmessage and EDpostfunctext to
% make this main function more readable. 
% 28 Sep. 2023 Combined the three EDmain_convexESIE, EDmain_convexESIEBEM,
% and EDmain_convex_time into one main function. Also made changes to use
% the "version 2" of functions, which reduces the code in this main
% function substantially.
% 29 Sep. 2023 Included also the ir-calculations (not included the day
% before). Included the ESIEBEM steps as well.
% 3 Oct. 2023 Added the optional output parameter EDres. Also added the
% saving of a settings hash, "EDsettingshash" in all results files.
% 5 Oct. 2023 Fixed a small mistake: the wrong time data was saved for the
% EDintegralequation_convex_tf call.
% 9 Oct. 2023 Fixed a small mistake: when the irhod cells were summed, the
% first irhod cell is empty because first-order diffraction is stored 
% separately.
% 12 Oct. 2023 Introduced the possibility for a free-field case
% 13 Oct. 2023 Small modification: only direct sound saved for the
% free-field case
% 27 Oct. 2023 v0.4: quite large changes, with using a single struct for S
% and R, instead of Sinputdata/Sdata. Also, entire structs are given as
% inputs to functions, so function calls are all modified.
% 1 Nov. 2023 Added an empty .offedges field to edgedata, when a
% freefieldcase is run.
% 29 Nov. 2023 Fixed small bug: the HOD irs lengths were not handled
% correctly when summed with irdirect etc.
% 24 March 2024 Changed the call syntax for the function EDedgeo (because
% one more input parameter was added to that function).
% 16 Apr. 2024 Created a new obligatory output file: _settings.mat
% which contains the six input structs (after EDcheckinputstructs has been
% run) and the EDversionnumber.
% 17 Apr 2024 Fixed a small bug that caused an error for IR of difforder 2
% or higher, and controlparameters.savealldifforders = 0.

[EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input data, assign default values if needed

if nargin < 6
    filehandlingparameters = struct('showtext',0);
    if nargin < 5
        error('ERROR: The first five input parameters to EDmain_convex must be specified')  
    end
end

if ispc == 1
   lineending = [13,10];
else
    lineending = 10;
end
[geoinputdata,Snewdata,Rnewdata,envdata,controlparameters,filehandlingparameters] ...
    = EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the input structs in the _settings file.

Sinputdata = Snewdata;
Rinputdata = Rnewdata;

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_settings.mat'];

eval(['save(''',desiredname,''',''geoinputdata'',''Sinputdata'',''Rinputdata'',''envdata'',''controlparameters'',''filehandlingparameters'',''EDversionnumber'');'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textline1 = ['EDmain_convex, v.',num2str(EDversionnumber),' (',lastsavedate,')'];
textline2 = ['filestem for results: ',filehandlingparameters.filestem];
textline3 = ['directory for results: ',filehandlingparameters.outputdirectory];
fid = EDmessage(filehandlingparameters,'sf',0,1,'',textline1,textline2,textline3);

if controlparameters.docalctf_ESIEBEM == 1
    textline1 = ['Calculating transfer functions with the ESIEBEM method'];
else
    if controlparameters.docalctf == 1
        textline1 = ['Calculating transfer functions with the ESIE method'];
    else
        if controlparameters.docalcir == 1
            textline1 = ['Calculating impulse responses with the repeated-order method'];
        end
    end
end
fid = EDmessage(filehandlingparameters,'sf',fid,1,'',textline1,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

if geoinputdata.freefieldcase == 0
   if isfield(geoinputdata,'geoinputfile')
        EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the CAD file']);
     [planedata,extraCATTdata,elapsedtimecadgeo] = EDreadcad(filehandlingparameters,geoinputdata.geoinputfile,0);
        textline1start = ['EDreadcad: '];
        t01 = elapsedtimecadgeo;
    else
        EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the input geometry matrices']);
        [planedata,elapsedtimegeomatrices] = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners,geoinputdata.planerefltypes);         
        textline1start = ['EDreadgeomatrices: '];
        t01 = elapsedtimegeomatrices;
   end
    %%%% disp('WARNING! Test for geometry convexity has been turned off!!')
    if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
        error('ERROR: EDmain_convex can only be used for convex scatterers, including a single thin plate')
    end
    ncorners = size(planedata.corners,1);
    nplanes = (size(planedata.planecorners,1));
    Timestring = [' Time: ',num2str(t01(1)),' s'];
    if length(t01) > 1, Timestring = [Timestring,' (Orig.: ',num2str(t01(2)), 's)'];  end
    textline1 = [textline1start,int2str(ncorners),' corners and ',int2str(nplanes),' planes.',Timestring];
    
    EDmessage(filehandlingparameters,'sf',fid,1,'',textline1);
    timingstruct = struct('geoinput',t01);
else
    planedata = struct('corners','','planecorners','');
    timingstruct = struct('geoinput',0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the ESIEBEM method should be used, store the original receiver
% points as "fieldpoints", and generate new
% receiver points at the surfaces of the polyhedron.

nsources = size(Snewdata.coordinates,1);

if controlparameters.docalctf_ESIEBEM == 1
    if nsources > 1
        error('ERROR: Unfortunately, ESIEBEM does nt yet support several source positions')
    end
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the new surface receiver positions']);
    t00 = clock;
    fieldpoints = Rinputdata.coordinates;
    doallSRcombinations_fieldpoints = Snewdata.doallSRcombinations;
    directsound_fieldpoints = controlparameters.directsound;
    
    [surfacerecs,surfacerecnvecs,surfacerecweights] = EDgensurfreceivers(planedata,...
        controlparameters.surfacegaussorder,EDversionnumber);
    Rnewdata.coordinates = surfacerecs;
    Rnewdata.nvecs = surfacerecnvecs;
    Rnewdata.weights = surfacerecweights;
    Snewdata.doallSRcombinations = 1;
    controlparameters.directsound = 1;
    t01 = etime(clock,t00);
    timingstruct.surfrecs = t01;
        
    textline1 = ['       Created ',int2str(size(Rnewdata.coordinates,1)),' surface receiver positions'];
    EDpostfunctext('EDgensurfreceivers',t01,'',...
        filehandlingparameters,fid,'',textline1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the planedata struct and create an edgedata struct

if geoinputdata.freefieldcase == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the edgedata struct']);
%     [edgedata,planedata,elapsedtimeedgeo,existingfilename]...
%         = EDedgeo(planedata,geoinputdata.firstcornertoskip,...
%         geoinputdata.listofcornerstoskip,geoinputdata.planeseesplanestrategy,...
%         EDversionnumber,filehandlingparameters);
    [edgedata,planedata,elapsedtimeedgeo,existingfilename]...
        = EDedgeo(planedata,geoinputdata,...
        EDversionnumber,filehandlingparameters);
    nedges = size(edgedata.edgecorners,1);
    EDpostfunctext('EDedgeo',elapsedtimeedgeo,existingfilename,...
        filehandlingparameters,fid,'',[int2str(nedges),' edges.'])
    timingstruct.edgedata = elapsedtimeedgeo(1);
else
     edgedata = struct('edgecorners','','offedges','');
     desiredname = [filehandlingparameters.outputdirectory,filesep,...
	    filehandlingparameters.filestem,'_eddata.mat'];
     elapsedtimeedgeo = 0;
     EDinputdatahash = '00000000';
     eval(['save(''',desiredname,''',''planedata'',''edgedata'',''EDinputdatahash'',''elapsedtimeedgeo'');'])
     timingstruct.edgedata = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Sdata struct = the Sinputdata struct + some related
% parameters, like visibility

if geoinputdata.freefieldcase == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the Sdata struct']);
    [Snewdata,elapsedtimeSgeo,existingfilename] = EDSorRgeo(planedata,edgedata,Snewdata,...
         'S',EDversionnumber,filehandlingparameters);
    if any(any(Snewdata.visplanesfroms)) == 0
        error('ERROR: No source can see any of the planes of the scattering object. ',...
            'Please check your geometry.') 
    end
    EDpostfunctext('EDSorRgeo',elapsedtimeSgeo,existingfilename,...
        filehandlingparameters,fid,'',[int2str(nsources),' sources.'])
    timingstruct.Sdata = elapsedtimeSgeo;
else
    Snewdata.visplanesfroms = [];
    Snewdata.vispartedgesfroms = [];
    Sdata = Snewdata;
     desiredname = [filehandlingparameters.outputdirectory,filesep,...
	    filehandlingparameters.filestem,'_Sdata.mat'];
     elapsedtimeSRgeo = 0;
     EDinputdatahash = '00000000';
     eval(['save(''',desiredname,''',''Sdata'',''EDinputdatahash'',''elapsedtimeSRgeo'');'])
    timingstruct.Sdata = 0;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Rdata struct

if geoinputdata.freefieldcase == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the Rdata struct']);
    [Rnewdata,elapsedtimeRgeo,existingfilename] = EDSorRgeo(planedata,edgedata,Rnewdata,...
         'R',EDversionnumber,filehandlingparameters);
    if any(any(Rnewdata.visplanesfromr)) == 0
       error('ERROR: No receiver can see any of the planes of the scattering object. ',...
             'Please check your geometry.') 
    end
    nreceivers = size(Rnewdata.coordinates,1);
    EDpostfunctext('EDSorRgeo',elapsedtimeRgeo,existingfilename,...
        filehandlingparameters,fid,'',[int2str(nreceivers),' receivers.'])
    timingstruct.Rdata = elapsedtimeRgeo;
else
    Rnewdata.visplanesfromr = [];
    Rnewdata.vispartedgesfromr = [];
    Rdata = Rnewdata;
     desiredname = [filehandlingparameters.outputdirectory,filesep,...
	    filehandlingparameters.filestem,'_Rdata.mat'];
     elapsedtimeSRgeo = 0;
     EDinputdatahash = '00000000';
     eval(['save(''',desiredname,''',''Rdata'',''EDinputdatahash'',''elapsedtimeSRgeo'');'])
    timingstruct.Rdata = 0;
end

DIF = int2str(controlparameters.difforder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the paths for direct sound, first-order specular, and first-order
% diffraction paths.

if controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Find the (first-order) GA paths']);
    [firstorderpathdata,elapsedtimefindpaths,existingfilename] = ...
        EDfindconvexGApaths(planedata,edgedata,Snewdata,Rnewdata,...
        controlparameters,EDversionnumber,filehandlingparameters);
    nspecs = size(firstorderpathdata.specrefllist,1);
    ndiffs = prod(size(firstorderpathdata.diffpaths));
    EDpostfunctext('EDfindconvexGApaths',elapsedtimefindpaths,existingfilename,...
        filehandlingparameters,fid,'',[int2str(nspecs),' spec. paths and ',...
        int2str(ndiffs),' diff. paths.'])
    timingstruct.findpaths = elapsedtimefindpaths;
else
    timingstruct.findpaths = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the first-order specular, and first-order
% diffraction tfs.

if ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)...
        && controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Generate the (first-order) GA and diff tfs']);
    nfrequencies = length(controlparameters.frequencies);
    [tfdirect,tfgeom,tfdiff,timingdata,elapsedtimemaketfs,existingfilename] ...
	    = EDmakefirstordertfs(firstorderpathdata,...
        planedata,edgedata,Snewdata,Rnewdata,envdata,controlparameters,...
        EDversionnumber,filehandlingparameters);  
    EDpostfunctext('EDmakefirstordertfs',elapsedtimemaketfs,existingfilename,...
        filehandlingparameters,fid,'',...
        [int2str(nfrequencies),' frequencies.'],...
        ['Direct sound: ',num2str(timingdata(1)),'s, spec.refl.: ',...
        num2str(timingdata(2)),'s, first-order diff.: ',num2str(timingdata(3)),'s'])
    timingstruct.maketfs = elapsedtimemaketfs;
else
    timingstruct.maketfs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the first-order specular, and first-order
% diffraction irs.

if controlparameters.docalcir == 1 && controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Generate the (first-order) GA and diff irs']);
    [irdirect,irgeom,irdiff,timingdata,elapsedtimemakeirs,existingfilename] ...
        = EDmakefirstorderirs(firstorderpathdata,...
        geoinputdata,edgedata,Snewdata,Rnewdata,envdata,controlparameters,...
        EDversionnumber,filehandlingparameters); 
    EDpostfunctext('EDmakefirstorderirs',elapsedtimemakeirs,existingfilename,...
        filehandlingparameters,fid,'',...
        ['fs = ',num2str(controlparameters.fs)],...
        ['Direct sound: ',num2str(timingdata(1)),'s, spec.refl.: ',...
        num2str(timingdata(2)),'s, first-order diff.: ',num2str(timingdata(3)),'s'])
    timingstruct.makeirs = elapsedtimemakeirs;
else
    irdirect = []; irgeom = []; irdiff = [];
    timingstruct.makeirs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add edge-to-edge visibility data

if controlparameters.difforder > 1
    % Temporarily add a field 'specorder'
    controlparametersmod = controlparameters;
    controlparametersmod.specorder = 1;
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the edgetoedgedata struct']);
    [edgetoedgedata,elapsedtimeed2geo,existingfilename] = EDed2geo(planedata,edgedata,Snewdata,...
        Rnewdata,controlparametersmod,EDversionnumber,1,filehandlingparameters);
    EDpostfunctext('EDed2geo',elapsedtimeed2geo,existingfilename,...
        filehandlingparameters,fid,'')
    timingstruct.edgetoedgedata = elapsedtimeed2geo;
else
    timingstruct.edgetoedgedata = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for the integral equation: set up the submatrix structure

if controlparameters.difforder > 1 && ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)
    EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the integral equation Hsubmatrixdata struct']);
    [Hsubmatrixdata,elapsedtimesubmatrix,existingfilename] = ...
	    EDinteg_submatrixstructure(edgedata,edgetoedgedata,...
        controlparameters,EDversionnumber,filehandlingparameters);
    edgeelemsizes = edgedata.edgelengthvec./Hsubmatrixdata.nedgeelems;
    meanelemsize = mean(edgeelemsizes);
    nonzeroelements = sum(prod(Hsubmatrixdata.nedgeelems(Hsubmatrixdata.edgetripletlist),2));   
    EDpostfunctext('EDinteg_submatrixstructure',elapsedtimesubmatrix,existingfilename,...
        filehandlingparameters,fid,'',...
        [int2str(Hsubmatrixdata.nuniquesubmatrices),' submatrices, out of ',int2str(size(Hsubmatrixdata.edgetripletlist,1)),', to compute.'],...
        ['Edges discretized with: ',int2str(min(Hsubmatrixdata.nedgeelems)),'-',...
        int2str(max(Hsubmatrixdata.nedgeelems)),' gauss points. ',int2str(Hsubmatrixdata.bigmatrixendnums(end)),' edge source signals to compute.'],...
        ['The IE matrix has ',int2str(nonzeroelements),' non-zero elements. Some may be identical due to symmetries'])    
    timingstruct.submatrixdata = elapsedtimesubmatrix;
else
    timingstruct.submatrixdata = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contribution with the integral equation

if controlparameters.difforder > 1 && ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)
    textline2 = ['      ',int2str(length(controlparameters.frequencies)),' frequencies. Diffraction order: ',DIF];
    EDmessage(filehandlingparameters,'s',fid,1,'',['Calculating the HOD contribution with the FD integral equation'],textline2);
    [tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
        EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,...
    	Snewdata,Rnewdata,envdata,controlparameters,EDversionnumber,filehandlingparameters);
    EDpostfunctext('EDintegralequation_convex_tf',elapsedtimehodtf,existingfilename,...
        filehandlingparameters,fid,'',...
        [int2str(nfrequencies),' frequencies. Diffraction order: ',DIF,'.'],...
        ['H-matrix: ',num2str(timingdata(1)),' s. Q_firstterm: ',num2str(timingdata(2)),' s (per freq.)'],...
        ['Qfinal: ',num2str(timingdata(3)),' s. P at receiver: ',num2str(timingdata(4)),' s.'])    
    timingstruct.integralequation = elapsedtimehodtf;
else
    tfinteqdiff = [];
    timingstruct.integralequation = [0 0 0 0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate the ESIEBEM surface pressure to the field points = the
% original receiver points.

if controlparameters.docalctf_ESIEBEM == 1     
    EDmessage(filehandlingparameters,'s',fid,1,'',['Calculating the field at the field points with the ESIEBEM method']);
    [tftot_ESIEBEM,elapsedtimeESIEBEMpropagate,existingfilename] = ...
        EDpropagateESIEBEM(tfdirect+tfgeom+tfdiff+tfinteqdiff,...
        Rnewdata,fieldpoints,Snewdata,directsound_fieldpoints,...
        controlparameters,envdata,EDversionnumber,filehandlingparameters);
    EDpostfunctext('EDpropagateESIEBEM',elapsedtimeESIEBEMpropagate,existingfilename,...
        filehandlingparameters,fid,'',...
        [int2str(nfrequencies),' frequencies. Diffraction order: ',DIF,'.'])
    timingstruct.esiebem = elapsedtimeESIEBEMpropagate;
else
    tftot_ESIEBEM = [];
    timingstruct.esiebem = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the higher-order diffraction paths for the ir calculation

if controlparameters.docalcir == 1 && controlparameters.difforder > 1
    EDmessage(filehandlingparameters,'s',fid,1,'',['   Finding the higher-order diffraction paths']);
    [hodpaths,hodpathsalongplane,elapsedtimehodpaths,existingfilename] ...
        = EDfindHODpaths(int8(sign(edgetoedgedata.edgeseespartialedge)),...
        sign(Snewdata.vispartedgesfroms),sign(Rnewdata.vispartedgesfromr),...
        controlparameters.difforder,EDversionnumber,filehandlingparameters);
    textline1 = [' Found higher-order diffraction paths up to order ',DIF];
    EDpostfunctext('EDfindHODpaths',elapsedtimehodpaths,existingfilename,...
        filehandlingparameters,fid,'',textline1);
    timingstruct.hodpaths = elapsedtimehodpaths;
else
    timingstruct.hodpaths = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the ir HOD contributions one order at a time

if controlparameters.difforder > 1 && controlparameters.docalcir == 1
   EDmessage(filehandlingparameters,'s',fid,1,'',['   Generating the higher-order diffraction irs ']);
     
%    elemsize = 2.^(-[0:controlparameters.difforder-1]);
%    elemsize = elemsize*2;
   [irhod,elapsedtimemakeirhod,existingfilename] = EDmakeHODirs(hodpaths,...
        hodpathsalongplane,edgedata,edgetoedgedata,Snewdata,Rnewdata,envdata,...
        controlparameters,EDversionnumber,filehandlingparameters);
    textline1 = [' Generated higher-order diffraction irs up to order ',DIF,'. '];
    EDpostfunctext('EDmakeHODirs',elapsedtimemakeirhod,existingfilename,...
        filehandlingparameters,fid,'',textline1);
    timingstruct.irhod = elapsedtimemakeirhod;
else
    irhod = zeros(size(irdirect));
    timingstruct.irhod = 0;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 0
    EDres = struct;
    EDres.timingstruct = timingstruct;
    if controlparameters.docalcir == 1
        if geoinputdata.freefieldcase == 0
            irtot = irdirect + irgeom + irdiff;
            [nold,nrec,nsou] = size(irtot);  
            if iscell(irhod)
                ncells = length(irhod);
                if controlparameters.savealldifforders == 0  % This should correspond to ncells == 1
                    irhodsum = irhod{1};
                else
                    irhodsum = irhod{2};
                    for ii = 3:ncells
                        irhodsum = irhodsum + irhod{ii};
                    end
                end                
            else
                irhodsum = irhod;
            end
            nnew = size(irhodsum,1);
            if nnew > nold
               irdirect = [irdirect;zeros(nnew-nold,nrec,nsou)];
               irgeom   = [irgeom;zeros(nnew-nold,nrec,nsou)];
               irdiff   = [irdiff;zeros(nnew-nold,nrec,nsou)];
               irtot    = [irtot;zeros(nnew-nold,nrec,nsou)];
            elseif nnew < nold
               irhodsum = [irhodsum;zeros(nold-nnew,nrec,nsou)];
            end
            irtot = irtot + irhodsum;
            EDres.irdirect = irdirect;
            EDres.irgeom   = irgeom;
            EDres.irdiff   = irdiff;
            if controlparameters.savealldifforders == 1
                EDres.irhod = irhod;
            else
                EDres.irhod    = irhodsum;
            end
            EDres.irtot    = irtot;
        else
            EDres.irdirect = irdirect;
            EDres.irgeom = [];
            EDres.irdiff = [];
            EDres.irhod = [];
            EDres.irtot = irdirect;
        end
    else
        EDres.irdirect = [];
        EDres.irgeom = [];
        EDres.irdiff = [];
        EDres.irhod = [];
        EDres.irtot = [];
    end
    if controlparameters.docalctf == 1
        if isempty(tfinteqdiff)
            tfinteqdiff = zeros(size(tfdirect));
        end
        tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
        EDres.tfdirect = tfdirect;
        EDres.tfgeom = tfgeom;
        EDres.tfdiff = tfdiff;
        EDres.tfinteqdiff = tfinteqdiff;
        EDres.tftot = tftot;
    else
        EDres.tfdirect = [];
        EDres.tfgeom = [];
        EDres.tfdiff = [];
        EDres.tfinteqdiff = [];
        EDres.tftot = [];
    end
    if controlparameters.docalctf_ESIEBEM == 1
        EDres.tftot_ESIEBEM = tftot_ESIEBEM;
    else
        EDres.tftot_ESIEBEM = [];
    end
    EDres.geoinputdata = geoinputdata;
    EDres.Sinputdata = Sinputdata;
    EDres.Rinputdata = Rinputdata;
    EDres.envdata = envdata;
    EDres.controlparameters = controlparameters;
    EDres.filehandlingparameters = filehandlingparameters;
    EDres.EDversionnumber = EDversionnumber;
end

if filehandlingparameters.savelogfile == 1
    fclose(fid);
end
