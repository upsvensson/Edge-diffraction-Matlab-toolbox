function EDres = EDmain_convex(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDmain_convex - Computes the sound field around a convex, rigid scattering object.
% The direct sound and the first-order specular reflections are computed
% with the image source method, and first-order diffraction with the
% explicit integral formulation. Second- and higher-order diffraction is
% computed with the edge source integral equation (ESIE) for transfer
% functions, or with the "repeated-order" approach for impulse responses.
% The method is chosen with controlparameters.generalmethod: 'tf', 'ir' or
% 'tf_ESIEBEM'.
%
% The input structs are checked by EDcheckinputstructs, which also sets the
% default values, and they are saved in <filestem>_settings.mat. The
% intermediate results are saved in stage files, which are reused by later
% runs when the data they depend on are the same (see help EDdatastructs).
% Each run writes a log file, <filestem>_log_<date and time>.txt.
%
% Input parameters:
%   geoinputdata, Sinputdata, Rinputdata, envdata, controlparameters,
%   filehandlingparameters
%                   Six structs. All fields, their default values and when
%                   they are used are given by the tables in
%                   EDcheckinputstructs. The most important fields:
%                   geoinputdata       .geoinputfile, or .corners and
%                                      .planecorners
%                   Sinputdata         .coordinates (or the piston fields),
%                                      .sourceamplitudes, .doaddsources
%                   Rinputdata         .coordinates
%                   envdata            .cair
%                   controlparameters  .generalmethod, .difforder
%                                      (required), .frequencies (for 'tf'),
%                                      .fs (for 'ir'), .ngauss
%                   filehandlingparameters  .outputdirectory, .filestem
%                   filehandlingparameters can be left out if a geometry
%                   file is given.
%
% Output parameters:
%   EDres           (optional) A struct with the results, the six input
%                   structs (each with the field .inputorigin), the timing
%                   and the diagnostics; see help EDdatastructs. The struct
%                   is also saved in <filestem>_EDres.mat (without the result
%                   matrices if filehandlingparameters.savefullEDresfile =
%                   0), and can be loaded later with EDgetresults.
%
% Uses the functions EDcheckinputstructs, EDinputhash, EDinputsummary,
% EDreadcad, EDreadgeomatrices, EDgensurfreceivers, EDedgeo, EDSorRgeo,
% EDfindconvexGApaths, EDmakefirstordertfs, EDmakefirstorderirs, EDed2geo,
% EDinteg_submatrixstructure, EDintegralequation_convex_tf, EDinteqsummary,
% EDpropagateESIEBEM, EDfindHODpaths, EDmakeHODirs, EDassembleresults,
% EDmessage, EDpostfunctext, EDgetversion from the EDtoolbox
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
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
% 10 Sep. 2026 Adapted to changes in EDintegralequation_convex by storing
% the inteqdiagnostics struct in EDres.
% 17 Sep. 2026 Introduced controlparameters.generalmethod instead of
% '.docalcir' etc.
% 11 Oct. 2026 Adapted to the new EDcheckinputstructs: the choice between
% geometry file and geometry matrices is made with ischar(.geoinputfile).
% Corrected the order of the input parameters in the call to EDreadcad.
% 11 Oct. 2026 The new output inputorigin from EDcheckinputstructs (which
% tells if each input value was given by the user, a default, derived or
% not applicable) is attached to each input struct, as the field
% .inputorigin, in the _settings file and in EDres.
% 11 Oct. 2026 The results struct is now always saved, in the file
% <filestem>_EDres.mat (filehandlingparameters.saveEDresfile, default 1; see
% 12 Oct. 2026),
% together with a hash of the complete input data (from EDinputhash).
% The hash is also saved in the _settings file.
% 12 Oct. 2026 The putting together of the results was moved to the new
% function EDassembleresults. The setting saveEDresfile was replaced by
% savefullEDresfile: the _EDres file is always saved, with or without the
% result matrices.
% 12 Oct. 2026 The log file now contains all the input data, with the
% origin of each value (via EDinputsummary), a summary of the solver
% diagnostics (via EDinteqsummary), and closing lines with the total time
% and the name of the results file.
% 12 Oct. 2026 The edge discretization is checked against the highest
% frequency (points per wavelength): a log line; a note in the log file below
% 6 points per wavelength; a warning on screen and in the log below 4.
% 12 Oct. 2026 EDmakefirstorderirs is given planedata instead of
% geoinputdata (the input parameter is called planedata).
% 14 Oct. 2026 filehandlingparameters is given to EDinputhash (its field
% .savediff2result is part of the run hash).
% 14 Oct. 2026 The log file name ends with the date and time of the run
% (see EDmessage), so reruns do not overwrite earlier logs. The name is
% stored in EDres.logfilename and shown in the closing lines.
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).

[EDversionnumber,lastsavedate,lastsavetime] = EDgetversion;

t00main = clock;

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
[geoinputdata,Snewdata,Rnewdata,envdata,controlparameters,filehandlingparameters,inputorigin] ...
    = EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the input structs in the _settings file.

Sinputdata = Snewdata;
Rinputdata = Rnewdata;

desiredname = [filehandlingparameters.outputdirectory,filesep,...
	filehandlingparameters.filestem,'_settings.mat'];

% The origin record (user/default/derived/notapplicable for each field) is
% attached to each struct, as the field .inputorigin, in the saved file
% (and in EDres, at the end of this function). The structs that are used in
% the calculations do not contain that field.
inputstructnames = {'geoinputdata','Sinputdata','Rinputdata','envdata','controlparameters','filehandlingparameters'};
settingsdata = struct;
settingsdata.geoinputdata = geoinputdata;
settingsdata.Sinputdata = Sinputdata;
settingsdata.Rinputdata = Rinputdata;
settingsdata.envdata = envdata;
settingsdata.controlparameters = controlparameters;
settingsdata.filehandlingparameters = filehandlingparameters;
for ii = 1:length(inputstructnames)
    settingsdata.(inputstructnames{ii}).inputorigin = inputorigin.(inputstructnames{ii});
end
settingsdata.EDversionnumber = EDversionnumber;

% One hash value for the complete input (the five structs that define the
% problem). It is saved in the settings file and in the results file, so
% that the results for a given set of input data can be found later.
EDruninputhash = EDinputhash(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);
settingsdata.EDruninputhash = EDruninputhash;

save(desiredname,'-struct','settingsdata');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textline1 = ['EDmain_convex, v.',num2str(EDversionnumber),' (',lastsavedate,')'];
textline2 = ['filestem for results: ',filehandlingparameters.filestem];
textline3 = ['directory for results: ',filehandlingparameters.outputdirectory];
fid = EDmessage(filehandlingparameters,'sf',0,1,'',textline1,textline2,textline3);
% The log file name contains the date and time (set by EDmessage)
if fid > 0
    [~,logfilestem,logfileext] = fileparts(fopen(fid));
    logfilename = [logfilestem,logfileext];
else
    logfilename = '';
end

%if controlparameters.docalctf_ESIEBEM == 1
if strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1
    textline1 = ['Calculating transfer functions with the ESIEBEM method'];
    textline2 = ['nfrequencies = ',int2str(length(controlparameters.frequencies)),', ngauss = ',int2str(controlparameters.ngauss),...
        ', cair = ',num2str(envdata.cair)];
else
%    if controlparameters.docalctf == 1
    if strcmp(controlparameters.generalmethod,'tf') == 1
        textline1 = ['Calculating transfer functions with the ESIE method'];
        textline2 = ['nfrequencies = ',int2str(length(controlparameters.frequencies)),', ngauss = ',int2str(controlparameters.ngauss),...
        ', cair = ',num2str(envdata.cair)];
    else
 %       if controlparameters.docalcir == 1
    if strcmp(controlparameters.generalmethod,'ir') == 1 
            textline1 = ['Calculating impulse responses with the repeated-order method'];
            textline2 = ['fs = ',num2str(controlparameters.fs),', cair = ',num2str(envdata.cair)];
        end
    end
end
fid = EDmessage(filehandlingparameters,'sf',fid,1,'',textline1,textline2,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write all the input data to the log file (not on the screen): one line
% per field, with its value and where the value came from (user, default,
% derived or not applicable).

if filehandlingparameters.savelogfile == 1
    EDmessage(filehandlingparameters,'f',fid,1,'',' ',['Hash of the input data: ',EDruninputhash],'');
    for ii = 1:length(inputstructnames)
        inputtextlines = EDinputsummary(settingsdata.(inputstructnames{ii}),inputstructnames{ii});
        EDmessage(filehandlingparameters,'f',fid,1,'',['Input data: ',inputstructnames{ii}],inputtextlines{:},'');
    end
end
clear settingsdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the input CAD-file, or input matrices, and create the planedata struct

if geoinputdata.freefieldcase == 0
   % After EDcheckinputstructs, geoinputdata.geoinputfile is either a file
   % name or NaN (when the geometry was given as matrices).
   if ischar(geoinputdata.geoinputfile)
        EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the CAD file']);
     [planedata,extraCATTdata,elapsedtimecadgeo] = EDreadcad(geoinputdata.geoinputfile,0,filehandlingparameters);
        textline1start = ['EDreadcad: '];
        t01 = elapsedtimecadgeo;
    else
        EDmessage(filehandlingparameters,'s',fid,1,'',['Creating the planedata struct from the input geometry matrices']);
        [planedata,elapsedtimegeomatrices] = EDreadgeomatrices(geoinputdata.corners,geoinputdata.planecorners,geoinputdata.planerefltypes);         
        textline1start = ['EDreadgeomatrices: '];
        t01 = elapsedtimegeomatrices;
   end
    if geoinputdata.allownonconvex == 0
        if isempty(strfind(planedata.modeltype,'convex_ext')) && isempty(strfind(planedata.modeltype,'singleplate'))
         error('ERROR: EDmain_convex can only be used for convex scatterers, including a single thin plate')
        end
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

%if controlparameters.docalctf_ESIEBEM == 1
if strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1
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
    % [edgedata,planedata,elapsedtimeedgeo,existingfilename,edgedatahash]...
    %     = EDedgeo(planedata,geoinputdata.firstcornertoskip,...
    %     geoinputdata.listofcornerstoskip,geoinputdata.planeseesplanestrategy,...
    %     EDversionnumber,filehandlingparameters);
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

%if ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)
if ( strcmp(controlparameters.generalmethod,'tf') == 1 || strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1)
    nfrequencies = length(controlparameters.frequencies);
    if controlparameters.skipfirstorder == 0
        EDmessage(filehandlingparameters,'s',fid,1,'',['Generate the (first-order) GA and diff tfs']);
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
        if nsources == 1
            tfdirect = zeros(nfrequencies,nreceivers);
            tfgeom = zeros(nfrequencies,nreceivers);
            tfdiff = zeros(nfrequencies,nreceivers);
        else
            if nreceivers == 1
                tfdirect = zeros(nfrequencies,nsources);
                tfgeom = zeros(nfrequencies,nsources);
                tfdiff = zeros(nfrequencies,nsources);
            else
                tfdirect = zeros(nfrequencies,nreceivers,nsources);
                tfgeom = zeros(nfrequencies,nreceivers,nsources);
                tfdiff = zeros(nfrequencies,nreceivers,nsources);
            end
        end
    end
else
    timingstruct.maketfs = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the first-order specular, and first-order
% diffraction irs.

%if controlparameters.docalcir == 1 && controlparameters.skipfirstorder == 0
if strcmp(controlparameters.generalmethod,'ir') == 1 && controlparameters.skipfirstorder == 0
    EDmessage(filehandlingparameters,'s',fid,1,'',['Generate the (first-order) GA and diff irs']);
    [irdirect,irgeom,irdiff,timingdata,elapsedtimemakeirs,existingfilename] ...
        = EDmakefirstorderirs(firstorderpathdata,...
        planedata,edgedata,Snewdata,Rnewdata,envdata,controlparameters,...
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

%if controlparameters.difforder > 1 && ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)
if controlparameters.difforder > 1 && ( strcmp(controlparameters.generalmethod,'tf') == 1 || strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1)
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
        int2str(max(Hsubmatrixdata.nedgeelems)),' gauss points. ',int2str(Hsubmatrixdata.Qvecendnums(end)),' edge source signals to compute.'],...
        ['The IE matrix has ',int2str(nonzeroelements),' non-zero elements. Some may be identical due to symmetries'])    
    timingstruct.submatrixdata = elapsedtimesubmatrix;

    % Check the edge discretization against the highest frequency: the
    % number of edge points per wavelength. The number of points on an
    % edge is ngauss*(edge length)/(longest edge length), so the points per
    % wavelength are (nearly) the same for all edges; the smallest value is
    % reported. Two levels, based on convergence tests (ngauss = 8 vs 16)
    % for a cube, Oct. 2026:
    %   >= 6 points per wavelength: the error is at a frequency-independent
    %        floor (about 0.4% relative to the largest amplitude).
    %   4-6: the overall accuracy is still at that floor, but the relative
    %        error at low-level receivers (deep minima, shadow regions)
    %        grows, to about 1-2%. A note is written to the log file.
    %   < 4: also the overall accuracy degrades. A warning is written to
    %        the log file and shown on the screen.
    pointsperwavelength_note = 6;
    pointsperwavelength_warning = 4;
    pointspermeter = min(double(Hsubmatrixdata.nedgeelems(:))./edgedata.edgelengthvec(:));
    maxfrequency = max(controlparameters.frequencies);
    pointsperwavelength = pointspermeter*envdata.cair/maxfrequency;
    frequencylimit_note = pointspermeter*envdata.cair/pointsperwavelength_note;
    frequencylimit_warning = pointspermeter*envdata.cair/pointsperwavelength_warning;
    EDmessage(filehandlingparameters,'f',fid,1,'',['   Edge discretization: ',num2str(pointsperwavelength,3),...
        ' points per wavelength at the highest frequency (',num2str(maxfrequency),' Hz). ',...
        num2str(pointsperwavelength_note),' points per wavelength up to ',num2str(round(frequencylimit_note)),' Hz, ',...
        num2str(pointsperwavelength_warning),' up to ',num2str(round(frequencylimit_warning)),' Hz.']);
    if pointsperwavelength < pointsperwavelength_warning
        EDmessage(filehandlingparameters,'sf',fid,0,'',['   WARNING: the edges have fewer than ',num2str(pointsperwavelength_warning),...
            ' points per wavelength above ',num2str(round(frequencylimit_warning)),' Hz (',num2str(pointsperwavelength,3),...
            ' at ',num2str(maxfrequency),' Hz), so the accuracy degrades there. Consider a higher controlparameters.ngauss.']);
    elseif pointsperwavelength < pointsperwavelength_note
        EDmessage(filehandlingparameters,'f',fid,1,'',['   NOTE: the edges have fewer than ',num2str(pointsperwavelength_note),...
            ' points per wavelength above ',num2str(round(frequencylimit_note)),' Hz. The overall accuracy is not affected, ',...
            'but relative errors at low-level receivers (deep minima, shadow regions) can reach 1-2% there.']);
    end
else
    timingstruct.submatrixdata = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HOD contribution with the integral equation

%if controlparameters.difforder > 1 && ( controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1)
if controlparameters.difforder > 1 && ( strcmp(controlparameters.generalmethod,'tf') == 1 || strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1)
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

    % A summary of the solver diagnostics goes to the log file. A warning
    % (if GMRES did not converge for some frequency) is also shown on the
    % screen.
    [inteqtextlines,inteqwarninglines] = EDinteqsummary(extraoutputdata.inteqdiagnostics);
    if ~isempty(inteqtextlines)
        EDmessage(filehandlingparameters,'f',fid,1,'','   Diagnostics:',inteqtextlines{:});
    end
    if ~isempty(inteqwarninglines)
        EDmessage(filehandlingparameters,'sf',fid,0,'',['   ',inteqwarninglines{1}]);
    end
%    timingstruct.integralequation = elapsedtimehodtf;
    timingstruct.integralequation = timingdata;
else
    tfinteqdiff = [];
    extraoutputdata = struct;
    extraoutputdata.inteqdiagnostics = [];
    timingstruct.integralequation = [0 0 0 0 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate the ESIEBEM surface pressure to the field points = the
% original receiver points.

%if controlparameters.docalctf_ESIEBEM == 1     
if strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1     
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

%if controlparameters.docalcir == 1 && controlparameters.difforder > 1
if strcmp(controlparameters.generalmethod,'ir') == 1 && controlparameters.difforder > 1
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

%if controlparameters.difforder > 1 && controlparameters.docalcir == 1
if strcmp(controlparameters.generalmethod,'ir') == 1 && controlparameters.difforder > 1
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

% Collect the results.
%
% EDmeta contains everything except the (possibly large) result matrices:
% timing data, diagnostics, the six input structs (each with its origin
% record), the version number and the hash of the input data.
% The result matrices are put together by EDassembleresults.
%
% The file <filestem>_EDres.mat is always saved:
%   savefullEDresfile = 1: the file contains the complete results struct.
%   savefullEDresfile = 0: the file contains EDmeta only; the result
%   matrices are then only in the _tf, _tfinteq, _ir and _irhod files, from
%   where EDgetresults can load them.

EDmeta = struct;
EDmeta.timingstruct = timingstruct;
if strcmp(controlparameters.generalmethod,'tf') == 1
    EDmeta.inteqdiagnostics = extraoutputdata.inteqdiagnostics;
end
EDmeta.geoinputdata = geoinputdata;
EDmeta.Sinputdata = Sinputdata;
EDmeta.Rinputdata = Rinputdata;
EDmeta.envdata = envdata;
EDmeta.controlparameters = controlparameters;
EDmeta.filehandlingparameters = filehandlingparameters;
EDmeta.EDversionnumber = EDversionnumber;
for ii = 1:length(inputstructnames)
    EDmeta.(inputstructnames{ii}).inputorigin = inputorigin.(inputstructnames{ii});
end
EDmeta.inputhash = EDruninputhash;
EDmeta.logfilename = logfilename;

buildfullresults = (nargout > 0) || (filehandlingparameters.savefullEDresfile == 1);
if buildfullresults
    resultparts = struct;
    if strcmp(controlparameters.generalmethod,'ir') == 1
        resultparts.irdirect = irdirect;
        resultparts.irgeom = irgeom;
        resultparts.irdiff = irdiff;
        resultparts.irhod = irhod;
    elseif strcmp(controlparameters.generalmethod,'tf') == 1
        resultparts.tfdirect = tfdirect;
        resultparts.tfgeom = tfgeom;
        resultparts.tfdiff = tfdiff;
        resultparts.tfinteqdiff = tfinteqdiff;
    end
    resultparts.tftot_ESIEBEM = tftot_ESIEBEM;

    EDresult = EDassembleresults(resultparts,controlparameters,geoinputdata);
    clear resultparts
    metanames = fieldnames(EDmeta);
    for ii = 1:length(metanames)
        EDresult.(metanames{ii}) = EDmeta.(metanames{ii});
    end
end

desiredname = [filehandlingparameters.outputdirectory,filesep,...
    filehandlingparameters.filestem,'_EDres.mat'];
resfiledata = struct;
resfiledata.EDruninputhash = EDruninputhash;
if filehandlingparameters.savefullEDresfile == 1
    resfiledata.EDres = EDresult;
    resfiledata.EDresisfull = 1;
else
    resfiledata.EDres = EDmeta;
    resfiledata.EDresisfull = 0;
    % Note the hashes of the files that hold the result matrices, so that
    % EDgetresults can detect if one of them has been overwritten by
    % another calculation.
    stagefilesuffixes = {'_tf','_tfinteq','_ir','_irhod','_tfESIEBEM'};
    EDstagefilehashes = struct;
    for ii = 1:length(stagefilesuffixes)
        stagefilename = [filehandlingparameters.outputdirectory,filesep,...
            filehandlingparameters.filestem,stagefilesuffixes{ii},'.mat'];
        if exist(stagefilename,'file') == 2
            variablesinfile = whos('-file',stagefilename);
            if any(strcmp({variablesinfile.name},'EDinputdatahash'))
                hashinfile = load(stagefilename,'EDinputdatahash');
                EDstagefilehashes.(stagefilesuffixes{ii}(2:end)) = hashinfile.EDinputdatahash;
            end
        end
    end
    resfiledata.EDstagefilehashes = EDstagefilehashes;
end
save(desiredname,'-struct','resfiledata');
clear resfiledata

if nargout > 0
    EDres = EDresult;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The closing lines in the log file

if filehandlingparameters.savefullEDresfile == 1
    textline2 = ['Results file: ',filehandlingparameters.filestem,'_EDres.mat (complete)'];
else
    textline2 = ['Results file: ',filehandlingparameters.filestem,'_EDres.mat (without the result matrices)'];
end
if isempty(logfilename)
    textline3 = 'No log file (savelogfile = 0)';
else
    textline3 = ['Log file: ',logfilename];
end
EDmessage(filehandlingparameters,'sf',fid,1,'',...
    ['EDmain_convex finished. Total time: ',num2str(etime(clock,t00main)),' s'],textline2,textline3);

if filehandlingparameters.savelogfile == 1
    fclose(fid);
end
