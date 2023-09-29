function [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDcheckinputstructs checks the input data structs to EDmain_convex
% and sets default values.
% 
% Input parameters:
%   geoinputdata     .geoinputfile        (specified, or a file open
%                                             window is given - unless the fields .corners and .planecorners are given)
%                   .corners             (optional; alternative to
%                                             .geoinputfile. But, the use of this option
%                                             requires that filehandlingparameters.outputdirectory and .filestem are specified)
%                   .planecorners        (optional; alternative to
%                                             .geoinputfile)
%                   .planerefltypes      (optional; gives the reflection
%                                         factors of the planes: 1,0,-1.
%                                         default: ones(nplanes,1) )
%                   .firstcornertoskip   (default: 1e6)
%   Sinputdata         .coordinates         (obligatory)
%                   .doaddsources        (default: 0 = no)
%                   .sourceamplitudes    Only used if doaddsources = 1
%                                        (default:
%                                         ones(nsources,nfrequencies))
%                   .doallSRcombinations  (default: 1 = yes)
%   Rinputdata         .coordinates         (obligatory)
%   envdata         .cair                (default: 344)
%                   .rhoair              (default: 1.21)
%   controlparameters   .fs              (default: 44100) Ignored by
%                                        EDmain_convexESIE, but used by
%                                        EDmain_convexESIE_ir
%                   .directsound         (default: 1 = yes)
%                   .difforder           (default: 15)
%                   .savealldifforders   (default: 0) Used only by
%                                        EDmain_convex_time
%                   .saveindividualfirstdiff (default: 0) Used only by
%                                        EDmain_convex_time
%                   .skipfirstorder      (default: 0)
%                   .docalctf            (default: 1 for EDmain_convexESIE,
%                                                  ignored by EDmain_convexESIE_ir)
%                   .docalcir            (default: 1 for EDMain_convexESIE_ir,
%                                                  ignored by EDmain_convexESIE)
%                   .Rstart              (default: 0)
%                   .frequencies         (obligatory for EDmain_convexESIE,
%                                         ignored by EDmain_convexESIE_ir)
%                   .discretizationtype  (default: 2 = G-L)
%                   .ngauss              (default: 16, but ignored by
%                                        EDmain_convex_time)
%                   .surfacegaussorder   (default: 5 for
%                                         EDmain_convexESIEBEM; ignored by other EDmain versions)
%   filehandlingparameters    .outputdirectory  (default: the folder of the geoinputfile)  
%                   .filestem        (default: name of the cad-file, with an underscore + running integer)
%                   .savecadgeofile      (default: 0)
%                   .saveSRdatafiles     (default: 1)
%                   .saveeddatafile      (default: 1)
%                   .saveed2datafile      (default: 1)
%                   .savesubmatrixdata    (default: 1)
%                   .saveinteqsousigs     (default: 0)
%                   .loadinteqsousigs     (default: 0)
%                   .savepathsfile        (default: 1)
%                   .saveISEStree         (default: 0) Not used by
%                                         EDmain_convexESIE
%                   .savelogfile          (default: 1)
%                   .savediff2result      (default: 0)
%                   .savehodpaths         (default: 0) Used only by
%                   EDmain_convex_time
%                   .showtext             (default: 1)
% 
% Peter Svensson 28 Sep. 2023 (peter.svensson@ntnu.no)
% 
% [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 First version
% 28 Nov. 2017 Cleaned code a bit
% 29 Nov. 2017 Corrected mistake: saveSRdatafiles was called
%              saveSRinputdatafiles. Adjusted to the new indata option:
%              specified corners and planecorners matrices instead of a CAD
%              file.
% 30 Nov. 2017 Fixed a bug where docalctf, instead of docalcir, was set to
%              0. Changed defaults to create eddata and SRfiles, because of
%              plotting the model. Changed one field from logfilename to
%              savelogfile. Removed the field lineending.
% 13 Dec. 2017 Added the input field sourceamplitudes
% 12 Jan. 2018 Forced doaddsources to be 1, if the number of sources = 1.
%              Also followed one yellow recommendation: numel instead of
%              prod.
% 17 Jan 2018  Added the input field planecornertype
% 18 Jan 2018  Fixed a bug with the sourceamplitudes; they were forced to
% one by mistake. Changed default for savelogfile to 1.
% 22 Jan 2018 Removed the field planecornertype in the struct geoinputdata
% 22 Jan 2018 Moved some controlparameter fields away from the convex TF
% case: fs, nedgepoints_visibility, docalcir, specorder.
% 22 Jan 2018 Changed the defaults for saving files.
% 26 Jan 2018 Changed the defaults for saving files. Also changed the
% default output directory to include "results" in the path.
% 26 Jan 2018 Introduced Sinputdata.doallSRcombinations. Default value 1.
% 28 Jan 2018 First version for EDmain_convexESIE_ir
% 31 Jan 2018 Corrected the handling of sourceamplitudes (the freq.
% dependence was not implemented correctly.
% 6 Feb 2018 Introduced a check if number of receivers/sources was zero
% 8 Feb 2018 v 0.109 Introduced a new parameter:
% controlparameters.skipfirstorder (default = 0).
% 14 Feb 2018 v0.112 Small change, assigning a default value to filestem if
% it was not given.
% 14 Feb 2018 Added a check if the functions 'DataHash.m' and 'lgwt.m' are
% available.
% 14 Feb 2018 Added a test if the source and receiver coordinates has the
% right number of columns (3).
% 15 Feb 2018 Removed the savesetupfile. Introduced the saveed2datafile,
% default = 1. Changed defaults for savepathsfile and savesubmatrixdata to
% 1. Stopped assigning default value to saveISEStree. Introduced the
% parameter suppressresultrecycling, default = 0.
% 1 Mar 2018 Introduced the version 3: EDmain_convexESIEBEM
% 1 Mar 2018 Stripped away irrelevant fields; otherwise the DataHash
% doesn't always recognize what is identical.
% 15 Mar 2018 Added the field filehandlingparameters.savehodpaths
% 16 Mar 2018 Combined two error messages
% 21 Mar 2018 Made a few changes to the controlparameters: gave default
% values for docalctf and docalcir. Removed specorder. Introduced the
% .savealldifforders parameter. 
% 21 Mar 2018 Changed savealldifforders to controlparameters instead of
% filehandlingparameters
% 7 Apr 2018 Introduced controlparameters.saveindividualfirstdiff
% 21 Apr 2018 Removed the extra "results" directory in the output
% directory.
% 28 May 2018 Added the input field geoinputdata.planerefltypes
% 22 May 2019 Fixed an error with sourceamplitudes; they did not have the
% right orientation.
% 3 June 2020 Fixed a but with sourceamplitudes that was found by EDdebug:
% If the user had specified a constant sourceamplitudes, it wasn't expanded
% to a matrix of size [nsources, nfrequencies].
% 14 March 2021 The section "% Check the struct Rinputdata" was moved to
% before "% Check the struct Sinputdata"
% 28 Sep. 2023 Adapted to EDmain_convex which does both tf and ir

% if nargin < 7
%     disp('ERROR: the input parameter EDmaincase was not specified')
%     return
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the needed non-EDtoolbox functions are available.

file1ok = (exist('DataHash.m','file') == 2);
file2ok = (exist('lgwt.m','file') == 2);

if file1ok == 0 && file2ok == 0
   error('ERROR: Matlab can not find the functions lgwt.m and DataHash.m. Please download them from Mathworks') 
end
if file1ok == 0
   error('ERROR: Matlab can not find the function DataHash.m. Please download it from Mathworks') 
end
if file2ok == 0
   error('ERROR: Matlab can not find the function lgwt.m. Please download it from Mathworks') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct geoinputdata

if ~isstruct(geoinputdata)

	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
    [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

    CADfile = [CADfilepath,CADfile];
    geoinputdata = struct ('geoinputfile',CADfile);

    [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
    end   
    if ~isfield(filehandlingparameters,'filestem')
        filehandlingparameters.filestem = CADfilestem;
    end
end
if isfield(geoinputdata,'geoinputfile')
    [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
    end    
    if ~isfield(filehandlingparameters,'filestem')
        filehandlingparameters.filestem = CADfilestem;
    end
else
    if ~isfield(geoinputdata,'corners') || ~isfield(geoinputdata,'planecorners')
    	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
        [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

        CADfile = [CADfilepath,CADfile];
        geoinputdata.geoinputfile = CADfile;
        [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);        
        if ~isfield(filehandlingparameters,'outputdirectory')
            filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
        end    
    else
        if isfield(filehandlingparameters,'outputdirectory') == 0 || isfield(filehandlingparameters,'filestem') == 0
            error('ERROR: When you give the geometry input in the form of data matrices, you must specify filehandlingparameters.outputdirectory and .filestem')            
        end
        if ~isfield(geoinputdata,'planerefltypes')
           nplanes = size(geoinputdata.planecorners,1);
           geoinputdata.planerefltypes = ones(nplanes,1);
        end
    end
end
if ~isfield(geoinputdata,'firstcornertoskip')
    geoinputdata.firstcornertoskip = 1e6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Rinputdata

if ~isstruct(Rinputdata)
    error('ERROR 1: receiver coordinates were not specified')
end
if ~isfield(Rinputdata,'coordinates')
    error('ERROR 2: receiver coordinates were not specified')
end
nreceivers = size(Rinputdata.coordinates,1);
if nreceivers == 0
     error('ERROR 3: receiver coordinates were not specified')            
end
ncolumns = size(Rinputdata.coordinates,2);
if ncolumns ~= 3
   error(['ERROR: check your receiver coordinates; there were ',int2str(ncolumns),' columns rather than 3']) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Sinputdata, but the field .sourceamplitudes is handled
% further down, since the values depend on whether it is a TD case or an FD
% case.

if ~isstruct(Sinputdata)
    error('ERROR 1: source coordinates were not specified')
end
if ~isfield(Sinputdata,'coordinates')
    error('ERROR 2: source coordinates were not specified')
end
nsources = size(Sinputdata.coordinates,1);
if nsources == 0
     error('ERROR 3: source coordinates were not specified')            
end
if ~isfield(Sinputdata,'doaddsources')
    Sinputdata.doaddsources = 0;
end
if nsources == 1
    Sinputdata.doaddsources = 1;
end    
if ~isfield(Sinputdata,'doallSRcombinations')
    Sinputdata.doallSRcombinations = 1;
else
    if Sinputdata.doallSRcombinations == 0 && nsources~=nreceivers
        disp(['   nsources = ',int2str(nsources)])
        disp(['   nreceivers = ',int2str(nreceivers)])
        
       error('ERROR: doallSRcombinations was set to 0, but the number of sources was not the same as the number of receivers'); 
    end
end
ncolumns = size(Sinputdata.coordinates,2);
if ncolumns ~= 3
   error(['ERROR: check your source coordinates; there were ',int2str(ncolumns),' columns rather than 3']) 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct envdata

if ~isstruct(envdata)
    error('ERROR: the struct envdata was not specified')
end
if ~isfield(envdata,'cair')
    envdata.cair = 344;
end
if ~isfield(envdata,'rhoair')
    envdata.rhoair = 1.21;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct controlparameters

if ~isstruct(controlparameters)
     error('ERROR: the struct controlparameters was not specified')           
end

% First some general parameters 

if ~isfield(controlparameters,'directsound')
    controlparameters.directsound = 1;
end
if ~isfield(controlparameters,'skipfirstorder')
    controlparameters.skipfirstorder = 0;
end
if ~isfield(controlparameters,'Rstart')
    controlparameters.Rstart = 0;
end
if ~isfield(controlparameters,'difforder')
    disp('WARNING: controlparameters.difforder was not specified. It is given the value 15.')
    controlparameters.difforder = 15;
end

% Then other parameters 

if ~isfield(controlparameters,'docalctf')
    controlparameters.docalctf = 0;
    disp('WARNING: controlparameters.docalctf was not specified. It is given the value 0.')
end
if ~isfield(controlparameters,'docalcir')
    controlparameters.docalcir = 0;
    disp('WARNING: controlparameters.docalcir was not specified. It is given the value 0.')
end
if ~isfield(controlparameters,'docalctf_BEM')
    controlparameters.docalctf_BEM = 0;
    disp('WARNING: controlparameters.docalctf_BEM was not specified. It is given the value 0.')
end

if controlparameters.docalctf == 1 || controlparameters.docalctf_BEM == 1
    if ~isfield(controlparameters,'frequencies')
        error('ERROR: controlparameters.frequencies were not specified')
    end
    if ~isfield(controlparameters,'ngauss')
        disp('WARNING! controlparameters.ngauss wasnt set; it is given the default value 16.')
        controlparameters.ngauss = 16;
    end
    if ~isfield(controlparameters,'discretizationtype')
        disp('WARNING! controlparameters.discretizationtype wasnt set; it is given the default value 2.')
        controlparameters.discretizationtype = 2;
    end
else
    controlparameters.frequencies = [];
    controlparameters.ngauss = 0;
    controlparameters.discretizationtype = 0;
end
nfrequencies = length(controlparameters.frequencies);

if controlparameters.docalctf_BEM == 1
    if ~isfield(controlparameters,'surfacegaussorder')
        disp('WARNING! controlparameters.surfacegaussorder wasnt set; it is given the default value 5.')
        controlparameters.surfacegaussorder = 5;
    end
else
    controlparameters.surfacegaussorder = 0;
end

if controlparameters.docalcir == 1
    if ~isfield(controlparameters,'fs')
        disp('WARNING! controlparameters.fs wasnt set; it is given the default value 44100.')
        controlparameters.fs = 44100;
    end   
else
   controlparameters.fs = 44100;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the field .sourceamplitudes of the struct Sinputdata

sourceamplitudes_default = 0;

% We need to treat the frequency-domain cases differently from the
% time-domain cases. But, if .sourceamplitudes was not specified, it gets
% the default value 1 for all sources (and frequencies).

if ~isfield(Sinputdata,'sourceamplitudes')
%     Sinputdata.sourceamplitudes = ones(nsources,1);    
    Sinputdata.sourceamplitudes = 1;    
    sourceamplitudes_default = 1;
end

[n1,n2] = size(Sinputdata.sourceamplitudes);

if n1 > 1
    if n1 ~= nsources
        error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one row, or one row per source, but it had the size ',int2str(n1),' by ',int2str(n2)])
    else % So, we know that there is one amplitude per source. Check if there is one value per freq.
        if controlparameters.docalctf == 1 || controlparameters.docalctf_BEM == 1      
            if n2 == 1
                Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(:,ones(1,nfrequencies));
            else
                if n2 ~= nfrequencies
                    error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one column, or one column per frequency, but it had the size ',int2str(n1),' by ',int2str(n2)])
                end
            end
        else
            if n2 > 1
               error(['ERROR: The Sinputdata.sourceamplitudes input parameter must, for TD calculations, have one column, but it had the size ',int2str(n1),' by ',int2str(n2)])
            end
        end
    end
else % Here we know that there is just one row. Check if there is one value per freq. (or a single value)
    if n2 > 1
        if controlparameters.docalctf == 1 || controlparameters.docalctf_BEM == 1      
            if n2 ~= nfrequencies
                error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have one column, or one column per frequency, but it had the size ',int2str(n1),' by ',int2str(n2)])
            else % Here we know that n2 == nfrequencies. Expand to number of sources
                Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1),:);
            end
        else
            error(['ERROR: The Sinputdata.sourceamplitudes input parameter must, for TD calculations, have one column, but it had the size ',int2str(n1),' by ',int2str(n2)])            
        end
    else % Here we know that n1 = 1 and n2 = 1
        if controlparameters.docalctf == 1 || controlparameters.docalctf_BEM == 1      
            Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1),ones(1,nfrequencies));
        else
            Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1));            
        end
    end
end

if controlparameters.docalcir == 1
    if ~isfield(controlparameters,'savealldifforders')
        controlparameters.savealldifforders = 0;
    end
    if ~isfield(controlparameters,'saveindividualfirstdiff')
        controlparameters.saveindividualfirstdiff = 0;
    end
else
    controlparameters.savealldifforders = 0;
    controlparameters.saveindividualfirstdiff = 0;
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct filehandlingparameters
    
if exist([filehandlingparameters.outputdirectory],'dir') ~=7
      mkdir([filehandlingparameters.outputdirectory])
end

if ~isfield(filehandlingparameters,'filestem')
    if exist('CADfilestem','var') == 1
        ifile = 1;
        filename1 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_tf.mat'];
        filename2 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_ir.mat'];
        if exist(filename1,'file') == 2 || exist(filename2,'file') == 2
            failedtofindfile = 0;
            while (failedtofindfile == 0)
                ifile = ifile + 1;
                filename1 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_tf.mat'];
                filename2 = [filehandlingparameters.outputdirectory,filesep,'results',filesep,CADfilestem,'_',int2str(ifile),'_ir.mat'];
                if exist(filename1,'file') ~= 2 && exist(filename2,'file') ~= 2
                    failedtofindfile = 1;
                end
            end
        end
        filehandlingparameters.filestem = [CADfilestem,'_',int2str(ifile)];
    else
       filehandlingparameters.filestem = Filestem;        
    end
end
% if ~isfield(filehandlingparameters,'savesetupfile')
%     filehandlingparameters.savesetupfile = 1;
% end
if ~isfield(filehandlingparameters,'suppressresultrecycling')
    filehandlingparameters.suppressresultrecycling = 0;
end
if ~isfield(filehandlingparameters,'showtext')
    filehandlingparameters.showtext = 1;
end
if ~isfield(filehandlingparameters,'savecadgeofile')
    filehandlingparameters.savecadgeofile = 0;
end
if ~isfield(filehandlingparameters,'saveSRdatafiles')
    filehandlingparameters.saveSRdatafiles = 1;
end
if ~isfield(filehandlingparameters,'saveeddatafile')
    filehandlingparameters.saveeddatafile = 1;
end
if ~isfield(filehandlingparameters,'saveed2datafile')
    filehandlingparameters.saveed2datafile = 1;
end
if ~isfield(filehandlingparameters,'savesubmatrixdata')
    filehandlingparameters.savesubmatrixdata = 1;
end
if ~isfield(filehandlingparameters,'saveinteqsousigs')
    filehandlingparameters.saveinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'loadinteqsousigs')
    filehandlingparameters.loadinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'savepathsfile')
    filehandlingparameters.savepathsfile = 1;
end

if controlparameters.docalcir == 1
    if ~isfield(filehandlingparameters,'saveISEStree')
        filehandlingparameters.saveISEStree = 0;
    end
    if ~isfield(filehandlingparameters,'savehodpaths')
        filehandlingparameters.savehodpaths = 0;
    end
else
    filehandlingparameters.savehodpaths = 0;
    filehandlingparameters.saveISEStree = 0;
end
if ~isfield(filehandlingparameters,'savediff2result')
    filehandlingparameters.savediff2result = 0;
end
if ~isfield(filehandlingparameters,'savelogfile')
    filehandlingparameters.savelogfile = 1;
end


    

    







