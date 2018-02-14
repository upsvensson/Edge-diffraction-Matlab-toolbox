function [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,EDmaincase)
% EDcheckinputstructs checks the input data structs to the various EDmain
% versions and sets default values.
% 
% Input parameters:
%   geoinputdata     .geoinputfile        (specified, or a file open
%                                             window is given - unless the fields .corners and .planecorners are given)
%                   .corners             (optional; alternative to
%                                             .geoinputfile. But, the use of this option
%                                             requires that filehandlingparameters.outputdirectory and .filestem are specified)
%                   .planecorners        (optional; alternative to
%                                             .geoinputfile)
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
%   controlparameters   .fs              (default: 44100) Ignored for
%                                        EDmain_convexESIE, but used by
%                                        EDmain_convexESIE_ir
%                   .directsound         (default: 1 = yes)
%                   .specorder           (ignored by EDmain_convexESIE and 
%                                        by EDmain_convexESIE_ir)
%                   .difforder           (default: 15)
%                   .skipfirstorder      (default: 0)
%                   .nedgepoints_visibility (default: 2) Ignored by
%                                        EDmain_convexESIE and by
%                                        EDmain_convexESIE_ir
%                   .docalctf            (default: 1 for EDmain_convexESIE,
%                                                  0 for EDmain_convexESIE_ir)
%                   .docalcir            (default: 0 for EDmain_convexESIE,
%                                                  1 for EDMain_convexESIEtime)
%                   .Rstart              (default: 0)
%                   .frequencies         (obligatory for EDmain_convexESIE,
%                                         ignored by EDmain_convexESIE_ir)
%                   .discretizationtype  (default: 2 = G-L)
%                   .ngauss              (default: 16)
%   filehandlingparameters    .outputdirectory  (default: the folder of the geoinputfile)  
%                   .filestem        (default: name of the cad-file, with an underscore + running integer)
%                   .savesetupfile       (default: 1)
%                   .savecadgeofile      (default: 0)
%                   .saveSRdatafiles     (default: 1)
%                   .saveeddatafile      (default: 1)
%                   .savesubmatrixdata   (default: 0)
%                   .saveinteqsousigs     (default: 0)
%                   .loadinteqsousigs     (default: 0)
%                   .savepathsfile        (default: 0)
%                   .saveISEStree         (default: 0)
%                   .savelogfile          (default: 1)
%                   .savediff2result      (default: 0)
%                   .showtext             (default: 1)
%   EDmaincase      1, for EDmain_convexESIE (frequency domain)
%                   2, for EDmain_convexESIE_ir (time domain)
% 
% Peter Svensson 8 Feb 2018 (peter.svensson@ntnu.no)
% 
% [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,EDmaincase);

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

if nargin < 7
    disp('ERROR: the input parameter EDmaincase was not specified')
    return
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
end
if isfield(geoinputdata,'geoinputfile')
    [infilepath,CADfilestem] = fileparts(geoinputdata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Sinputdata

sourceamplitudes_default = 0;

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
if ~isfield(Sinputdata,'sourceamplitudes')
    Sinputdata.sourceamplitudes = ones(1,nsources);    
    sourceamplitudes_default = 1;
end
if nsources == 1
    Sinputdata.doaddsources = 1;
end    
if ~isfield(Sinputdata,'doallSRcombinations')
    Sinputdata.doallSRcombinations = 1;
else
    if Sinputdata.doallSRcombinations == 0 && nsources~=nreceivers
       error('ERROR: doallSRcombinations was set to 0, but the number of sources was not the same as the number of receivers'); 
    end
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
if ~isfield(controlparameters,'directsound')
    controlparameters.directsound = 1;
end
if ~isfield(controlparameters,'difforder')
    controlparameters.difforder = 15;
end
if ~isfield(controlparameters,'Rstart')
    controlparameters.Rstart = 0;
end
if ~isfield(controlparameters,'discretizationtype')
    controlparameters.discretizationtype = 2;
end
if ~isfield(controlparameters,'ngauss')
    controlparameters.ngauss = 16;
end
if ~isfield(controlparameters,'skipfirstorder')
    controlparameters.skipfirstorder = 0;
end

if EDmaincase == 1
    if ~isfield(controlparameters,'docalctf')
        controlparameters.docalctf = 1;
    end
    if isfield(controlparameters,'frequencies') == 0 
        if controlparameters.docalctf == 1
            error('ERROR: The frequencies were not specified')
        else
           controlparameters.frequencies = []; 
        end
    end
      
    nfrequencies = length(controlparameters.frequencies);
    nsources = size(Sinputdata.coordinates,1);
    [n1,n2] = size(Sinputdata.sourceamplitudes);
    if  n1 ~= nsources || n2 ~= nfrequencies
        if (n1 == 1 && n2 == 1) || sourceamplitudes_default == 1
            Sinputdata.sourceamplitudes = ones(nsources,nfrequencies);
        else
            error(['ERROR: The Sinputdata.sourceamplitudes input parameter must have the size [1,1] or [nsources,nfrequencies], but it had the size ',int2str(n1),' by ',int2str(n2)])
        end
    end   
end

if EDmaincase == 2
    if ~isfield(controlparameters,'fs')
        controlparameters.fs = 44100;
    end    
    if ~isfield(controlparameters,'docalcir')
        controlparameters.docalcir = 1;
    end
end

if EDmaincase > 2
     if ~isfield(controlparameters,'nedgepoints_visibility')
        controlparameters.nedgepoints_visibility = 2;
     end   
    if ~isfield(controlparameters,'specorder')
        controlparameters.specorder = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct filehandlingparameters
    
if exist([filehandlingparameters.outputdirectory,filesep,'results'],'dir') ~=7
      mkdir([filehandlingparameters.outputdirectory,filesep,'results'])
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
if ~isfield(filehandlingparameters,'savesetupfile')
    filehandlingparameters.savesetupfile = 1;
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
if ~isfield(filehandlingparameters,'savesubmatrixdata')
    filehandlingparameters.savesubmatrixdata = 0;
end
if ~isfield(filehandlingparameters,'saveinteqsousigs')
    filehandlingparameters.saveinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'loadinteqsousigs')
    filehandlingparameters.loadinteqsousigs = 0;
end
if ~isfield(filehandlingparameters,'savepathsfile')
    filehandlingparameters.savepathsfile = 0;
end
if ~isfield(filehandlingparameters,'saveISEStree')
    filehandlingparameters.saveISEStree = 0;
end
if ~isfield(filehandlingparameters,'savediff2result')
    filehandlingparameters.savediff2result = 0;
end
if ~isfield(filehandlingparameters,'savelogfile')
    filehandlingparameters.savelogfile = 1;
end


    

    







