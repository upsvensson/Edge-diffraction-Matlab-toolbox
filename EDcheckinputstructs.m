function [geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters,EDmaincase)
% EDcheckinputstructs checks the input data structs to the various EDmain
% versions and sets default values.
% 
% Input parameters:
%   geofiledata         .geoinputfile        (specified or file open
%                                             window - unless the fields .corners and .planecorners are given)
%                       .corners             (optional; alternative to
%                                             .geoinputfile. But, the use of this option
%                                             requires that filehandlingparameters.outputdirectory and .filestem are specified)
%                       .planecorners        (optional; alternative to
%                                             .geoinputfile)
%                       .firstcornertoskip   (default: 1e6)
%   Sindata             .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%                       .sourceamplitudes    Only used if doaddsources = 1
%                                            (default:
%                                             ones(nsources,nfrequencies))
%   Rindata             .coordinates         (obligatory)
%   envdata             .cair                (default: 344)
%                       .rhoair              (default: 1.21)
%   controlparameters   .fs                  (default: 44100) Irrelevant for
%                                            EDmain_convexESIE
%                       .directsound         (default: 1 = yes)
%                       .difforder           (default: 15)
%                       .nedgepoints_visibility (default: 2) Irrelevant for
%                                            EDmain_convexESIE
%                       .docalctf            (default: 1)
%                       .docalcir            (default: 0) Irrelevant for
%                                            EDmain_convexESIE
%                       .Rstart              (default: 0)
%                       .frequencies         (obligatory)
%                       .discretizationtype  (default: 2 = G-L)
%                       .ngauss              (default: 16)
%   filehandlingparameters    .outputdirectory  (default: /result, in the folder of the geoinputfile)  
%                       .filestem        (default: name of the cad-file, with an underscore + running integer)
%                       .savesetupfile       (default: 1)
%                       .savecadgeofile      (default: 0)
%                       .saveSRdatafiles     (default: 1)
%                       .saveeddatafile      (default: 1)
%                       .savesubmatrixdata   (default: 0)
%                       .saveinteqsousigs     (default: 0)
%                       .loadinteqsousigs     (default: 0)
%                       .savepathsfile        (default: 0)
%                       .saveISEStree         (default: 0)
%                       .savelogfile          (default: 1)
%                       .savediff2result      (default: 0)
%                       .showtext             (default: 1)
%   EDmaincase          1, if convexESIE (frequency domain)
%                       2, if convexESIE (time domain)
% 
% Peter Svensson 22 Jan. 2018 (peter.svensson@ntnu.no)
% 
% [geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters,EDmaincase);

% 24 Nov. 2017 First version
% 28 Nov. 2017 Cleaned code a bit
% 29 Nov. 2017 Corrected mistake: saveSRdatafiles was called
%              saveSRindatafiles. Adjusted to the new indata option:
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
% 22 Jan 2018 Removed the field planecornertype in the struct geofiledata
% 22 Jan 2018 Moved some controlparameter fields away from the convex TF
% case: fs, nedgepoints_visibility, docalcir, specorder.

if nargin < 7
    disp('ERROR: the input parameter EDmaincase was not specified')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct geofiledata

if ~isstruct(geofiledata)

	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
    [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

    CADfile = [CADfilepath,CADfile];
    geofiledata = struct ('geoinputfile',CADfile);

    [infilepath,CADfilestem] = fileparts(geofiledata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = infilepath;
    end    
end
if isfield(geofiledata,'geoinputfile')
    [infilepath,CADfilestem] = fileparts(geofiledata.geoinputfile);
    if ~isfield(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = infilepath;
    end    
else
    if ~isfield(geofiledata,'corners') || ~isfield(geofiledata,'planecorners')
    	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
        [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

        CADfile = [CADfilepath,CADfile];
        geofiledata.geoinputfile = CADfile;
        [infilepath,CADfilestem] = fileparts(geofiledata.geoinputfile);        
        if ~isfield(filehandlingparameters,'outputdirectory')
            filehandlingparameters.outputdirectory = infilepath;
        end    
    else
        if isfield(filehandlingparameters,'outputdirectory') == 0 || isfield(filehandlingparameters,'filestem') == 0
            error('ERROR: When you give the geometry input in the form of data matrices, you must specify filehandlingparameters.outputdirectory and .filestem')            
        end
    end
end
if ~isfield(geofiledata,'firstcornertoskip')
    geofiledata.firstcornertoskip = 1e6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Sindata

if ~isstruct(Sindata)
    error('ERROR 1: source coordinates were not specified')
end
if ~isfield(Sindata,'coordinates')
    error('ERROR 2: source coordinates were not specified')
end
if ~isfield(Sindata,'doaddsources')
    Sindata.doaddsources = 0;
end
if ~isfield(Sindata,'sourceamplitudes')
    Sindata.sourceamplitudes = 1;
end
nsources = size(Sindata.coordinates,1);
if nsources == 1
    Sindata.doaddsources = 1;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Rindata

if ~isstruct(Rindata)
    error('ERROR 1: receiver coordinates were not specified')
end
if ~isfield(Rindata,'coordinates')
    error('ERROR 2: receiver coordinates were not specified')
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
if EDmaincase == 1
    if ~isfield(controlparameters,'difforder')
        controlparameters.difforder = 15;
    end
    if ~isfield(controlparameters,'docalctf')
        controlparameters.docalctf = 1;
    end
    if ~isfield(controlparameters,'Rstart')
        controlparameters.Rstart = 0;
    end
    if isfield(controlparameters,'frequencies') == 0 
        if controlparameters.docalctf == 1
            error('ERROR: The frequencies were not specified')
        else
           controlparameters.frequencies = []; 
        end
    end
    if ~isfield(controlparameters,'discretizationtype')
        controlparameters.discretizationtype = 2;
    end
    if ~isfield(controlparameters,'ngauss')
        controlparameters.ngauss = 16;
    end
      
    nfrequencies = length(controlparameters.frequencies);
    nsources = size(Sindata.coordinates,1);
    [n1,n2] = size(Sindata.sourceamplitudes);
    if  n1 ~= nsources || n2 ~= nfrequencies
        if n1 == 1 && n2 == 1
            Sindata.sourceamplitudes = ones(nsources,nfrequencies);
        else
            error(['ERROR: The Sindata.sourceamplitudes input parameter must have the size [1,1] or [nsources,nfrequencies], but it had the size ',int2str(n1),' by ',int2str(n2)])
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


    

    







