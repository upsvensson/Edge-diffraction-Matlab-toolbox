function [geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters] = EDcheckinputstructs(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters,EDmaincase)
% EDcheckinputstructs checks the input data structs to the various EDmain
% versions and sets default values.
% 
% Input parameters:
%   geofiledata         .geoinputfile        (specified or file open
%                                             window)
%                       .firstcornertoskip   (default: 1e6)
%   Sindata             .coordinates         (obligatory)
%                       .doaddsources        (default: 0 = no)
%   Rindata             .coordinates         (obligatory)
%   envdata             .cair                (default: 344)
%                       .rhoair              (default: 1.21)
%   controlparameters   .fs                  (default: 44100)
%                       .directsound         (default: 1 = yes)
%                       .difforder           (default: 15)
%                       .nedgepoints_visibility (default: 2)
%                       .docalctf            (default: 1)
%                       .docalcir            (default: 0)
%                       .Rstart              (default: 0)
%                       .frequencies         (obligatory)
%                       .discretizationtype  (default: 2 = G-L)
%                       .ngauss              (default: 16)
%   filehandlingparameters    .outputdirectory  (default: /result, in the folder of the geoinputfile)  
%                       .filestem        (default: name of the cad-file, with an underscore + running integer)
%                       .savesetupfile       (default: 1)
%                       .showtext        (default: 1)
%                       .savecadgeofile      (default: 0)
%                       .saveSRindatafiles     (default: 0)
%                       .saveeddatafile      (default: 0)
%                       .savesubmatrixdata   (default: 0)
%                       .saveinteqsousigs     (default: 0)
%                       .loadinteqsousigs     (default: 0)
%                       .savediff2result      (default:0)
%                       .logfilename         (default: '')
%   EDmaincase          1, if convexESIE
% 
% Peter Svensson 28 Nov. 2017 (peter.svensson@ntnu.no)
% 
% [geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters,EDmaincase);

% 24 Nov. 2017 First version
% 28 Nov. 2017 Cleaned code a bit

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
    geofiledata.geoinputfile = CADfile;

end
if ~isfield(geofiledata,'geoinputfile')
	[CADfile,CADfilepath] = uigetfile('*.*','Please select the cadfile');
    [~,CADfile,~] = fileparts([CADfilepath,CADfile]);

    CADfile = [CADfilepath,CADfile];
    geofiledata.geoinputfile = CADfile;

end
if ~isfield(geofiledata,'firstcornertoskip')
    geofiledata.firstcornertoskip = 1e6;
end

[infilepath,CADfilestem] = fileparts(geofiledata.geoinputfile);
if exist([infilepath,filesep,'results'],'dir') ~=7
      mkdir([infilepath,filesep,'results'])
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Rindata

if ~isstruct(Rindata)
    error('ERROR: receiver coordinates were not specified')
end
if ~isfield(Rindata,'coordinates')
    error('ERROR: receiver coordinates were not specified')
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

if EDmaincase == 1
    if ~isstruct(controlparameters)
         error('ERROR: the struct controlparameters was not specified')           
    end
    if ~isfield(controlparameters,'fs')
        controlparameters.fs = 44100;
    end
    if ~isfield(controlparameters,'directsound')
        controlparameters.directsound = 1;
    end
%     if ~isfield(controlparameters,'specorder')
%         controlparameters.specorder = 1;
%     end
    if ~isfield(controlparameters,'difforder')
        controlparameters.difforder = 15;
    end
    if ~isfield(controlparameters,'nedgepoints_visibility')
        controlparameters.nedgepoints_visibility = 2;
    end
    if ~isfield(controlparameters,'docalctf')
        controlparameters.docalctf = 1;
    end
    if ~isfield(controlparameters,'docalcir')
        controlparameters.docalctf = 0;
    end
    if ~isfield(controlparameters,'Rstart')
        controlparameters.Rstart = 0;
    end
    if isfield(controlparameters,'frequencies') == 0 && controlparameters.docalctf == 1
       disp('ERROR: The frequencies were not specified')
       return
    end
    if ~isfield(controlparameters,'discretizationtype')
        controlparameters.discretizationtype = 2;
    end
    if ~isfield(controlparameters,'ngauss')
        controlparameters.ngauss = 16;
    end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct filehandlingparameters
    
if ~isfield(filehandlingparameters,'outputdirectory')
    filehandlingparameters.outputdirectory = infilepath;
end
if ~isfield(filehandlingparameters,'filestem')
    
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
if ~isfield(filehandlingparameters,'saveSRindatafiles')
    filehandlingparameters.saveSRindatafiles = 0;
end
if ~isfield(filehandlingparameters,'saveeddatafile')
    filehandlingparameters.saveeddatafile = 0;
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
if ~isfield(filehandlingparameters,'savediff2result')
    filehandlingparameters.savediff2result = 0;
end
if ~isfield(filehandlingparameters,'logfilename')
    filehandlingparameters.logfilename = '';
end
compstr = computer;
compstr = lower(compstr(1:3));
if compstr == 'mac'  
	filehandlingparameters.lineending = 13;
elseif compstr == 'sun' || compstr == 'sol'            
	filehandlingparameters.lineending = 10;
elseif compstr == 'pcw'
	filehandlingparameters.lineending = [13,10];
else
    error('ERROR: Not implemented for this computer type yet')	
end


    

    







