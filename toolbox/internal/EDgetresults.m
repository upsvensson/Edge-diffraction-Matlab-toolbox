function EDres = EDgetresults(inpar1,inpar2,geoinputdata,Rindata,Sindata,envdata,EDversionnumber)
% EDgetresults reads the relevant files in a directory and returns the data
% as a struct.
%
% Input parameters
%   inpar1          For version 1:  outputdirectory
%                   For version 2: filehandlingparameters
%   inpar2          For version 1: filestem
%                   For version 2: controlparameters
%   geofiledata
%   Rindata
%   Sindata
%   envdata
%   EDversionnumber
%   nfft            NOT IMPL. YET (optional) the fft size to use for the conversion of
%                   impulse responses to tfXXX_fft. Default is the shortest
%                   possible
% 
% Output parameters
%   EDres   A struct with these fields:
%               irdirect
%               irgeom
%               irdiff
%               irhod
%               irtot
%               irdifforders  NOT IMPL. YET  
%               tfdirect
%               tfgeom
%               tfdiff
%               tftot
%               tftot_ESIEBEM
%               tfdirect_fft  NOT IMPL. YET 
%               tfgeom_fft  NOT IMPL. YET 
%               tfdiff_fft  NOT IMPL. YET 
%               tftot_fft   NOT IMPL. YET 
%               fvec_fft    NOT IMPL. YET 
%               nfft        NOT IMPL. YET 
%               geoinputdata                Input data 
%               Sinputdata                  Input data
%               Rinputdata                  Input data
%               envdata                     Input data
%               controlparameters           Input data
%               filehandlingsparameters     Input data
%               EDversionnumber 
% 
% Peter Svensson 16 Apr. 2024   (peter.svensson@ntnu.no)
% 
% EDres =
% EDgetresults(inpar1,inpar2,geofiledata,Rindata,Sindata,envdata,EDversionnumber);

% 1 July 2021 First work
% 16 Sep. 2021 A few more output fields added
% 29 Sep. 2023 Fixed a bug with the size of the ir matrices
% 3 Oct 2023 Made a version 2 with more input parameters
% 16 Apr. 2024 Added the input data structs to EDres

if nargin == 2   % Version 1 of this function
    functionversion = 1;
    outputdirectory = inpar1;
    filestem = inpar2;
else % nargin must be 7 -> Version 2 of this function
    functionversion = 2;
    filehandlingparameters = inpar1;
    controlparameters = inpar2;
    outputdirectory = filehandlingparameters.outputdirectory;
    filestem = filehandlingparameters.filestem;
end

if outputdirectory(end) ~= filesep
    outputdirectory = [outputdirectory,filesep];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a hash for the settings

% We have a special case for the Sindata.sourceamplitudes: if this field
% was not specified in the settings script, then the field was generated
% with default values in the EDcheckinputstructs. Therefore, if the field
% is missing here, the field is generated like in EDcheckinputstructs.

if ~isfield(Sindata,'sourceamplitudes')
    nsources = size(Sindata.coordinates,1);
    Sindata.sourceamplitudes = 1;    
    if controlparameters.docalctf == 1 || controlparameters.docalctf_ESIEBEM == 1      
        nfrequencies = length(controlparameters.frequencies);
        Sindata.sourceamplitudes = Sindata.sourceamplitudes(ones(nsources,1),ones(1,nfrequencies));
    else
        Sindata.sourceamplitudes = Sindata.sourceamplitudes(ones(nsources,1));            
    end
end

% Another special case: controlparameters.surfacegaussorder might not have
% been set in the settings script, so if that field is missing here, we
% generate the field with its default value, like in EDcheckinputstructs.

if ~isfield(controlparameters,'surfacegaussorder')
    controlparameters.surfacegaussorder = 5;
end

if functionversion == 2
    if controlparameters.docalcir == 1
        EDsettingsdatastructload = struct('geoinputdata_corners',geoinputdata.corners,...
            'Sindata_coordinates',Sindata.coordinates,...
            'Sindata_sourceamplitudes',Sindata.sourceamplitudes,...
            'Rindata_coordinates',Rindata.coordinates,...
            'fs',controlparameters.fs,...
            'difforder',controlparameters.difforder,...
            'Rstart',controlparameters.Rstart,...
            'cair',envdata.cair,...
            'EDversionnumber',EDversionnumber);
        EDsettingshashtolookfor = DataHash(EDsettingsdatastructload);        
    end
    if controlparameters.docalctf == 1
        EDsettingsdatastruct = struct('geoinputdata_corners',geoinputdata.corners,...
            'Sindata_coordinates',Sindata.coordinates,...
            'Sindata_sourceamplitudes',Sindata.sourceamplitudes,...
            'Rindata_coordinates',Rindata.coordinates,...
            'frequencies',controlparameters.frequencies,...
            'ngauss',controlparameters.ngauss,...
            'difforder',controlparameters.difforder,...
            'Rstart',controlparameters.Rstart,...
            'cair',envdata.cair,...
            'EDversionnumber',EDversionnumber);
        EDsettingshashtolookfor = DataHash(EDsettingsdatastruct);
    end
    if controlparameters.docalctf_ESIEBEM == 1
        EDsettingsdatastruct = struct('geoinputdata_corners',geoinputdata.corners,...
            'Sindata_coordinates',Sindata.coordinates,...
            'Sindata_sourceamplitudes',Sindata.sourceamplitudes,...
            'Rindata_coordinates',Rindata.coordinates,...
            'frequencies',controlparameters.frequencies,...
            'ngauss',controlparameters.ngauss,...
            'surfacegaussorder',controlparameters.surfacegaussorder,...
            'difforder',controlparameters.difforder,...
            'Rstart',controlparameters.Rstart,...
            'cair',envdata.cair,...
            'EDversionnumber',EDversionnumber);
        EDsettingshashtolookfor = DataHash(EDsettingsdatastruct);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
EDres = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The ir files

if controlparameters.docalcir == 1
    filetolookfor = [outputdirectory,filestem,'_ir.mat'];
    if exist(filetolookfor,'file') == 2
        foundmatch = 1;
        if functionversion == 2
           hashinfile = load(filetolookfor,'EDsettingshash');
            if isfield(hashinfile,'EDsettingshash')
                if ~strcmp(hashinfile.EDsettingshash,EDsettingshashtolookfor)
                    foundmatch = 0;
                end  
            else
                foundmatch = 0;
            end        
        end
        if foundmatch == 1
            eval(['load ',filetolookfor])
            EDres.irdirect = full(irdirect);
            EDres.irgeom = full(irgeom);
            EDres.irdiff = full(irdiff);
            irtot = EDres.irdirect + EDres.irgeom + EDres.irdiff;
            EDres.irtot = irtot;
            [nold,nrec,nsou] = size(EDres.irtot);
        end
    end    
    
    filetolookfor = [outputdirectory,filestem,'_irhod.mat'];
    if exist(filetolookfor,'file') == 2
        foundmatch = 1;
        if functionversion == 2
           hashinfile = load(filetolookfor,'EDsettingshash');
            if isfield(hashinfile,'EDsettingshash')
                if ~strcmp(hashinfile.EDsettingshash,EDsettingshashtolookfor)
                    foundmatch = 0;
                end 
            else
                foundmatch = 0;
            end        
        end
        if foundmatch == 1
            eval(['load ',filetolookfor])
            if iscell(irhod)
                ncells = length(irhod);
                irhodsum = full(irhod{1});
                for ii = 2:ncells
                    irhodsum = irhodsum + full(irhod{ii});
                end
            else
                irhodsum = full(irhod);
            end
            nnew = size(irhodsum,1);
            if nnew > nold
                irdirect = [full(irdirect);zeros(nnew-nold,nrec,nsou)];
                irgeom   = [full(irgeom);zeros(nnew-nold,nrec,nsou)];
                irdiff   = [full(irdiff);zeros(nnew-nold,nrec,nsou)];
                irtot    = [full(irtot);zeros(nnew-nold,nrec,nsou)];
            elseif nnew < nold
                irhodsum = [irhodsum;zeros(nold-nnew,nrec,nsou)];
            end
            irtot = irtot + irhodsum;
            EDres.irdirect = irdirect;
            EDres.irgeom   = irgeom;
            EDres.irdiff   = irdiff;
            EDres.irhod    = irhodsum;
            EDres.irtot    = irtot;
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The tf files

if controlparameters.docalctf == 1
    filetolookfor = [outputdirectory,filestem,'_tf.mat'];
    if exist(filetolookfor,'file') == 2
        foundmatch = 1;
        if functionversion == 2
           hashinfile = load(filetolookfor,'EDsettingshash');
            if isfield(hashinfile,'EDsettingshash')
                if ~strcmp(hashinfile.EDsettingshash,EDsettingshashtolookfor)
                    foundmatch = 0;
                end  
            else
                foundmatch = 0;
            end        
        end
        if foundmatch == 1
            eval(['load ',filetolookfor])
            EDres.tfdirect = full(tfdirect);
            EDres.tfgeom = full(tfgeom);
            EDres.tfdiff = full(tfdiff);
            EDres.tftot = full(tfdirect) + full(tfgeom) + full(tfdiff);
        end
    end    
    
    filetolookfor = [outputdirectory,filestem,'_tfinteq.mat'];
    if exist(filetolookfor,'file') == 2
        foundmatch = 1;
        if functionversion == 2
           hashinfile = load(filetolookfor,'EDsettingshash')
            if isfield(hashinfile,'EDsettingshash')
                if ~strcmp(hashinfile.EDsettingshash,EDsettingshashtolookfor)
                    foundmatch = 0;
                end    
            else
                foundmatch = 0;
            end        
        end
        if foundmatch == 1
            eval(['load ',filetolookfor])
            EDres.tfinteqdiff = full(tfinteqdiff);
            EDres.tftot = EDres.tftot + full(tfinteqdiff);
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The tf_ESIEBEM file

filetolookfor = [outputdirectory,filestem,'_tfESIEBEM.mat'];
if exist(filetolookfor,'file') == 2
    foundmatch = 1;
    if functionversion == 2
       hashinfile = load(filetolookfor,'EDsettingshash')
        if isfield(hashinfile,'EDsettingshash')
            if ~strcmp(hashinfile.EDsettingshash,EDsettingshashtolookfor)
                foundmatch = 0;
            end    
        else
            foundmatch = 0;
        end        
    end
    if foundmatch == 1
        eval(['load ',filetolookfor])
        EDres.tftot_ESIEBEM = full(tftot_ESIEBEM);
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The settings file

filetolookfor = [outputdirectory,filestem,'_settings.mat'];
if exist(filetolookfor,'file') == 2
   eval(['load ',filetolookfor])
   EDres.geoinputdata = geoinputdata;
   EDres.Sinputdata = Sinputdata;
   EDres.Rinputdata = Rinputdata;
   EDres.envdata = envdata;
   EDres.controlparameters = controlparameters;
   EDres.filehandlingparameters = filehandlingparameters;
   EDres.EDversionnumber = EDversionnumber;   
end 






