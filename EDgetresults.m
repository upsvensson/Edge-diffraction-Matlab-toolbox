function EDres = EDgetresults(outputdirectory,filestem)
% EDgetresults reads the relevant files in a directoty and returns the data
% as a struct.
%
% Input parameters
%   outputdirectory 
%   filestem
%   nfft            (optional) the fft size to use for the conversion of
%                   impulse responses to tfXXX_fft. Default is the shortest
%                   possible
% 
% Output parameters
%   EDres   A struct with these fields:
%               irdirect
%               irgeom
%               irdiff
%               irtot
%               irdifforders
%               tfdirect
%               tfgeom
%               tfdiff
%               tftot
%               tfdirect_fft
%               tfgeom_fft
%               tfdiff_fft
%               tftot_fft
%               fvec_fft
%               nfft
%               controlparameters   A struct
% 
% Peter Svensson 1 July 2021   (peter.svensson@ntnu.no)
% 
% EDres = EDgetresults(outputdirectory,filestem);

if outputdirectory(end) ~= filesep
    outputdirectory = [outputdirectory,filesep];
end
    
EDres = struct;

filetolookfor = [outputdirectory,filestem,'_ir.mat'];
if exist(filetolookfor,'file') == 2
    eval(['load ',filetolookfor])
    EDres.irdirect = full(irdirect);
    EDres.irgeom = full(irgeom);
    EDres.irdiff = full(irdiff);
end    



