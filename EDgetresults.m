function EDres = EDgetresults(outputdirectory,filestem)
% EDgetresults reads the relevant files in a directoty and returns the data
% as a struct. PRELIMINARY VERSION
%
% Input parameters
%   outputdirectory 
%   filestem
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
%               tfdirect_fft  NOT IMPL. YET 
%               tfgeom_fft  NOT IMPL. YET 
%               tfdiff_fft  NOT IMPL. YET 
%               tftot_fft   NOT IMPL. YET 
%               fvec_fft    NOT IMPL. YET 
%               nfft        NOT IMPL. YET 
%               controlparameters   A struct NOT IMPL. YET 
% 
% Peter Svensson 16 Sep. 2021   (peter.svensson@ntnu.no)
% 
% EDres = EDgetresults(outputdirectory,filestem);

% 1 July 2021 First work
% 16 Sep. 2021 A few more output fields added

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
    EDres.irtot = EDres.irdirect + EDres.irgeom + EDres.irdiff;
end    
filetolookfor = [outputdirectory,filestem,'_irhod.mat'];
if exist(filetolookfor,'file') == 2
    eval(['load ',filetolookfor])
    EDres.irhod = full(irhod);
    nnew = length(cell2mat(EDres.irhod));
    nold = length(EDres.irtot);
    EDres.irtot = [EDres.irtot;zeros(nnew-nold,1)];
    EDres.irtot = EDres.irtot + cell2mat(EDres.irhod);
    if nnew > nold
       EDres.irdirect = [EDres.irdirect;zeros(nnew-nold,1)];
       EDres.irgeom   = [EDres.irgeom;zeros(nnew-nold,1)];
       EDres.irdiff   = [EDres.irdiff;zeros(nnew-nold,1)];
    end
end    

filetolookfor = [outputdirectory,filestem,'_tf.mat'];
if exist(filetolookfor,'file') == 2
    eval(['load ',filetolookfor])
    EDres.tfdirect = full(tfdirect);
    EDres.tfgeom = full(tfgeom);
    EDres.tfdiff = full(tfdiff);
end    

filetolookfor = [outputdirectory,filestem,'_tf.mat'];
if exist(filetolookfor,'file') == 2
    eval(['load ',filetolookfor])
    EDres.tfdirect = full(tfdirect);
    EDres.tfgeom = full(tfgeom);
    EDres.tfdiff = full(tfdiff);
end    

filetolookfor = [outputdirectory,filestem,'_tfinteq.mat'];
if exist(filetolookfor,'file') == 2
    eval(['load ',filetolookfor])
    EDres.tfinteqdiff = full(tfinteqdiff);
end    





