function EDres = EDassembleresults(resultparts,controlparameters,geoinputdata)
% EDassembleresults puts together the result fields of the results struct
% EDres from the individual sound field components: it sums the components
% to the total field and, for impulse responses, gives the components the
% same length.
%
% This function is used by EDmain_convex (with the components that were
% just computed) and by EDgetresults (with the components loaded from the
% _tf, _tfinteq, _ir and _irhod files), so that the two give the same
% result.
%
% Input parameters:
%   resultparts         A struct with the sound field components. The
%                       fields that are relevant depend on
%                       controlparameters.generalmethod:
%                       'tf':   .tfdirect, .tfgeom, .tfdiff, .tfinteqdiff
%                               (.tfinteqdiff can be empty; then it is
%                               replaced by zeros)
%                       'ir':   .irdirect, .irgeom, .irdiff, .irhod
%                               (.irhod can be a cell variable, as given by
%                               EDmakeHODirs)
%                       'tf_ESIEBEM':   .tftot_ESIEBEM
%                       Fields that are not relevant can be left out.
%   controlparameters   Struct; the fields .generalmethod and
%                       .savealldifforders are used.
%   geoinputdata        Struct; the field .freefieldcase is used.
%
% Output parameters:
%   EDres               A struct with the fields
%                       irdirect, irgeom, irdiff, irhod, irtot,
%                       tfdirect, tfgeom, tfdiff, tfinteqdiff, tftot,
%                       tftot_ESIEBEM
%                       The fields that are not relevant for the
%                       generalmethod are empty.
%
% Peter Svensson (peter.svensson@ntnu.no) 12 Oct. 2026
%
% EDres = EDassembleresults(resultparts,controlparameters,geoinputdata);

% 12 Oct. 2026 First version. The code was moved here from EDmain_convex,
% by Claude, prompted by PS.

% Fields that were not given are set to empty

partnames = {'tfdirect','tfgeom','tfdiff','tfinteqdiff','tftot_ESIEBEM','irdirect','irgeom','irdiff','irhod'};
for ii = 1:length(partnames)
    if ~isfield(resultparts,partnames{ii})
        resultparts.(partnames{ii}) = [];
    end
end

EDres = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulse responses

irdirect = resultparts.irdirect;
irgeom   = resultparts.irgeom;
irdiff   = resultparts.irdiff;
irhod    = resultparts.irhod;

if strcmp(controlparameters.generalmethod,'ir') == 1
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer functions

if strcmp(controlparameters.generalmethod,'tf') == 1
    tfinteqdiff = resultparts.tfinteqdiff;
    if isempty(tfinteqdiff)
        tfinteqdiff = zeros(size(resultparts.tfdirect));
    end
    EDres.tfdirect = resultparts.tfdirect;
    EDres.tfgeom = resultparts.tfgeom;
    EDres.tfdiff = resultparts.tfdiff;
    EDres.tfinteqdiff = tfinteqdiff;
    EDres.tftot = resultparts.tfdirect + resultparts.tfgeom + resultparts.tfdiff + tfinteqdiff;
else
    EDres.tfdirect = [];
    EDres.tfgeom = [];
    EDres.tfdiff = [];
    EDres.tfinteqdiff = [];
    EDres.tftot = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer functions with the ESIEBEM method

if strcmp(controlparameters.generalmethod,'tf_ESIEBEM') == 1
    EDres.tftot_ESIEBEM = resultparts.tftot_ESIEBEM;
else
    EDres.tftot_ESIEBEM = [];
end
