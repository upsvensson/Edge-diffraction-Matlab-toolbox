function inputhash = EDinputhash(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDinputhash computes one hash value for the complete input to a
% calculation: the five input structs that define the problem, after they
% have been through EDcheckinputstructs (so that all default values are
% filled in).
%
% Two runs get the same hash if, and only if, all the values in the five
% structs are the same (and, if filehandlingparameters is given, the
% field .savediff2result, which changes what the results contain).
% The hash does not depend on
%   - the other fields of filehandlingparameters (where and what to save)
%   - the field .inputorigin (if a value was given by the user or was a
%     default value)
%   - geoinputdata.geoinputfile (the geometry is represented by .corners
%     and .planecorners, so a renamed geometry file gives the same hash)
%   - the order of the fields in the structs
%   - the version of the EDtoolbox
%
% This hash is saved in the results file (_EDres.mat) by EDmain_convex, and
% is used by EDgetresults to find the results file for a given set of
% input data.
%
% Input parameters:
%   geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters
%                   The five structs, as returned by EDcheckinputstructs.
%   filehandlingparameters  (optional) Only the field .savediff2result is
%                   used.
%
% Output parameters:
%   inputhash       A text string.
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% inputhash = EDinputhash(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 11 Oct. 2026 First version, by Claude, prompted by PS
% 14 Oct. 2026 Optional input parameter filehandlingparameters; its field
% .savediff2result is part of the hash.
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

structlist = {geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters};
fieldstoskip = {'inputorigin','geoinputfile'};

for ii = 1:length(structlist)
    onestruct = structlist{ii};
    onestruct = rmfield(onestruct,fieldstoskip(isfield(onestruct,fieldstoskip)));
    structlist{ii} = orderfields(onestruct);
end

hashdata = struct;
hashdata.geoinputdata      = structlist{1};
hashdata.Sinputdata        = structlist{2};
hashdata.Rinputdata        = structlist{3};
hashdata.envdata           = structlist{4};
hashdata.controlparameters = structlist{5};
if nargin >= 6 && isstruct(filehandlingparameters) && isfield(filehandlingparameters,'savediff2result')
    hashdata.savediff2result = filehandlingparameters.savediff2result;
end

inputhash = DataHash(hashdata);
