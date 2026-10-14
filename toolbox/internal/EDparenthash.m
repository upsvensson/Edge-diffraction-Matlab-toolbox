function parenthash = EDparenthash(datastruct)
% EDparenthash returns the stage hash that was stored in a data struct
% (edgedata, Sdata, Rdata, firstorderpathdata, edgetoedgedata) by the
% function that created it. Used for chained hashes: a stage includes
% the hash(es) of the stage(s) it builds on in its own EDinputdatahash,
% so that a change in an earlier stage (including a change of its
% stageversion) invalidates all later stage files.
%
% If the struct has no field EDinputdatahash (e.g., the dummy structs of
% the free-field case), a hash of the complete struct is returned instead.
%
% Input parameters:
%   datastruct      A struct, normally with the field .EDinputdatahash
%
% Output parameters:
%   parenthash      A hash string
%
% Uses the function DataHash from Matlab Central
%
% Peter Svensson (peter.svensson@ntnu.no) 12 Oct. 2026
%
% parenthash = EDparenthash(datastruct);

% 12 Oct. 2026 First version

if isfield(datastruct,'EDinputdatahash')
    parenthash = datastruct.EDinputdatahash;
else
    parenthash = DataHash(datastruct);
end
