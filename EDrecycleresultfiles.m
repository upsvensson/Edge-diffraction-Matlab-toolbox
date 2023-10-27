function [foundmatch,existingfilename] = EDrecycleresultfiles(filedirectory,fileending,hashtocomparewith)
% EDrecycleresultfiles looks through all files with a certain fileending in
% the specified directory. From each file, only the hash named
% EDinputdatahash is loaded and it is compared with a reference hash.
% If there is a match this function returns the match.
%
% Input parameters:
%   filedirectory
%   fileending
%   hashtocomparewith
% 
% Output parameters:
%   foundmatch
%   existingfilename
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Oct. 2023
% 
% [foundmatch,existingfilename] =
% EDrecycleresultfiles(outputdirectory,fileending,hashtocomparewith);

% 8 Feb 2018 First version
% 27 Oct. 2023 Fixed a mistake: the 'existingfilename' returned the last
% checked file, whether foundmatch = 1 or 0. After the fix, 
% existingfilename is empty if foundmatch = 0.

listoffilestocheck = dir([filedirectory,filesep,'*',fileending,'.mat']);

nfiles = size(listoffilestocheck,1);
foundmatch = 0;
existingfilename = [];

warning('off','MATLAB:load:variableNotFound');

if nfiles > 0
    ii = 0;
    while foundmatch == 0 && ii < nfiles
        ii = ii + 1;
        existingfilename = [filedirectory,filesep,listoffilestocheck(ii).name];
        hashinfile = load(existingfilename,'EDinputdatahash'); 
        if isfield(hashinfile,'EDinputdatahash')
            if strcmp(hashinfile.EDinputdatahash,hashtocomparewith)
                foundmatch = 1;
            end    
        end
    end
end
if foundmatch == 0
    existingfilename = [];
end

