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
% Peter Svensson (peter.svensson@ntnu.no) 8 Feb 2018
% 
% [foundmatch,existingfilename] =
% EDrecycleresultfiles(outputdirectory,fileending,hashtocomparewith);

% 8 Feb 2018 First version

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
