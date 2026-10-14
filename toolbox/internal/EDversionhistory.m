function EDversionhistory
% EDversionhistory prints the EDtoolbox version history on the screen.
% The history is kept in the file CHANGELOG.md, in the top folder of the
% EDtoolbox (newest version first). That file can also be read in any
% text editor, or on GitHub.
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
%
% EDversionhistory

% 28 Jan 2018 First version, just after version 0.101
% 14 Oct. 2026 The version history was moved to the file CHANGELOG.md;
% this function now prints that file. (By Claude, prompted by PS.)

% Look for CHANGELOG.md in the folder of this function, and up to three
% folders above it.
searchfolder = fileparts(mfilename('fullpath'));
changelogfile = '';
for ii = 1:4
    candidate = fullfile(searchfolder,'CHANGELOG.md');
    if exist(candidate,'file') == 2
        changelogfile = candidate;
        break
    end
    parentfolder = fileparts(searchfolder);
    if strcmp(parentfolder,searchfolder)
        break
    end
    searchfolder = parentfolder;
end

if isempty(changelogfile)
    disp('EDversionhistory: the file CHANGELOG.md could not be found in the EDtoolbox folders.')
    return
end

type(changelogfile)
