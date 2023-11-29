function fid = EDmessage(filehandlingparameters,showprint,fid,showtextlevelneeded,...
    existingfilename,varargin)
% EDmessage prints out one or more lines of text on screen and/or
% writes into a log-file.
%
% Input parameters:
%   filehandlingparameters      The struct with parameter settings
%   showprint                   A text string, 's' or 'f' or 'sf'
%                               if 's' is included: show text on screen
%                               if 'f' is included: write to logfile
%   fid                         The file id, if the file has already been opened.
%                               If fid is zero, a file will be opened (if
%                               writing to file is required, i.e., if 'f' is
%                               included in showprint).
%   showtextlevelneeded         A value 0,1,2,3,... The lowest value of
%                               filehandlingparameters.showtext that is
%                               needed for the text to be printed on
%                               screen.
%   existingfilename            A string: '' or a file name. If not empty,
%                               then the text "Recycled and duplicated"
%                               will be printed, followed by the file name.
%   varargin                    One or more text strings, separated with commas,
%                               that will be printed on one line each. If
%                               the last one is an empty string, a blank
%                               line will be printed. If the two last ones
%                               are empty, only one blank will be printed.
%
% Output parameters
%   fid                         The file id: either 0 or the value for the
%                               actually opened file
%
% Peter Svensson 6 Oct. 2023 (peter.svensson@ntnu.no)
%
% fid = EDmessage(filehandlingparameters,showprint,fid,showtextlevelneeded,...
%    existingfilename,varargin);

% 7 Sep. 2022 First version. Made with the purpose to have fewer lines of
% code in the EDmain functions.
% 8 Sep. 2022 Added the foundmatch and existingfilename input parameters
% to move even more code from the calling function.
% 9 Sep. 2022 Removed the foundmatch; existingfilename is enough
% 6 Feb. 2023 Added the optional input parameter addspace
% 8 Aug. 2023 Stopped using the int2bit function since it requires the
% Communications toolbox
% 26 Sep. 2023 Introduced the varargin input parameter so that 1-N lines of
% text can be specified (and transfered as a cell variable).
% 2 Oct. 2023 If the two last input strings are empty, only one will be
% shown.
% 6 Oct. 2023 Fixed error: when a file is recycled, the detailed timing
% data wasnt written to the log file.

ntextstrings = nargin - 5;

if ntextstrings > 0
    textstr1 = setstr(varargin{1});
    if ntextstrings >= 3
        if isempty(varargin{ntextstrings}) && isempty(varargin{ntextstrings-1})
            ntextstrings = ntextstrings - 1;
        end
    end
end

if ispc == 1
   lineending = [13,10];
else
    lineending = 10;
end

if ~isempty( strfind(showprint,'s') )
    showtext = 1;
else
    showtext = 0;
end
if ~isempty( strfind(showprint,'f') )
    writetofile = 1;
else
    writetofile = 0;
end

% if textstr1(1:2) == 'ED'
%     pretext = '####################################################################';
% else
    pretext = '';
%end

if showtext == 1 & filehandlingparameters.showtext >= showtextlevelneeded
    if ~isempty(pretext)
        disp(' ')
        disp(pretext)
    end
    disp(textstr1)
    if ~isempty(existingfilename)
        disp(['   Recycled and duplicated ',existingfilename])
       for ii = 2:ntextstrings
            disp(['   (',setstr(varargin{ii}),')'])
        end
    else
        for ii = 2:ntextstrings
            disp(['   ',setstr(varargin{ii})])
        end
    end
    if ~isempty(pretext)
        disp(' ')
    end
end

if writetofile == 1
    if filehandlingparameters.savelogfile == 1
        logfilename = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_log.txt'];
        if fid == 0
            fid = fopen(logfilename,'w');
            if fid == -1
                error('ERROR: The planned logfile is not possible to open - check that it isn''t opened by any other program!')
    	        return
            end
        end
    
        if ~isempty(pretext)
            fwrite(fid,[pretext,lineending],'char');
        end
        fwrite(fid,[textstr1,lineending],'char');
        if ~isempty(existingfilename)
            fwrite(fid,['   (Recycled and duplicated ',existingfilename,')',lineending],'char');
            for ii = 2:ntextstrings
                fwrite(fid,['   (',setstr(varargin{ii}),')',lineending],'char');
            end
        else
            for ii = 2:ntextstrings
                fwrite(fid,['   ',setstr(varargin{ii}),lineending],'char');
            end
        end
        if ~isempty(pretext)
            fwrite(fid,[' ',lineending],'char');
        end
    end
end
