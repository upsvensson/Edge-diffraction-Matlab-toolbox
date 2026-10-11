function EDpostfunctext(functype,timingvalue,existingfilename,...
    filehandlingparameters,fid,TTT,varargin)
% EDpostfunctext will print various text on the screen and/or to a file
% after a function has been called by EDmain_xxx functions.
%
% The log file gets all the text. On the screen, only the first line (the
% function name and the time) is shown, unless
% filehandlingparameters.showtext is 2 or higher; then all the text is
% shown on the screen too.
%
% Input parameters
%   functype      A text string which selects which function has been
%                   called: 'EDedgeo' etc
%   timingvalue     A value in seconds which will be printed out with some
%                   text which is specified inside this function. The 
%                   timingvalue might be a list of two values, where the 
%                   second one is the original (in case of reuse).
%   existingfilename 
%   filehandlingparameters
%   fid
%   TTT             A text string of spaces, which indents the text
%   varargin        One or more text strings, separated with commas,
%                   that will be printed on one line each. The first one is
%                   put on the same line as the function name and the time.
%                   The following ones are the detail lines.
%
% Uses the function EDmessage
%
% Peter Svensson 12 Oct. 2026
%
% EDpostfunctext(functype,timingvalue,existingfilename,...
%    filehandlingparameters,fid,TTT,varargin)

% 7 March 2023 First version, to remove code from the EDmain_ functions
% 27 Sep. 2023 Polished it, for inclusion in the EDtoolbox
% 12 Oct. 2026 The detail lines, and the name of a reused file, are written
% to the log file only, unless showtext >= 2 (by Claude, prompted by PS).
% Removed the use of eval.

ntextstrings = nargin - 6;
if ntextstrings > 0
    textstr1 = setstr(varargin{1});
else
    textstr1 = [];
end

Timestring = [' Time: ',num2str(timingvalue(1)),' s'];
if length(timingvalue) > 1
    Timestring = [Timestring,' (Orig.: ',num2str(timingvalue(2)), 's)'];  
end

if isempty(textstr1)
    textline1 = [TTT,functype,': ',Timestring];
else
    textline1 = [TTT,functype,': ',textstr1,Timestring];
end

if ~isempty(existingfilename)
    [~,existingfilename,~] = fileparts(existingfilename); 
end

detaillines = cell(1,max(ntextstrings-1,0));
for ii = 2:ntextstrings
    detaillines{ii-1} = char(varargin{ii});
end

% The log file gets everything: the first line, the name of a reused file,
% and the detail lines.
EDmessage(filehandlingparameters,'f',fid,1,existingfilename,textline1,detaillines{:});

% The screen gets only the first line (function name and time), unless
% showtext is 2 or higher; then the details are shown on the screen too.
if filehandlingparameters.showtext >= 2
    EDmessage(filehandlingparameters,'s',fid,1,existingfilename,textline1,detaillines{:});
else
    EDmessage(filehandlingparameters,'s',fid,1,'',textline1);
end
