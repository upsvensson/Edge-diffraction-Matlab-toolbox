function EDpostfunctext(functype,timingvalue,existingfilename,...
    filehandlingparameters,fid,TTT,varargin)
% EDpostfunctext will print various text on the screen and/or to a file
% after a function has been called by EDmain_xxx functions.
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
%                   that will be printed on one line each.
%
% Uses the function EDmessage
%
% Peter Svensson 27 September 2023
%
% EDpostfunctext(functype,timingvalue,existingfilename,...
%    filehandlingparameters,fid,TTT,varargin)

% 7 March 2023 First version, to remove code from the EDmain_ functions
% 27 Sep. 2023 Polished it, for inclusion in the EDtoolbox

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

commandstring = ['EDmessage(filehandlingparameters,''sf'',fid,1,existingfilename,textline1'];
for ii = 2:ntextstrings
    eval(['textstring',int2str(ii),' = setstr(varargin{ii});'])
    commandstring = [commandstring,',textstring',int2str(ii)];
end

commandstring = [commandstring,');'];
eval(commandstring)
