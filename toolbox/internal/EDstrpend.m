function Filenameout = EDstrpend(Filenamein,striptext)
% EDstrpend - Removes a specified ending from a filename.
% First, any extension is removed.
%
% Input parameters:
%   Filenemain	A text string with the filename
%   stripext    A text string with the ending that should 
%               be removed.
%
% Output parameters:
%   Filenameout  A text string with the extension-stripped filename
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% Filenameout = EDstrpend(Filenamein,strptext);

% 27 Nov. 2017 Copied without changes from ESIE2toolbox
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

%Filenameout = ESIE2strpext(Filenamein);
[Filepath,Filenameout,~] = fileparts(Filenamein);
Filenameout = [Filepath,filesep,Filenameout];

str1 = lower(Filenameout);
str2 = lower(striptext);
n1 = length(str1);
n2 = length(str2);
if n1 >= n2
	if str1(n1-n2+1:n1) == str2
		Filenameout = Filenameout(1:n1-n2);
	end
end

