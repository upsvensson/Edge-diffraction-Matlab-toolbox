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
% Uses no special subroutines
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 27 Nov. 2017
%
% Filenameout = EDstrpend(Filenamein,strptext);

% 27 Nov. 2017 Copied without changes from ESIE2toolbox

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

