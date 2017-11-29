function numvect = EDextrnums(Str,delim)
% EDextrnums - Extracts the numerical values in a text-string, separated by a given delimiter.
% Periods (.) are interpreted as decimal points. 
%
% Input parameters:
%   Str         Text string, [nchars,1].
%   delim       The delimiter to look for, e.g., comma, space, slash etc.
%               Can be specified as the ASCII value or as the character:
%               '/'
%
% Output parameters:
%   numvect     A list, [1,nnums], of the numerical values that were found.
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
% Peter Svensson (peter.svensson@ntnu.no) 24 Nov. 2017
%
% numvect = EDextrnums(Str,delim);

% 24 Nov. 2017 Copied, without changes, from ESIE2extrnums

if isempty(Str)
	numvect = [];
	return
end
if nargin < 2
	iv = find( (Str>=48 & Str<=57) | (Str==46) | (Str==45) | (Str==101) | (Str==69)  == 1);
else
	iv = find( Str ~= delim );
end
diffiv = diff(iv);

if ~isempty(diffiv)
	startindices = find(diffiv ~= 1);
else
	startindices = [];
end

if length(startindices) >= 1
	startpos = [iv(1) iv(startindices+1)];

	numvect = zeros(1,length(startpos));
	for ii = 1:length(startpos)-1
		shortstr = Str(startpos(ii):startpos(ii+1)-1);
		if nargin < 2
			shortstr = shortstr( find( (shortstr>=48 & shortstr<=57) | (shortstr==46) | (shortstr==45) | (shortstr==101) | (shortstr==69) == 1));
		else
			shortstr = shortstr( find (shortstr ~= delim) );		
        end
 		numvect(ii) = str2double(shortstr);
	end
	shortstr = Str(startpos(ii+1):length(Str));
	if nargin < 2
		shortstr = shortstr( find( (shortstr>=48 & shortstr<=57) | (shortstr==46) | (shortstr==45) | (shortstr==101) | (shortstr==69) == 1));
	else
		shortstr = shortstr( find (shortstr ~= delim) );		
	end
	numvect(ii+1) = str2double(shortstr);
elseif length(iv) >= 1
	shortstr = Str(iv(1):length(Str));
	if nargin < 2
		shortstr = shortstr( find( (shortstr>=48 & shortstr<=57) | (shortstr==46) | (shortstr==45) | (shortstr==101) | (shortstr==69) == 1));
	else
		shortstr = shortstr( find (shortstr ~= delim) );		
    end
	numvect = str2double(shortstr);
else
	numvect = [];
end

