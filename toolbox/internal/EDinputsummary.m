function textlines = EDinputsummary(inputstruct,structname)
% EDinputsummary makes a text summary of one of the input structs, for the
% log file: one line per field, with the value and the origin of the value
% (user/default/derived/not applicable).
%
% Small values are written out in full. For large matrices, the size and
% the range of values are given.
%
% Input parameters:
%   inputstruct     One of the six input structs, after
%                   EDcheckinputstructs, and with the origin record
%                   attached as the field .inputorigin (as in EDres and in
%                   the _settings file). If the field .inputorigin is
%                   missing, the origins are left out.
%   structname      The name of the struct, e.g. 'controlparameters'.
%
% Output parameters:
%   textlines       A cell array with one text string per field, e.g.
%                       'ngauss = 16   (default)'
%                       'frequencies = [1x400], 50 to 20000   (user)'
%                       'fs : not applicable'
%
% Peter Svensson (peter.svensson@ntnu.no) 12 Oct. 2026
%
% textlines = EDinputsummary(inputstruct,structname);

% 12 Oct. 2026 First version, by Claude, prompted by PS

if isfield(inputstruct,'inputorigin') && isstruct(inputstruct.inputorigin)
    origin = inputstruct.inputorigin;
    inputstruct = rmfield(inputstruct,'inputorigin');
else
    origin = struct;
end

names = fieldnames(inputstruct);
nnames = length(names);
textlines = cell(1,nnames);

% Align the equal signs
namewidth = max(cellfun(@length,names));

for ii = 1:nnames
    name = names{ii};
    paddedname = [name,blanks(namewidth-length(name))];

    if isfield(origin,name) && ischar(origin.(name))
        originvalue = origin.(name);
    else
        originvalue = '';
    end

    if strcmp(originvalue,'notapplicable')
        textlines{ii} = [paddedname,' : not applicable'];
    else
        textline = [paddedname,' = ',edvaluetotext(inputstruct.(name))];
        if ~isempty(originvalue)
            textline = [textline,'   (',originvalue]; %#ok<AGROW>
            % For a geometry that was read from a file: tell which file
            if strcmp(originvalue,'derived') && strcmp(structname,'geoinputdata') && ...
                    any(strcmp(name,{'corners','planecorners'})) && ...
                    isfield(inputstruct,'geoinputfile') && ischar(inputstruct.geoinputfile)
                [~,geofilename,geofileext] = fileparts(inputstruct.geoinputfile);
                textline = [textline,' from ',geofilename,geofileext]; %#ok<AGROW>
            end
            textline = [textline,')']; %#ok<AGROW>
        end
        textlines{ii} = textline;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valuetext = edvaluetotext(value)
% Returns a short text which describes a value: the value itself if it is
% small; otherwise the size and (for numbers) the range.

maxnumelinfull = 6;

if ischar(value)
    valuetext = ['''',value,''''];
elseif isnumeric(value) || islogical(value)
    sizetext = ['[',regexprep(int2str(size(value)),'\s+','x'),']'];
    if isempty(value)
        valuetext = '[]';
    elseif isscalar(value)
        valuetext = num2str(double(value),10);
    elseif numel(value) <= maxnumelinfull && ismatrix(value)
        valuetext = mat2str(double(value),6);
    elseif isreal(value)
        valuetext = [sizetext,', ',num2str(double(min(value(:))),6),' to ',num2str(double(max(value(:))),6)];
    else
        valuetext = [sizetext,', complex'];
    end
elseif iscell(value)
    valuetext = ['cell [',regexprep(int2str(size(value)),'\s+','x'),']'];
elseif isstruct(value)
    valuetext = 'struct';
else
    valuetext = class(value);
end
