function EDinputdoc(structname)
% EDinputdoc - Describes the fields of the input structs of EDmain_convex.
% For each field: the default value (or "required"), and a description,
% which also says when the field is used. The text is made from the same
% tables that EDcheckinputstructs uses to check the input, so it always
% describes the present version of the toolbox.
%
% Input parameters:
%   structname      (optional) One of 'geoinputdata', 'Sinputdata',
%                   'Rinputdata', 'envdata', 'controlparameters',
%                   'filehandlingparameters'. If it is not given, all six
%                   structs are described.
%
% Output parameters:
%   (none) The text is shown on the screen.
%
% See also EDmain_convex, EDcheckinputstructs, EDdatastructs
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% EDinputdoc(structname);

% 15 Oct. 2026 First version (by Claude, prompted by PS)

specs = EDcheckinputstructs('paramspecs');
allnames = {'geoinputdata','Sinputdata','Rinputdata','envdata','controlparameters','filehandlingparameters'};
if nargin < 1 || isempty(structname)
    structnames = allnames;
else
    if ~ismember(structname,allnames)
        error('EDinputdoc:unknownstruct','%s is not one of the input structs: %s',structname,strjoin(allnames,', '))
    end
    structnames = {structname};
end

textwidth = 76;
indent = repmat(' ',1,8);

for ii = 1:numel(structnames)
    name = structnames{ii};
    fprintf('\n%s\n%s\n',name,repmat('-',1,length(name)));
    textnote = '';
    if isfield(specs.notes,[name,'_text'])
        textnote = specs.notes.([name,'_text']){1};
    end
    if ~isempty(textnote) && ~strcmp(name,'Sinputdata')
        edprintwrapped(textnote,'',textwidth);
    end
    table = specs.(name);
    for jj = 1:size(table,1)
        fieldname = table{jj,1};
        defaultvalue = table{jj,3};
        isrequired = table{jj,5};
        if isrequired
            defaulttext = 'required';
        else
            defaulttext = ['default: ',edvaluetext(defaultvalue)];
        end
        fprintf('  .%s   (%s)\n',fieldname,defaulttext);
        edprintwrapped(table{jj,6},indent,textwidth);
    end
    if isfield(specs.notes,name)
        notes = specs.notes.(name);
        for jj = 1:size(notes,1)
            fprintf('  .%s\n',notes{jj,1});
            edprintwrapped(notes{jj,2},indent,textwidth);
        end
    end
    if ~isempty(textnote) && strcmp(name,'Sinputdata')
        edprintwrapped(textnote,'',textwidth);
    end
end
fprintf('\n');

end

%--------------------------------------------------------------------------
function textstring = edvaluetext(value)
% A short text for a default value.

if isa(value,'function_handle')
    textstring = 'computed, see below';
elseif ischar(value)
    textstring = ['''',value,''''];
elseif isempty(value)
    textstring = '[]';
elseif isnumeric(value) && isscalar(value)
    textstring = num2str(value);
else
    textstring = mat2str(value);
end

end

%--------------------------------------------------------------------------
function edprintwrapped(textstring,indent,textwidth)
% Prints a text, wrapped at word boundaries, with an indent.

words = strsplit(textstring,' ');
oneline = indent;
for ii = 1:numel(words)
    if length(oneline) + length(words{ii}) + 1 > textwidth && length(oneline) > length(indent)
        fprintf('%s\n',oneline);
        oneline = indent;
    end
    if length(oneline) > length(indent)
        oneline = [oneline,' ',words{ii}]; %#ok<AGROW>
    else
        oneline = [oneline,words{ii}]; %#ok<AGROW>
    end
end
if length(oneline) > length(indent)
    fprintf('%s\n',oneline);
end

end
