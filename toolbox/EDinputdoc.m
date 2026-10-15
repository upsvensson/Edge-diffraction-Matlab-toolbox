function EDinputdoc(structname,latexfile)
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
%                   structs are described. '' or [] also gives all six.
%   latexfile       (optional) The name of a .tex file. If it is given, the
%                   description is written to this file, in LaTeX form,
%                   instead of to the screen. The file has only content,
%                   given with the macros below, which the document
%                   that inputs the file must define (see
%                   EDtoolbox_cheatsheet.tex for an example):
%                     \EDstruct{name}            start of an input struct
%                     \EDstructtext{text}        a general note on a struct
%                     \EDfield{name}{default}{description}
%                     \EDval{value}              a numerical value
%                     \EDstr{text}               a text value, like
%                                                GMRES_fast (without quotes)
%                     \EDfile{name}              a function or file name
%                     \EDrequired                written as the default of
%                                                required fields
%                     \EDderived                 written as the default when
%                                                it is computed from other
%                                                fields
%
% Output parameters:
%   (none) The text is shown on the screen, or written to latexfile.
%
% See also EDmain_convex, EDcheckinputstructs, EDdatastructs
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% EDinputdoc(structname,latexfile);

% 15 Oct. 2026 First version (by Claude, prompted by PS)
% 15 Oct. 2026 New optional input parameter latexfile: the description is
% written as LaTeX, for the cheatsheet and the manual.
% 15 Oct. 2026 No date in the generated LaTeX file.
% 15 Oct. 2026 The patterns for field and function names were rewritten
% without nested groups (Matlab did not set the field names in typewriter
% font). They were then replaced by a scan of the characters
% (edmarknames), since Matlab still did not find the field names.

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

if nargin >= 2 && ~isempty(latexfile)
    edwritelatex(specs,structnames,latexfile);
    return
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

%--------------------------------------------------------------------------
function edwritelatex(specs,structnames,latexfile)
% Writes the description to a .tex file, with the macros that are described
% in the header.

fid = fopen(latexfile,'w');
if fid == -1
    error('EDinputdoc:cannotopen','The file %s could not be opened for writing',latexfile)
end
closefile = onCleanup(@() fclose(fid));

try
    versionnumber = num2str(EDgetversion);
catch
    versionnumber = '?';
end
% No date in the file: it is then unchanged as long as the specification
% is unchanged (EDverify, test 19, compares it with a new one).
fprintf(fid,'%% Generated by EDinputdoc, EDtoolbox version %s.\n',versionnumber);
fprintf(fid,'%% Do not edit this file: change the descriptions in EDcheckinputstructs\n');
fprintf(fid,'%% (the local functions edparamspec and ednotes), and run EDinputdoc again.\n');

for ii = 1:numel(structnames)
    name = structnames{ii};
    fprintf(fid,'\n\\EDstruct{%s}\n',name);
    textnote = '';
    if isfield(specs.notes,[name,'_text'])
        textnote = specs.notes.([name,'_text']){1};
    end
    if ~isempty(textnote) && ~strcmp(name,'Sinputdata')
        fprintf(fid,'\\EDstructtext{%s}\n',edlatextext(textnote));
    end
    table = specs.(name);
    for jj = 1:size(table,1)
        if table{jj,5}
            defaulttext = '\EDrequired';
        else
            defaulttext = edlatexvalue(table{jj,3});
        end
        fprintf(fid,'\\EDfield{%s}{%s}{%s}\n',edlatexescape(table{jj,1}),defaulttext,...
            edlatextext(table{jj,6}));
    end
    if isfield(specs.notes,name)
        notes = specs.notes.(name);
        for jj = 1:size(notes,1)
            % A note text can start with "Default: <value>. ", which then
            % goes to the default column.
            descriptiontext = notes{jj,2};
            defaulttext = '';
            tokens = regexp(descriptiontext,'^Default: (\S+)\.\s+(.*)$','tokens','once');
            if ~isempty(tokens)
                defaulttext = ['\EDval{',edlatexescape(tokens{1}),'}'];
                descriptiontext = tokens{2};
            elseif ~isempty(strfind(descriptiontext,'Required'))
                defaulttext = '\EDrequired';
            end
            fprintf(fid,'\\EDfield{%s}{%s}{%s}\n',edlatexescape(notes{jj,1}),defaulttext,...
                edlatextext(descriptiontext));
        end
    end
    if ~isempty(textnote) && strcmp(name,'Sinputdata')
        fprintf(fid,'\\EDstructtext{%s}\n',edlatextext(textnote));
    end
end

end

%--------------------------------------------------------------------------
function textstring = edlatexvalue(value)
% A default value, in LaTeX form.

if isa(value,'function_handle')
    textstring = '\EDderived';
elseif ischar(value)
    textstring = ['\EDstr{',edlatexescape(value),'}'];
elseif isempty(value)
    textstring = '\EDval{[]}';
elseif isnumeric(value) && isscalar(value)
    if value >= 1e5 && value == 10^round(log10(value))
        % 1000000 -> 1e6
        textstring = ['\EDval{1e',int2str(round(log10(value))),'}'];
    else
        textstring = ['\EDval{',edlatexescape(num2str(value)),'}'];
    end
else
    textstring = ['\EDval{',edlatexescape(mat2str(value)),'}'];
end

end

%--------------------------------------------------------------------------
function textstring = edlatextext(textstring)
% A description text in LaTeX form: the special characters are escaped,
% quoted values ('GMRES_fast') are given as \EDstr{...}, function names
% (EDmain_convex) and stage files (_eddata) as \EDfile{...}, and field
% names (.corners) are set in typewriter font.

textstring = edlatexescape(textstring);
% Quoted values: 'tf' -> \EDstr{tf}
textstring = regexprep(textstring,'''([^''\s]+)''','\\EDstr{$1}');
% Function names (EDmain_convex) and stage files (_eddata), and field
% names (.corners), are found by a scan of the characters (the regular
% expressions that were used first gave different results in Matlab and
% Octave).
textstring = edmarknames(textstring);

end

%--------------------------------------------------------------------------
function textstring = edlatexescape(textstring)
% Escapes the LaTeX special characters.

textstring = strrep(textstring,'\','\textbackslash ');
textstring = regexprep(textstring,'([{}_%&#$])','\\$1');
textstring = strrep(textstring,'\textbackslash ','\textbackslash{}');
textstring = strrep(textstring,'^','\^{}');
textstring = strrep(textstring,'~','\~{}');

end

%--------------------------------------------------------------------------
function outstring = edmarknames(textstring)
% Marks names in an (escaped) description text:
%   .corners        -> \texttt{.corners}  (a field name: a dot, after a
%                      blank, '(' or '/', followed by a letter)
%   EDmain\_convex  -> \EDfile{EDmain\_convex}  (after a blank or '(')
%   \_eddata        -> \EDfile{\_eddata}        (after a blank or '(')
% A name consists of letters and digits, and \_ followed by a letter or
% digit.

n = length(textstring);
outstring = '';
k = 1;
while k <= n
    if k == 1
        previous = ' ';
    else
        previous = textstring(k-1);
    end
    isfield = textstring(k) == '.' && any(previous == ' (/') && ...
        k < n && edisletter(textstring(k+1));
    isfunction = any(previous == ' (') && k+1 <= n && strcmp(textstring(k:k+1),'ED');
    isstagefile = any(previous == ' (') && k+2 <= n && strcmp(textstring(k:k+1),'\_') && ...
        edisletter(textstring(k+2));
    if ~(isfield || isfunction || isstagefile)
        outstring = [outstring,textstring(k)]; %#ok<AGROW>
        k = k + 1;
        continue
    end
    % Find the end of the name
    if isfield
        j = k + 1;
    elseif isfunction
        j = k + 2;
    else
        j = k + 2;
    end
    while j <= n
        if edisletter(textstring(j)) || edisdigit(textstring(j))
            j = j + 1;
        elseif j+2 <= n && strcmp(textstring(j:j+1),'\_') && ...
                (edisletter(textstring(j+2)) || edisdigit(textstring(j+2)))
            j = j + 2;
        else
            break
        end
    end
    if isfield
        outstring = [outstring,'\texttt{',textstring(k:j-1),'}']; %#ok<AGROW>
    else
        outstring = [outstring,'\EDfile{',textstring(k:j-1),'}']; %#ok<AGROW>
    end
    k = j;
end

end

function yesno = edisletter(c)
yesno = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
end

function yesno = edisdigit(c)
yesno = c >= '0' && c <= '9';
end
