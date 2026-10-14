function report = EDcheckheaders(toolboxfolder,outputfile)
% EDcheckheaders checks that the header of every function in the EDtoolbox
% follows the standard structure, and lists the deviations.
%
% The standard structure of a function file is:
%   function [out1,out2] = EDname(in1,in2)
%   % EDname - A short description (the H1 line, used by lookfor).
%   % More description, one or more paragraphs.
%   %
%   % Input parameters:
%   %   ...
%   %
%   % Output parameters:
%   %   ...
%   %
%   % Peter Svensson (peter.svensson@ntnu.no) <date of the last change>
%   %
%   % [out1,out2] = EDname(in1,in2);
%   <an empty line, without %, so that the history is not shown by help>
%   % <date> <change>
%   % <date> <change>
%
% These are checked for each file:
%   - the function name is the same as the file name
%   - the first help line starts with the function name
%   - there are sections "Input parameters" and "Output parameters" (if
%     the function has input/output parameters)
%   - there is an author line with a date
%   - after the author line, there is a call syntax with the same number of
%     input and output parameters as the function line
%   - after the help text, there is an empty line and then a history block
%   - the date on the author line is not older than the last date in the
%     history block
%
% Input parameters:
%   toolboxfolder   (optional) The top folder of the EDtoolbox. All .m
%                   files in this folder and its subfolders are checked.
%                   Default: the folder above the folder of this function.
%   outputfile      (optional) A text file which the report is written to.
%                   If not given, the report is shown on the screen.
%
% Output parameters:
%   report          A struct array with one element per file that has
%                   remarks, with the fields .file and .remarks (a cell
%                   array of text strings).
%
% Dates are recognized in the form "14 Oct. 2026" (day, month name, year).
% Other date formats in the history are ignored by the date check.
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
%
% report = EDcheckheaders(toolboxfolder,outputfile);

% 14 Oct. 2026 First version (by Claude, prompted by PS)

if nargin < 1 || isempty(toolboxfolder)
    toolboxfolder = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2
    outputfile = '';
end

filelist = dir(fullfile(toolboxfolder,'**','*.m'));
report = struct('file',{},'remarks',{});
nchecked = 0;

for ifile = 1:numel(filelist)
    filename = fullfile(filelist(ifile).folder,filelist(ifile).name);
    [~,name] = fileparts(filename);
    remarks = edcheckonefile(filename,name);
    nchecked = nchecked + 1;
    if ~isempty(remarks)
        report(end+1).file = strrep(filename,[toolboxfolder,filesep],''); %#ok<AGROW>
        report(end).remarks = remarks;
    end
end

textlines = {sprintf('EDcheckheaders: %d files checked, %d without remarks.',nchecked,nchecked-numel(report))};
for ii = 1:numel(report)
    textlines{end+1} = ''; %#ok<AGROW>
    textlines{end+1} = report(ii).file; %#ok<AGROW>
    for jj = 1:numel(report(ii).remarks)
        textlines{end+1} = ['   - ',report(ii).remarks{jj}]; %#ok<AGROW>
    end
end

if isempty(outputfile)
    fprintf('%s\n',textlines{:});
else
    fid = fopen(outputfile,'w');
    if fid == -1
        error('EDcheckheaders:file','The file %s could not be opened for writing',outputfile)
    end
    fprintf(fid,'%s\n',textlines{:});
    fclose(fid);
    disp(['   The report was written to ',outputfile])
end

end

%--------------------------------------------------------------------------
function remarks = edcheckonefile(filename,name)

remarks = {};
lines = splitlines(fileread(filename));

% The function line (possibly continued with ...)
ii = 1;
while ii <= numel(lines) && isempty(strtrim(lines{ii}))
    ii = ii + 1;
end
if ii > numel(lines) || ~startsWith(strtrim(lines{ii}),'function')
    remarks{end+1} = 'a script, or no function line (not checked)';
    return
end
signature = strtrim(lines{ii});
while endsWith(signature,'...') && ii < numel(lines)
    ii = ii + 1;
    signature = [signature(1:end-3),' ',strtrim(lines{ii})];
end
ii = ii + 1;
[fname,ins,outs] = edparsecall(signature,'');
if ~strcmp(fname,name)
    remarks{end+1} = sprintf('the function name "%s" is not the same as the file name',fname);
end

% The help text = the comment lines directly after the function line
helplines = {};
while ii <= numel(lines) && startsWith(strtrim(lines{ii}),'%')
    oneline = strtrim(lines{ii});
    helplines{end+1} = oneline(2:end); %#ok<AGROW>
    ii = ii + 1;
end
if isempty(helplines)
    remarks{end+1} = 'no help text';
    return
end
h1 = strtrim(helplines{1});
if ~startsWith(h1,name)
    remarks{end+1} = sprintf('the first help line does not start with the function name: "%s"',h1(1:min(60,end)));
end
helptext = strjoin(helplines,newline);
if ~isempty(ins) && ~isequal(ins,{'varargin'}) && isempty(regexpi(helptext,'input parameter','once'))
    remarks{end+1} = 'no "Input parameters" section';
end
if ~isempty(outs) && isempty(regexpi(helptext,'output parameter','once'))
    remarks{end+1} = 'no "Output parameters" section';
end

% The author line, and the call syntax after it
authorindex = find(contains(helplines,'Svensson'),1,'last');
headerdate = [];
if isempty(authorindex)
    remarks{end+1} = 'no author line with a date';
else
    headerdates = edfinddates(helplines{authorindex});
    if isempty(headerdates)
        remarks{end+1} = 'the author line has no readable date';
    else
        headerdate = headerdates(1);
    end
    aftertext = strjoin(helplines(authorindex+1:end),newline);
    aftertext = regexprep(aftertext,'\.\.\.\s*\n',' ');
    callline = regexp(aftertext,['(\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?',name,'\s*\([^)]*\)'],'match','once');
    if isempty(callline)
        if ~isempty(ins)
            remarks{end+1} = 'no call syntax after the author line';
        end
    else
        [~,callins,callouts] = edparsecall(callline,name);
        if numel(callins) ~= numel(ins)
            remarks{end+1} = sprintf('the call syntax has %d input parameters, the function has %d',numel(callins),numel(ins));
        end
        if ~isempty(outs) && numel(callouts) ~= numel(outs)
            remarks{end+1} = sprintf('the call syntax has %d output parameters, the function has %d',numel(callouts),numel(outs));
        end
    end
end

% The history block: an empty line, then comment lines
historylines = {};
if ii <= numel(lines) && isempty(strtrim(lines{ii}))
    jj = ii;
    while jj <= numel(lines) && isempty(strtrim(lines{jj}))
        jj = jj + 1;
    end
    while jj <= numel(lines) && startsWith(strtrim(lines{jj}),'%')
        historylines{end+1} = lines{jj}; %#ok<AGROW>
        jj = jj + 1;
    end
end
if isempty(historylines)
    remarks{end+1} = 'no history block (an empty line, then comment lines) after the help text';
elseif ~isempty(headerdate)
    historydates = edfinddates(strjoin(historylines,newline));
    if ~isempty(historydates) && max(historydates) > headerdate
        remarks{end+1} = sprintf('the date on the author line (%s) is older than the last history entry (%s)',...
            datestr(headerdate,'dd mmm yyyy'),datestr(max(historydates),'dd mmm yyyy')); %#ok<DATST>
    end
end

end

%--------------------------------------------------------------------------
function [fname,ins,outs] = edparsecall(textline,knownname)
% Splits a function line, or a call syntax, into the name, the input
% parameter names and the output parameter names.

textline = regexprep(textline,'^\s*function\s*','');
textline = regexprep(textline,';\s*$','');
outs = {};
eqpos = strfind(textline,'=');
if ~isempty(eqpos)
    outpart = textline(1:eqpos(1)-1);
    textline = textline(eqpos(1)+1:end);
    outs = regexp(outpart,'\w+','match');
end
tokens = regexp(textline,'^\s*(\w+)\s*(\(([^)]*)\))?','tokens','once');
if isempty(tokens)
    fname = knownname;
    ins = {};
    return
end
fname = tokens{1};
ins = {};
if numel(tokens) >= 3 && ~isempty(tokens{3})
    ins = strtrim(strsplit(tokens{3},','));
    ins = ins(~cellfun(@isempty,ins));
end

end

%--------------------------------------------------------------------------
function datenumbers = edfinddates(textstring)
% Finds dates of the form "14 Oct. 2026" in a text, as datenums.

monthnames = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
tokens = regexp(textstring,'\<(\d{1,2})\s*[ .]\s*([A-Za-z]{3,9})\.?,?\s+(\d{4})\>','tokens');
datenumbers = [];
for ii = 1:numel(tokens)
    monthnumber = find(strncmpi(tokens{ii}{2},monthnames,3),1);
    daynumber = str2double(tokens{ii}{1});
    if ~isempty(monthnumber) && daynumber >= 1 && daynumber <= 31
        datenumbers(end+1) = datenum(str2double(tokens{ii}{3}),monthnumber,daynumber); %#ok<AGROW,DATNM>
    end
end

end
