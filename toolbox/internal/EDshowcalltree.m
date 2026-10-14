function [treelines,directdeps] = EDshowcalltree(rootfunction,outputfile)
% EDshowcalltree shows the call tree of a function: which functions it
% calls, which functions those call, and so on. Only functions outside the
% Matlab installation are shown (EDtoolbox functions, and others such as
% DataHash and lgwt). The tree is generated from the code itself, with
% Matlab's dependency analysis, so it is always up to date.
%
% A function that has already been shown with its subtree is shown again
% further down with the text '(see above)', and its subtree is not
% repeated.
%
% Input parameters:
%   rootfunction    Name of the function at the top of the tree, e.g.
%                   'EDmain_convex'.
%   outputfile      (optional) Name of a text file which the tree is written
%                   to. If not given, the tree is shown on the screen.
%
% Output parameters:
%   treelines       Cell array with one text line per line of the tree.
%   directdeps      A containers.Map: for each function name, a cell array
%                   with the names of the functions it calls directly.
%                   Can be used to find all the functions that a given
%                   function depends on.
%
% Limitation: the dependency analysis only finds functions that are called
% directly in the code. Calls via function handles stored in variables, or
% via eval, are not found.
%
% Uses the Matlab function matlab.codetools.requiredFilesAndProducts
%
% Peter Svensson (peter.svensson@ntnu.no) 14 Oct. 2026
%
% [treelines,directdeps] = EDshowcalltree(rootfunction,outputfile);

% 14 Oct. 2026 First version (by Claude, prompted by PS)

if nargin < 2
    outputfile = '';
end

rootfile = which(rootfunction);
if isempty(rootfile)
    error('EDshowcalltree:notfound','The function %s can not be found on the Matlab path',rootfunction)
end

directdeps = containers.Map;    % function name -> names of called functions
filepaths  = containers.Map;    % function name -> full file name
expanded   = containers.Map;    % names whose subtree has been shown
treelines  = {};

addnode(rootfile,0,{});

if isempty(outputfile)
    fprintf('%s\n',treelines{:});
else
    fid = fopen(outputfile,'w');
    if fid == -1
        error('EDshowcalltree:file','The file %s could not be opened for writing',outputfile)
    end
    fprintf(fid,'Call tree for %s, generated %s\n\n',rootfunction,char(datetime('now')));
    fprintf(fid,'%s\n',treelines{:});
    fclose(fid);
    disp(['   The call tree was written to ',outputfile])
end

%--------------------------------------------------------------------------
    function addnode(filename,level,ancestors)
        [~,name] = fileparts(filename);
        indent = repmat('   ',1,level);
        if any(strcmp(ancestors,name))
            treelines{end+1,1} = [indent,name,'  (recursive call)'];
            return
        end
        calledfunctions = getdirectdeps(filename,name);
        if isKey(expanded,name)
            if isempty(calledfunctions)
                treelines{end+1,1} = [indent,name];
            else
                treelines{end+1,1} = [indent,name,'  (see above)'];
            end
            return
        end
        treelines{end+1,1} = [indent,name];
        expanded(name) = true;
        for ii = 1:numel(calledfunctions)
            addnode(filepaths(calledfunctions{ii}),level+1,[ancestors,{name}]);
        end
    end

%--------------------------------------------------------------------------
    function calledfunctions = getdirectdeps(filename,name)
        if isKey(directdeps,name)
            calledfunctions = directdeps(name);
            return
        end
        calledfunctions = {};
        [~,~,fileext] = fileparts(filename);
        if strcmpi(fileext,'.m')
            try
                filelist = matlab.codetools.requiredFilesAndProducts(filename,'toponly');
            catch
                filelist = {};
            end
            for jj = 1:numel(filelist)
                if strcmp(filelist{jj},filename) || strncmp(filelist{jj},matlabroot,length(matlabroot))
                    continue
                end
                [~,depname] = fileparts(filelist{jj});
                filepaths(depname) = filelist{jj};
                calledfunctions{end+1} = depname; %#ok<AGROW>
            end
            calledfunctions = sort(calledfunctions);
        end
        directdeps(name) = calledfunctions;
    end

end