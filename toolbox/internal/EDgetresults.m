function [EDres,resultfile] = EDgetresults(varargin)
% EDgetresults - Loads the results struct, EDres, from a results file
% (<filestem>_EDres.mat) that was saved by EDmain_convex.
%
% The results file can be specified in two ways:
%
% 1. By name:
%       EDres = EDgetresults(outputdirectory,filestem);
%
% 2. By the input data. Then the results file which was computed with
%    exactly these input data is found, whatever filestem it has:
%       EDres = EDgetresults(geoinputdata,Sinputdata,Rinputdata,envdata,...
%               controlparameters,filehandlingparameters);
%    These are the same six structs, in the same order, as for
%    EDmain_convex, so the settings script for a calculation can be reused.
%    Only the fields that you would give to EDmain_convex need to be given;
%    the default values are filled in by EDcheckinputstructs, in the same
%    way as when the calculation was run.
%    The directory filehandlingparameters.outputdirectory is searched.
%    The other fields in filehandlingparameters do not influence which
%    file is found.
%    If several results files match, the one with the filestem given in
%    filehandlingparameters.filestem is chosen, if there is one; otherwise
%    the most recent one.
%
% Input parameters:
%   outputdirectory,filestem    Text strings (alternative 1)
%   geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,
%   filehandlingparameters      Structs (alternative 2)
%
% Output parameters:
%   EDres           The results struct, as returned by EDmain_convex. It
%                   contains the results, the six input structs (each with
%                   the field .inputorigin), the EDtoolbox version number
%                   and the hash of the input data. If no results file was
%                   found, EDres is empty. If the calculation has been
%                   run, but did not finish (so that there is a
%                   _settings file but no _EDres file), a message says so.
%   resultfile      The name, with path, of the file that was loaded (empty
%                   if no file was found).
%
% The file _EDres.mat can have been saved without the result matrices
% (filehandlingparameters.savefullEDresfile = 0), to save disk space. Then
% this function loads the result matrices from the _tf, _tfinteq, _ir and
% _irhod files with the same filestem, and puts the results together with
% EDassembleresults, so the returned EDres is the same in both cases.
% Those files must then still exist, and must not have been overwritten
% by another calculation; this is checked.
%
% Results that were computed with EDtoolbox versions before the
% _EDres.mat file was introduced (Oct. 2026) can not be loaded with this
% function.
%
% See also EDmain_convex, EDdatastructs
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [EDres,resultfile] = EDgetresults(outputdirectory,filestem);
% [EDres,resultfile] = EDgetresults(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 1 July 2021 First work
% 16 Sep. 2021 A few more output fields added
% 29 Sep. 2023 Fixed a bug with the size of the ir matrices
% 3 Oct 2023 Made a version 2 with more input parameters
% 16 Apr. 2024 Added the input data structs to EDres
% 11 Oct. 2026 Rewritten (by Claude, prompted by PS). The function now
% loads the _EDres.mat file that EDmain_convex saves, instead of
% collecting the results from the separate files. A results file can be
% found from the input data, via the hash from EDinputhash. The old
% "version 2" call, with seven input parameters in another order, was
% removed.
% 12 Oct. 2026 Handles an _EDres file without the result matrices
% (savefullEDresfile = 0), by loading them from the stage files.
% 14 Oct. 2026 filehandlingparameters is given to EDinputhash, since its
% field .savediff2result is now part of the run hash.
% 15 Oct. 2026 The header was updated (by Claude,
% prompted by PS).

EDres = [];
resultfile = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative 1: the file is specified by its name

if nargin == 2 && ischar(varargin{1}) && ischar(varargin{2})
    outputdirectory = varargin{1};
    filestem = varargin{2};

    filetolookfor = fullfile(outputdirectory,[filestem,'_EDres.mat']);
    if exist(filetolookfor,'file') ~= 2
        disp(['   EDgetresults: the file ',filetolookfor,' was not found.'])
        return
    end
    resultfile = filetolookfor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative 2: the file is specified by the input data

elseif nargin == 6 && isstruct(varargin{1})
    geoinputdata           = varargin{1};
    Sinputdata             = varargin{2};
    Rinputdata             = varargin{3};
    envdata                = varargin{4};
    controlparameters      = varargin{5};
    filehandlingparameters = varargin{6};

    if ~isstruct(filehandlingparameters) || ~isfield(filehandlingparameters,'outputdirectory')
        error('ERROR: filehandlingparameters.outputdirectory must be given; that is the directory which is searched for results files')
    end
    outputdirectory = filehandlingparameters.outputdirectory;
    if exist(outputdirectory,'dir') ~= 7
        disp(['   EDgetresults: the directory ',outputdirectory,' does not exist.'])
        return
    end

    % The filestem is not needed for finding the file, but
    % EDcheckinputstructs requires one.
    filestemgiven = isfield(filehandlingparameters,'filestem') && ischar(filehandlingparameters.filestem);
    if ~filestemgiven
        filehandlingparameters.filestem = 'unspecified';
    end

    % Fill in all the default values, in the same way as when the
    % calculation was run, and compute the hash of the input data.

    [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
        EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

    hashtolookfor = EDinputhash(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

    % Go through all the results files in the directory, and read only the
    % hash from each file.

    filelist = dir(fullfile(outputdirectory,'*_EDres.mat'));
    nfiles = length(filelist);
    ismatch = false(nfiles,1);
    for ii = 1:nfiles
        onefile = fullfile(outputdirectory,filelist(ii).name);
        variablesinfile = whos('-file',onefile);
        if any(strcmp({variablesinfile.name},'EDruninputhash'))
            hashinfile = load(onefile,'EDruninputhash');
            ismatch(ii) = strcmp(hashinfile.EDruninputhash,hashtolookfor);
        end
    end

    filelist = filelist(ismatch);
    nmatches = length(filelist);

    if nmatches == 0
        % No results file. Check if the calculation has been started, but
        % did not finish: the _settings file is saved at the start of a
        % calculation, and it has the same hash.
        settingsfilelist = dir(fullfile(outputdirectory,'*_settings.mat'));
        runswithoutresfile = {};
        for ii = 1:length(settingsfilelist)
            onefile = fullfile(outputdirectory,settingsfilelist(ii).name);
            variablesinfile = whos('-file',onefile);
            if any(strcmp({variablesinfile.name},'EDruninputhash'))
                hashinfile = load(onefile,'EDruninputhash');
                if strcmp(hashinfile.EDruninputhash,hashtolookfor)
                    runswithoutresfile{end+1} = settingsfilelist(ii).name(1:end-length('_settings.mat')); %#ok<AGROW>
                end
            end
        end
        if isempty(runswithoutresfile)
            disp(['   EDgetresults: found no results with these input data, among the ',int2str(nfiles),...
                ' results files in ',outputdirectory])
        else
            disp(['   EDgetresults: a calculation with these input data has been started, with the filestem ''',...
                runswithoutresfile{1},''', but it has no results file (_EDres.mat), so it did not finish.'])
            disp('   Run EDmain_convex again to get the results. Intermediate results that were saved will be reused.')
        end
        return
    end

    % If there are several matching files: choose the one with the given
    % filestem, if there is one; otherwise the most recent one.

    [~,sortvec] = sort([filelist.datenum],'descend');
    filelist = filelist(sortvec);
    ivchoice = 1;
    if filestemgiven
        ivsamestem = find(strcmp({filelist.name},[filehandlingparameters.filestem,'_EDres.mat']));
        if ~isempty(ivsamestem)
            ivchoice = ivsamestem(1);
        end
    end
    resultfile = fullfile(outputdirectory,filelist(ivchoice).name);

    if nmatches > 1
        disp(['   EDgetresults: ',int2str(nmatches),' results files have these input data. Loading ',filelist(ivchoice).name])
    end

else
    error(['ERROR: EDgetresults must be called with either (outputdirectory,filestem) or ',...
        '(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the results struct

filecontents = load(resultfile);
if ~isfield(filecontents,'EDres')
    error(['ERROR: the file ',resultfile,' does not contain the variable EDres'])
end
EDres = filecontents.EDres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the results file does not contain the result matrices (it was saved
% with filehandlingparameters.savefullEDresfile = 0), they are loaded from
% the _tf, _tfinteq, _ir, _irhod files with the same filestem, and put
% together in the same way as in EDmain_convex.

if isfield(filecontents,'EDresisfull') && filecontents.EDresisfull == 0

    [filedirectory,filename] = fileparts(resultfile);
    stagefilestem = fullfile(filedirectory,filename(1:end-length('_EDres')));
    if isfield(filecontents,'EDstagefilehashes')
        stagefilehashes = filecontents.EDstagefilehashes;
    else
        stagefilehashes = struct;
    end

    controlparameters = EDres.controlparameters;
    resultparts = struct;

    if strcmp(controlparameters.generalmethod,'tf')
        if controlparameters.difforder > 1
            stagedata = edloadstagefile(stagefilestem,'tfinteq',{'tfinteqdiff'},stagefilehashes);
            resultparts.tfinteqdiff = stagedata.tfinteqdiff;
        else
            resultparts.tfinteqdiff = [];
        end
        if controlparameters.skipfirstorder == 0
            stagedata = edloadstagefile(stagefilestem,'tf',{'tfdirect','tfgeom','tfdiff'},stagefilehashes);
            resultparts.tfdirect = stagedata.tfdirect;
            resultparts.tfgeom   = stagedata.tfgeom;
            resultparts.tfdiff   = stagedata.tfdiff;
        else
            % The first-order terms were skipped, so they are zeros. The
            % sizes are the same as in EDmain_convex.
            nfrequencies = length(controlparameters.frequencies);
            nsources = size(EDres.Sinputdata.coordinates,1);
            nreceivers = size(EDres.Rinputdata.coordinates,1);
            if nsources == 1
                zeromatrix = zeros(nfrequencies,nreceivers);
            elseif nreceivers == 1
                zeromatrix = zeros(nfrequencies,nsources);
            else
                zeromatrix = zeros(nfrequencies,nreceivers,nsources);
            end
            resultparts.tfdirect = zeromatrix;
            resultparts.tfgeom   = zeromatrix;
            resultparts.tfdiff   = zeromatrix;
        end

    elseif strcmp(controlparameters.generalmethod,'tf_ESIEBEM')
        stagedata = edloadstagefile(stagefilestem,'tfESIEBEM',{'tftot_ESIEBEM'},stagefilehashes);
        resultparts.tftot_ESIEBEM = stagedata.tftot_ESIEBEM;

    elseif strcmp(controlparameters.generalmethod,'ir')
        if controlparameters.skipfirstorder == 0
            stagedata = edloadstagefile(stagefilestem,'ir',{'irdirect','irgeom','irdiff'},stagefilehashes);
            resultparts.irdirect = stagedata.irdirect;
            resultparts.irgeom   = stagedata.irgeom;
            resultparts.irdiff   = stagedata.irdiff;
        else
            resultparts.irdirect = [];
            resultparts.irgeom   = [];
            resultparts.irdiff   = [];
        end
        if controlparameters.difforder > 1
            stagedata = edloadstagefile(stagefilestem,'irhod',{'irhod'},stagefilehashes);
            resultparts.irhod = stagedata.irhod;
        else
            resultparts.irhod = zeros(size(resultparts.irdirect));
        end
    end
    clear stagedata

    resultfields = EDassembleresults(resultparts,controlparameters,EDres.geoinputdata);
    clear resultparts
    resultnames = fieldnames(resultfields);
    for ii = 1:length(resultnames)
        EDres.(resultnames{ii}) = resultfields.(resultnames{ii});
    end
    clear resultfields
end

% Tell the user if the results were computed with another version of the
% EDtoolbox than the present one.

if isfield(EDres,'EDversionnumber')
    EDversionnumber = EDgetversion;
    if ~isequal(EDres.EDversionnumber,EDversionnumber)
        disp(['   EDgetresults: NB! The results in ',resultfile,' were computed with EDtoolbox v. ',...
            num2str(EDres.EDversionnumber),'. The present version is ',num2str(EDversionnumber),'.'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stagedata = edloadstagefile(stagefilestem,suffix,variablenames,stagefilehashes)
% Loads the listed variables from one stage file, <stagefilestem>_<suffix>.mat.
% If a hash for that file was noted when the calculation was run, it is
% checked that the file still has the same hash (i.e., that it has not been
% overwritten by another calculation).

stagefilename = [stagefilestem,'_',suffix,'.mat'];

if exist(stagefilename,'file') ~= 2
    error('EDgetresults:missingstagefile', ...
        ['The results file was saved without the result matrices (savefullEDresfile = 0), and the file ',...
        '%s, which should contain them, does not exist. Run EDmain_convex again to recreate it.'],stagefilename)
end

if isfield(stagefilehashes,suffix)
    variablesinfile = whos('-file',stagefilename);
    hashisthesame = false;
    if any(strcmp({variablesinfile.name},'EDinputdatahash'))
        hashinfile = load(stagefilename,'EDinputdatahash');
        hashisthesame = isequal(hashinfile.EDinputdatahash,stagefilehashes.(suffix));
    end
    if ~hashisthesame
        error('EDgetresults:changedstagefile', ...
            ['The file %s has been changed since the results file was saved (probably by another calculation ',...
            'with the same filestem), so the results can not be put together. Run EDmain_convex again.'],stagefilename)
    end
end

stagedata = load(stagefilename,variablenames{:});
for ii = 1:length(variablenames)
    if ~isfield(stagedata,variablenames{ii})
        error('EDgetresults:missingvariable','The file %s does not contain the variable %s',stagefilename,variablenames{ii})
    end
end
