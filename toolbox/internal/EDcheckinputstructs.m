function [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,inputorigin] = ...
    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDcheckinputstructs checks the input data structs to EDmain_convex
% and sets default values.
% 
% Input parameters:
%   geoinputdata    Exactly one geometry source must be given (unless
%                   .freefieldcase = 1): either .geoinputfile, or
%                   .corners and .planecorners.
%                   .geoinputfile        Name of a geometry file (must
%                                        exist; there is no file open window)
%                   .corners             (alternative to .geoinputfile. The
%                                        use of this option requires that
%                                        filehandlingparameters.outputdirectory
%                                        and .filestem are specified)
%                   .planecorners        (alternative to .geoinputfile)
%                   .planerefltypes      (optional; gives the reflection
%                                         factors of the planes: 1,0,-1.
%                                         default: ones(nplanes,1) )
%                   .firstcornertoskip   (optional, default: 1e6)
%                   .listofcornerstoskip (optional, default [])
%                   .freefieldcase       (optional, default: 0)
%                   .planeseesplanestrategy  (optional, defaul: 0)
%   Sinputdata         .coordinates         (obligatory)
%                   .doaddsources        (default: 0 = no)
%                   .sourceamplitudes    Only used if doaddsources = 1
%                                        (default:
%                                         ones(nsources,nfrequencies))
%                   .doallSRcombinations  (default: 1 = yes)
%                   .sourcetype           'monopole' (default) or
%                                         'polygonpiston'
%                   .pistoncornercoordinates (obligatory if .sourcetype =
%                   'polygonpiston'). Matrix size[npistoncorners,3].
%                   .pistoncornernumbers   (obligatory if .sourcetype =
%                   'polygonpiston'). Matrix size [npistons,maxncornersperpiston].
%                   .pistonplanes          (obligatory if .sourcetype =
%                   'polygonpiston'). List of size [npistons,1].
%                   .nedgesubs           (optional, default 2)
%   Rinputdata         .coordinates         (obligatory)
%                   .nedgesubs           (optional, default 2)
%   envdata         .cair                (default: 344)
%                   .rhoair              (default: 1.21)
%   controlparameters   .fs              (default: 44100) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .directsound         (default: 1 = yes)
%                   .difforder           (default: 10)
%                   .savealldifforders   (default: 0) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .saveindividualfirstdiff (default: 0) Only relevant if
%                                        controlparameters.docalcir = 1.
%                   .skipfirstorder      (default: 0)
%             Maximum one of .docalctf,.docalcir,docalctf_ESIEBEM can be 1.
%                   .docalctf            (default: 0) Only on 
%                   .docalcir            (default: 0)
%                   .docalctf_ESIEBEM    (default: 0)
%                   .Rstart              (default: 0)
%                   .frequencies         (obligatory if .docalctf = 1 or
%                                        docalctf_ESIEBEM = 1)
%                   .discretizationtype  (default: 2 = G-L)
%                   .ngauss              (default: 16, but ignored if
%                                        .docalcir = 1)
%                   .surfacegaussorder   (default: 5, but ignored if 
%                                         .docalctf_ESIEBEM = 0)
%                   .HODirelemsize       (default:
%                                       2*2.^(-[0:controlparameters.difforder-1]))
%   filehandlingparameters    .outputdirectory  (default: the folder
%                                        "results" in the folder of the
%                                        geoinputfile; obligatory if no
%                                        geoinputfile is given)  
%                   .filestem            (default: name of the geoinputfile,
%                                        without extension; obligatory if
%                                        no geoinputfile is given)
%                   .saveSRdatafiles     (default: 1)
%                   .saveeddatafile      (default: 1)
%                   .saveed2datafile      (default: 1)
%                   .savesubmatrixdata    (default: 1)
%                   .saveinteqsousigs     (default: 0)
%                   .loadinteqsousigs     (default: 0)
%                   .savepathsfile        (default: 1)
%                   .saveISEStree         (default: 0) Not used by
%                                         EDmain_convexESIE
%                   .savelogfile          (default: 1)
%                   .savediff2result      (default: 0)
%                   .savehodpaths         (default: 0) Used only if
%                                         .docalcir = 1
%                   .showtext             (default: 1)
% 
% Output parameters:
%   The six input structs, where all fields have got values, and
%   inputorigin     A struct which tells where each value came from. It has
%                   one field per input struct, e.g.,
%                   inputorigin.controlparameters.ngauss, with one of the
%                   values
%                     'user'          given by the user
%                     'default'       not given; the default value was used
%                     'derived'       not given; computed from other input
%                                     (e.g., .corners read from a geometry
%                                     file)
%                     'notapplicable' has no meaning for the current
%                                     settings; the value is NaN
%
% A field name which is not known gives an error (to catch misspelt
% names). A field which has been removed from the toolbox gives a warning
% and is ignored.
%
% Peter Svensson 11 Oct. 2026 (peter.svensson@ntnu.no)
% 
% [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters] = ...
% EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

% 24 Nov. 2017 First version
% 28 Nov. 2017 Cleaned code a bit
% 29 Nov. 2017 Corrected mistake: saveSRdatafiles was called
%              saveSRinputdatafiles. Adjusted to the new indata option:
%              specified corners and planecorners matrices instead of a CAD
%              file.
% 30 Nov. 2017 Fixed a bug where docalctf, instead of docalcir, was set to
%              0. Changed defaults to create eddata and SRfiles, because of
%              plotting the model. Changed one field from logfilename to
%              savelogfile. Removed the field lineending.
% 13 Dec. 2017 Added the input field sourceamplitudes
% 12 Jan. 2018 Forced doaddsources to be 1, if the number of sources = 1.
%              Also followed one yellow recommendation: numel instead of
%              prod.
% 17 Jan 2018  Added the input field planecornertype
% 18 Jan 2018  Fixed a bug with the sourceamplitudes; they were forced to
% one by mistake. Changed default for savelogfile to 1.
% 22 Jan 2018 Removed the field planecornertype in the struct geoinputdata
% 22 Jan 2018 Moved some controlparameter fields away from the convex TF
% case: fs, nedgepoints_visibility, docalcir, specorder.
% 22 Jan 2018 Changed the defaults for saving files.
% 26 Jan 2018 Changed the defaults for saving files. Also changed the
% default output directory to include "results" in the path.
% 26 Jan 2018 Introduced Sinputdata.doallSRcombinations. Default value 1.
% 28 Jan 2018 First version for EDmain_convexESIE_ir
% 31 Jan 2018 Corrected the handling of sourceamplitudes (the freq.
% dependence was not implemented correctly.
% 6 Feb 2018 Introduced a check if number of receivers/sources was zero
% 8 Feb 2018 v 0.109 Introduced a new parameter:
% controlparameters.skipfirstorder (default = 0).
% 14 Feb 2018 v0.112 Small change, assigning a default value to filestem if
% it was not given.
% 14 Feb 2018 Added a check if the functions 'DataHash.m' and 'lgwt.m' are
% available.
% 14 Feb 2018 Added a test if the source and receiver coordinates has the
% right number of columns (3).
% 15 Feb 2018 Removed the savesetupfile. Introduced the saveed2datafile,
% default = 1. Changed defaults for savepathsfile and savesubmatrixdata to
% 1. Stopped assigning default value to saveISEStree. Introduced the
% parameter suppressresultrecycling, default = 0.
% 1 Mar 2018 Introduced the version 3: EDmain_convexESIEBEM
% 1 Mar 2018 Stripped away irrelevant fields; otherwise the DataHash
% doesn't always recognize what is identical.
% 15 Mar 2018 Added the field filehandlingparameters.savehodpaths
% 16 Mar 2018 Combined two error messages
% 21 Mar 2018 Made a few changes to the controlparameters: gave default
% values for docalctf and docalcir. Removed specorder. Introduced the
% .savealldifforders parameter. 
% 21 Mar 2018 Changed savealldifforders to controlparameters instead of
% filehandlingparameters
% 7 Apr 2018 Introduced controlparameters.saveindividualfirstdiff
% 21 Apr 2018 Removed the extra "results" directory in the output
% directory.
% 28 May 2018 Added the input field geoinputdata.planerefltypes
% 22 May 2019 Fixed an error with sourceamplitudes; they did not have the
% right orientation.
% 3 June 2020 Fixed a but with sourceamplitudes that was found by EDdebug:
% If the user had specified a constant sourceamplitudes, it wasn't expanded
% to a matrix of size [nsources, nfrequencies].
% 14 March 2021 The section "% Check the struct Rinputdata" was moved to
% before "% Check the struct Sinputdata"
% 28 Sep. 2023 Adapted to EDmain_convex which does both tf and ir
% 2 Oct. 2023 Changed so that only one of docalftf, docalcir, docalctf_ESIEBEM
% can be set to zero.
% 4 Oct. 2023 Changed default difforder to 10
% 12 Oct. 2023 Introduced the field geoinputdata.freefieldcase
% 26 Oct. 2023 Introduced new fields in Sinputdata which define piston
% sources.
% 29 Oct. 2023 Introduced the new fields .listofcornerstoskip and 
% .planeseesplanestrategy in geoinputdata. Introduced nedgesubs for Sdata
% and Rdata. Added the new field .HODirelemsize
% 30 Oct. 2023 Introduced the piston source checking
% 4 Nov. 2023 Calculates the piston areas.
% 22 Nov. 2023 Fixed a bug: regular polygons had not got the gaussweights
% as a vector, only as a single value.
% 29 Nov. 2023 Changed the name of the filed pistongausspoints to
% pistongaussorder.
% 22 March 2024 Updated the function header (there were references to
% EDmain_convexESIE etc, which was removed).
% 22 March 2024 Removed the requirement that one of the calculation options
% must be 1. If no calculation option is set to 1, then only geometry files
% will be generated, and that is very useful when the geometry is being
% developed.
% 22 March 2024 Added the field listofedgestoskip to geoinputdata
% 17 Sep. 2026 Tried to systematize, using suggestion by Claude. Also
% introduced field .generalmethod = 'tf' or 'ir' or 'tf_ESIEBEM'
% and geoinputdata.allownonconvex
% 20 Sep. 2026 Introduced the controlparameters.ESIEpropagationimplementation
% 24 Sep. 2026 Introduced cp.diff1method and cp.diff1GLorder 
% 8 Oct 2026 Introduced cp.ESIEsourcetermimplementation
% 10 Oct. 2026 Step 2 fixes (by Claude, prompted by PS): validator typos
% (xx, isnumerix, missing comma); .frequencies also applicable (and
% required) for 'tf_ESIEBEM'; .fs default back to 44100; freefieldcase no
% longer blocked by required .corners/.planecorners; piston sources no
% longer blocked by required Sinputdata.coordinates; scalar/row/column
% .sourceamplitudes are expanded again; nreceivers defined; 'saveHODpaths'
% -> 'savehodpaths'; empty S/R coordinates rejected.
% 10 Oct. 2026 For difforder = 2, 'iterative_fixeddifforder' is enforced
% (with a warning if the user had chosen a GMRES method).
% 11 Oct. 2026 Restructuring (by Claude, prompted by PS): the five copies
% of the table-processing loop were replaced by one local function,
% edapplyparamspec, so all six structs are handled in the same way
% (non-applicable fields get the value NaN in all structs). Geometry input:
% exactly one of .geoinputfile and .corners/.planecorners must be given;
% the file open window was removed; .outputdirectory and .filestem get
% defaults from the geometry file name. Removed .savecadgeofile and the
% old code that was commented out.
% 11 Oct. 2026 A field with the value NaN is now treated as not given
% (local function edisgiven), so that a struct which has already been
% through this check can be sent in again.
% 11 Oct. 2026 Step 8 (by Claude, prompted by PS): new output inputorigin
% (user/default/derived/notapplicable for each field). Unknown field names
% give an error; removed fields give a warning. A geometry file is now read
% here, so .corners and .planecorners are always filled in (as 'derived').
% .skipreceiverpropagation is now in the table. The old .docalctf etc
% fields are removed from the struct after they have been translated.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the needed non-EDtoolbox functions are available.

file1ok = (exist('DataHash.m','file') == 2);
file2ok = (exist('lgwt.m','file') == 2);

if file1ok == 0 && file2ok == 0
   error('ERROR: Matlab can not find the functions lgwt.m and DataHash.m. Please download them from Mathworks') 
end
if file1ok == 0
   error('ERROR: Matlab can not find the function DataHash.m. Please download it from Mathworks') 
end
if file2ok == 0
   error('ERROR: Matlab can not find the function lgwt.m. Please download it from Mathworks') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace controlparameters.docalctf = 1 etc with 
% controlparameters.generalmethod = 'tf' or 'ir' or 'tf_ESIEBEM'

if ~isfield(controlparameters,'generalmethod') && any(isfield(controlparameters,{'docalctf','docalcir','docalctf_ESIEBEM'}))
    tf         = isfield(controlparameters,'docalctf')         && controlparameters.docalctf==1;
    ir         = isfield(controlparameters,'docalcir')         && controlparameters.docalcir==1;
    tf_ESIEBEM = isfield(controlparameters,'docalctf_ESIEBEM') && controlparameters.docalctf_ESIEBEM==1;

    if sum([tf,ir,tf_ESIEBEM]) ~= 1
        error('EDinputcheck:calcmethod', ...
            'Exactly one of .docalctf, .docalcir, .docalctf_ESIEBEM must be set to 1')
    elseif tf
        controlparameters.generalmethod = 'tf';
    elseif ir
        controlparameters.generalmethod = 'ir';
    else
        controlparameters.generalmethod = 'tf_ESIEBEM';
    end

    warning('EDinputcheck:deprecated', ...
        ['controlparameters.docalctf/.docalcir/.docalctf_ESIEBEM are deprecated. ', ...
        'This run will keep working via automatic translation, but please update your script to set ', ...
        'controlparameters.generalmethod = ''%s'' directly.'], controlparameters.generalmethod)
end
% The old fields are not needed after this point
oldmethodfields = {'docalctf','docalcir','docalctf_ESIEBEM'};
controlparameters = rmfield(controlparameters,oldmethodfields(isfield(controlparameters,oldmethodfields)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify allowed values and default values for the controlparameters struct 

if ~isstruct(controlparameters)
     error('ERROR: the struct controlparameters was not specified')           
end

% The tables below have these columns:
%    name, validator, default, applicability, required
% and they are processed by the local function edapplyparamspec (at the end
% of this file). The default can be a function handle, which then takes the
% struct as input. The applicability is a function handle which takes the
% struct as input. A field which is not applicable gets the value NaN.

paramspec = {
    'generalmethod',       @(x) ismember(x,{'ir','tf','tf_ESIEBEM'}),                             'tf', @(cp) true, true;
    'frequencies',         @(x) isnumeric(x) && isreal(x) && ~isempty(x) && all(x(:) >= 0),       NaN, @(cp) ismember(cp.generalmethod,{'tf','tf_ESIEBEM'}), true; 
    'fs',                  @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 44100, @(cp) strcmp(cp.generalmethod,'ir'), false; 
    'Rstart',              @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 0,               0,   @(cp) true, false;
    'directsound',         @(x) isnumeric(x) && (x == 0 | x == 1),                                1,   @(cp) true, false;
    'skipfirstorder',      @(x) isnumeric(x) && (x == 0 | x == 1),                                0,   @(cp) true, false;
    'difforder',           @(x) isnumeric(x) && isscalar(x) && x>=0 && x==round(x),               10,  @(cp) true, true;
    'diff1method',         @(x) ismember(x,{'GL','accurate'}),                                   'accurate', @(cp) cp.difforder>=1, false;
    'diff1GLorder',        @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0 && x==round(x), 16, @(cp) strcmp(cp.diff1method,'GL'), false;
    'inteqsolutionmethod', @(x) ismember(x,{'GMRES_fast','GMRES_large','iterative_fixeddifforder'}), 'GMRES_fast', @(cp) cp.difforder>=2, false;
    'inteqtargettol',      @(x) isnumeric(x) && isscalar(x) && x>=0,                              1e-10, @(cp) strcmp(cp.inteqsolutionmethod,'GMRES_fast')||strcmp(cp.inteqsolutionmethod,'GMRES_large'), false;
    'inteqgmresmaxiter',   @(x) isnumeric(x) && isscalar(x) && x>=0 && round(x)==x,               100, @(cp) strcmp(cp.inteqsolutionmethod,'GMRES_fast')||strcmp(cp.inteqsolutionmethod,'GMRES_large'), false;
    'ESIEsourcetermimplementation',  @(x) ismember(x,{'legacy','accurate'}),                     'legacy', @(cp) cp.difforder>=2, false;
    'ESIEpropagationimplementation',  @(x) ismember(x,{'legacy','accurate'}),                     'legacy', @(cp) cp.difforder>=2, false;
    'skipreceiverpropagation',  @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),              0, @(cp) cp.difforder>=2, false;
    'ngauss',              @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 16,  @(cp) true, false;
    'discretizationtype',  @(x) isnumeric(x) && (x == 0 | x == 1 | x == 2),                       2,   @(cp) true, false;
    'surfacegaussorder',   @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 5,   @(cp) strcmp(cp.generalmethod,'tf_ESIEBEM'), false; 
    'HODelemsize',         @(x) isnumeric(x) && isreal(x) && ~isempty(x) && all(x(:) >= 0),       @(cp) 2*2.^(-[0:cp.difforder-1]), @(cp) strcmp(cp.generalmethod,'ir'), false;
    'savealldifforders',   @(x) isnumeric(x) && (x == 0 | x == 1),                                0,   @(cp) strcmp(cp.generalmethod,'ir'), false;
    'saveindividualfirstdiff',   @(x) isnumeric(x) && (x == 0 | x == 1),                                0,   @(cp) strcmp(cp.generalmethod,'ir'), false;
};


inputorigin = struct;

removedfields = {'inteqnconsecutive','inteqrestimate','intequsemaxcheck','nedgepoints_visibility','specorder'};
[controlparameters,inputorigin.controlparameters] = ...
    edapplyparamspec(controlparameters,paramspec,'controlparameters',{},removedfields);

% For difforder = 2, there is no integral equation to solve: the result is
% the source term only. Then GMRES must not be run (the H-matrix is not
% even built), so the iterative method is enforced.
if controlparameters.difforder == 2 && strncmp(controlparameters.inteqsolutionmethod,'GMRES',5)
    if strcmp(inputorigin.controlparameters.inteqsolutionmethod,'user')
        warning('EDinputcheck:methodoverride', ...
            ['controlparameters.inteqsolutionmethod = ''%s'' can not be used for difforder = 2. ', ...
            'It was changed to ''iterative_fixeddifforder''.'], controlparameters.inteqsolutionmethod)
    end
    controlparameters.inteqsolutionmethod = 'iterative_fixeddifforder';
    controlparameters.inteqtargettol = NaN;
    controlparameters.inteqgmresmaxiter = NaN;
    inputorigin.controlparameters.inteqsolutionmethod = 'derived';
    inputorigin.controlparameters.inteqtargettol = 'notapplicable';
    inputorigin.controlparameters.inteqgmresmaxiter = 'notapplicable';
end

nfrequencies = numel(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct envdata

if ~isstruct(envdata)
    error('ERROR: the struct envdata was not specified')
end

paramspec = {
    'cair',         @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 344,  @(s) true, false; 
    'rhoair',       @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 1.21, @(s) true, false; 
};

[envdata,inputorigin.envdata] = edapplyparamspec(envdata,paramspec,'envdata',{},{});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct geoinputdata
%
% Exactly one geometry source must be given, unless freefieldcase = 1:
% either a geometry file (.geoinputfile) or the two matrices .corners and
% .planecorners.
%
% If a geometry file is given, it is read here, and .corners and
% .planecorners are filled in (marked as 'derived' in inputorigin). So,
% after this check, .corners and .planecorners always describe the
% geometry. A struct which has been through this check before has both the
% file name and the matrices. That is accepted if the matrices are
% identical to those in the file.

if ~isstruct(geoinputdata)
    error('ERROR: the struct geoinputdata was not specified')
end

isfreefieldcase = isfield(geoinputdata,'freefieldcase') && isequal(geoinputdata.freefieldcase,1);
usesgeofile     = ~isfreefieldcase && edisgiven(geoinputdata,'geoinputfile');
usesgeomatrices = ~isfreefieldcase && (edisgiven(geoinputdata,'corners') || edisgiven(geoinputdata,'planecorners'));

if ~isfreefieldcase
    if usesgeofile && usesgeomatrices
        % Could be a struct that has been through this check before. Then
        % the matrices must be the same as those in the file.
        [corners_fromfile,planecorners_fromfile] = edreadgeometryfile(geoinputdata.geoinputfile);
        if edisgiven(geoinputdata,'corners') && edisgiven(geoinputdata,'planecorners') && ...
                isequal(geoinputdata.corners,corners_fromfile) && isequal(geoinputdata.planecorners,planecorners_fromfile)
            geoinputdata = rmfield(geoinputdata,{'corners','planecorners'});
            if isfield(geoinputdata,'planerefltypes')
                geoinputdata = rmfield(geoinputdata,'planerefltypes');
            end
            usesgeomatrices = false;
        else
            error('EDinputcheck:geometrysource', ...
                ['Both geoinputdata.geoinputfile and geoinputdata.corners/.planecorners were given, and they do not ',...
                'describe the same geometry. Please give only one geometry source.'])
        end
    end
    if ~usesgeofile && ~usesgeomatrices
        error('EDinputcheck:geometrysource', ...
            'No geometry was given. Please specify either geoinputdata.geoinputfile, or geoinputdata.corners and .planecorners.')
    end
end

paramspec = {
    'freefieldcase',          @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),                  0, @(g) true, false;
    'geoinputfile',           @(x) (ischar(x) || (isstring(x) && isscalar(x))) && exist(char(x),'file')==2,  NaN, @(g) usesgeofile, true; 
    'corners',                @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3,         NaN, @(g) usesgeomatrices, true; 
    'planecorners',           @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,1)>=2,         NaN, @(g) usesgeomatrices, true; 
    'planerefltypes',         @(x) isnumeric(x) && isreal(x) && numel(x)==size(geoinputdata.planecorners,1),   @(g) ones(size(g.planecorners,1),1), @(g) usesgeomatrices, false;
    'firstcornertoskip',      @(x) isnumeric(x) && isreal(x) && isscalar(x) && x>0,                  1e6, @(g) true, false;
    'listofcornerstoskip',    @(x) isnumeric(x) && isreal(x) && all(x(:) >= 0),                      [], @(g) true, false;
    'listofedgestoskip',      @(x) isnumeric(x) && isreal(x) && size(x,2)<=2 && all(x(:)>0),         [], @(g) true, false;
    'planeseesplanestrategy', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),                    0, @(g) true, false;
    'allownonconvex',         @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),                    0, @(g) true, false;
    };

[geoinputdata,inputorigin.geoinputdata] = edapplyparamspec(geoinputdata,paramspec,'geoinputdata',{},{});

if usesgeofile
    geoinputdata.geoinputfile = char(geoinputdata.geoinputfile);
    [geoinputdata.corners,geoinputdata.planecorners] = edreadgeometryfile(geoinputdata.geoinputfile);
    inputorigin.geoinputdata.corners = 'derived';
    inputorigin.geoinputdata.planecorners = 'derived';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Rinputdata

if ~isstruct(Rinputdata)
    error('ERROR: the struct Rinputdata was not specified')
end

paramspec = {
    'coordinates',         @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3,       NaN, @(s) true, true; 
    'nedgesubs',           @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 2 && round(x) == x, 2, @(s) true, false; 
    };

[Rinputdata,inputorigin.Rinputdata] = edapplyparamspec(Rinputdata,paramspec,'Rinputdata',{},{});

nreceivers = size(Rinputdata.coordinates,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct Sinputdata, but the field .sourceamplitudes is handled
% further down, since its size depends on the number of sources, and on
% whether it is a TD case or an FD case.

if ~isstruct(Sinputdata)
    error('ERROR: the struct Sinputdata was not specified')
end

% NB! 'sourcetype' must come before 'coordinates', since the applicability
% of 'coordinates' depends on it (for piston sources, the coordinates are
% computed further down).
paramspec = {
    'sourcetype',          @(x) ismember(x,{'monopole','polygonpiston'}),                       'monopole', @(s) true, false;
    'coordinates',         @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3,       NaN, @(s) strcmp(s.sourcetype,'monopole'), true; 
    'doaddsources',        @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),                  0, @(s) true, false;
    'doallSRcombinations', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1),                  1, @(s) true, false;
    'nedgesubs',           @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 2 && round(x) == x, 2, @(s) true, false; 
    };

% These fields are not in the table, but are handled further down: the
% piston input fields, the fields that are computed for a piston, and
% .sourceamplitudes.
pistoninputfields   = {'pistoncornercoordinates','pistoncornernumbers','pistonplanes','pistongaussorder'};
pistonderivedfields = {'pistongausscoordinates','pistongaussweights','ncornersperpiston','pistonareas'};

% For a piston, the coordinates are computed further down. If they are
% present (as in a struct which has been through this check before), they
% are removed here.
if isfield(Sinputdata,'sourcetype') && isequal(Sinputdata.sourcetype,'polygonpiston') && isfield(Sinputdata,'coordinates')
    Sinputdata = rmfield(Sinputdata,'coordinates');
end

sourceamplitudesgiven = edisgiven(Sinputdata,'sourceamplitudes');
pistongaussordergiven = edisgiven(Sinputdata,'pistongaussorder');
doaddsourcesgiven     = edisgiven(Sinputdata,'doaddsources');

[Sinputdata,inputorigin.Sinputdata] = edapplyparamspec(Sinputdata,paramspec,'Sinputdata',...
    [pistoninputfields,pistonderivedfields,{'sourceamplitudes'}],{});

if strcmp(Sinputdata.sourcetype,'monopole')
    nsources = size(Sinputdata.coordinates,1);
    Sinputdata.pistoncornercoordinates = [];
    Sinputdata.pistoncornernumbers = [];
    Sinputdata.pistonplanes = []; 
    Sinputdata.pistongaussorder = 0;
    for ii = 1:length(pistoninputfields)
        inputorigin.Sinputdata.(pistoninputfields{ii}) = 'notapplicable';
    end
    
else   % Sinputdata.sourcetype is 'polygonpiston'
   if ~isfield(Sinputdata,'pistoncornercoordinates') || ~isfield(Sinputdata,'pistoncornernumbers') || ~isfield(Sinputdata,'pistonplanes')
        error('ERROR: when you specify a polygonpiston sources, you must specify corner coordinates and numbers, and piston plane(s).')
    end
    npistons = size(Sinputdata.pistoncornernumbers,1);
    [npoints,ncolumns] = size(Sinputdata.pistoncornercoordinates);
    if npistons == 0
        error('ERROR: You must define at least one piston source by specifying the corner numbers')
    else
        nsources = npistons;
    end
    if npoints < 3
        error('ERROR: You must define at least three pistoncornercoordinates')
    end
    if ncolumns ~= 3
        error('ERROR: There should be 3 columns for the pistoncornercoordinates')
    end        
    if ~isfield(Sinputdata,'pistongaussorder')
        Sinputdata.pistongaussorder = 3;
    end

    pistonmidpointcoordinates = zeros(npistons,3);   
    ncornersperpiston = sum(Sinputdata.pistoncornernumbers.'>0).';
    pistonareavec = zeros(npistons,1);

    Sinputdata.pistongausscoordinates = cell(npistons,1);
    Sinputdata.pistongaussweights = cell(npistons,1);

    for ii = 1:npistons
        corners_onepiston = Sinputdata.pistoncornercoordinates...
            (Sinputdata.pistoncornernumbers(ii,1:ncornersperpiston(ii)),:);
        pistonmidpointcoordinates(ii,:) = mean(corners_onepiston);

        % Calculate the area of each piston. Implemented for triangles
        % and rectangles, using a version of Heron's formula for the
        % triangle: A = 1/4*sqrt( 4*a^2*b^2 - (a^2 + b^2 - c^2)^2 )
        % where a,b,c = triangle side lengths
        corners_onepiston_expanded = [corners_onepiston;corners_onepiston(1,:)];
        sidevectors = diff(corners_onepiston_expanded);
        sidelengthsquared = sum(sidevectors.^2,2);

        if ncornersperpiston(ii) == 3
            pistonarea = 1/4*sqrt( 4*sidelengthsquared(1)*sidelengthsquared(2) ...
                - (sidelengthsquared(1) + sidelengthsquared(2) - sidelengthsquared(3))^2  );
            switch Sinputdata.pistongaussorder
                case {0,1}
                   gaussvalues = [1/3 1/3 1]; 
                case {2}
                    gaussvalues = [2/3 1/6 1/3;...
                                   1/6 2/3 1/3;
                                   1/6 1/6 1/3];
                case {3}
                    gaussvalues = [1/3 1/3 -9/32; ...
                        0.2 0.6 25/96; ...
                        0.2 0.2 25/96; ...
                        0.6 0.2 25/96];
                case {4}
                    gaussvalues = [0.108103018168070  0.445948490915965 0.223381589678011; ...
                                   0.445948490915965  0.108103018168070 0.223381589678011; ...
                                   0.445948490915965  0.445948490915965 0.223381589678011; ...
                                   0.816847572980459  0.091576213509771 0.109951743655322; ...
                                   0.091576213509771  0.816847572980459 0.109951743655322; ...
                                   0.091576213509771  0.091576213509771 0.109951743655322];
                case {5}   % 7 gauss points
                    gaussvalues = [1/3  1/3  0.225;...
                                   0.059715871789770 0.470142064105115 0.132394152788506;...
                                   0.470142064105115 0.059715871789770 0.132394152788506;...
                                   0.470142064105115 0.470142064105115 0.132394152788506;...
                                   0.797426985353087 0.101286507323456 0.125939180544827;...
                                   0.101286507323456 0.797426985353087 0.125939180544827;...
                                   0.101286507323456 0.101286507323456 0.125939180544827];
                case{6}  % 12 gauss points
                    gaussvalues = [0.501426509658179 0.249286745170910 0.116786275726379;...
                                   0.249286745170910 0.501426509658179 0.116786275726379;...
                                   0.249286745170910 0.249286745170910 0.116786275726379;...
                                   0.873821971016996 0.063089014491502 0.050844906370207;...
                                   0.063089014491502 0.873821971016996 0.050844906370207;...
                                   0.063089014491502 0.063089014491502 0.050844906370207;...                           
                                   0.053145049844817 0.310352451033784 0.082851075618374;...
                                   0.310352451033784 0.053145049844817 0.082851075618374;...
                                   0.053145049844817 0.636502499121399 0.082851075618374;...
                                   0.636502499121399 0.053145049844817 0.082851075618374;...
                                   0.636502499121399 0.310352451033784 0.082851075618374;...
                                   0.310352451033784 0.636502499121399 0.082851075618374];
                case{7,8} % 16 gauss points 
                      gaussvalues = [1/3 1/3 0.144315607677787;...
                                     0.081414823414554 0.459292588292723 0.095091634267285;...
                                     0.459292588292723 0.081414823414554  0.095091634267285;...
                                     0.459292588292723 0.459292588292723 0.095091634267285;...
                                     0.658861384496480 0.170569307751760 0.103217370534718;...
                                     0.170569307751760 0.658861384496480 0.103217370534718;...
                                     0.170569307751760 0.170569307751760 0.103217370534718;...
                                     0.898905543365938 0.050547228317031 0.032458497623198;...
                                     0.050547228317031 0.898905543365938 0.032458497623198;...
                                     0.050547228317031 0.050547228317031 0.032458497623198;...
                                     0.008394777409958 0.263112829634638 0.027230314174435;...
                                     0.728492392955404 0.008394777409958 0.027230314174435;...
                                     0.008394777409958 0.728492392955404 0.027230314174435;...
                                     0.728492392955404 0.263112829634638 0.027230314174435;...
                                     0.263112829634638 0.008394777409958 0.027230314174435;...
                                     0.263112829634638 0.728492392955404 0.027230314174435];
                                                               
                otherwise
                    error(['ERROR: triangle quadrature not implemented for this number: ',int2str(Sinputdata.pistongaussorder)])
            end
            n2 = size(gaussvalues,1);  
%             corners_x = planedata.corners(planedata.planecorners(planenumber,1:3),1);
%             corners_y = planedata.corners(planedata.planecorners(planenumber,1:3),2);
%             corners_z = planedata.corners(planedata.planecorners(planenumber,1:3),3);
            corners_x = corners_onepiston(:,1);
            corners_y = corners_onepiston(:,2);
            corners_z = corners_onepiston(:,3);
            
            % Convert the normalized triangle coordinates to the 3D coordinates
            
            psi=[gaussvalues(:,1) gaussvalues(:,2) 1-gaussvalues(:,1)-gaussvalues(:,2)];
            
            xgausspoints = psi*corners_x;
            ygausspoints = psi*corners_y;
            zgausspoints = psi*corners_z;
            Sinputdata.pistongausscoordinates{ii} = [xgausspoints ygausspoints zgausspoints];

            Sinputdata.pistongaussweights{ii} = gaussvalues(:,3);

        elseif ncornersperpiston(ii) == 4
            extradistsquared = sum( (corners_onepiston_expanded(2,:) - corners_onepiston_expanded(4,:)).^2 );
            triarea1 =  1/4*sqrt( 4*sidelengthsquared(1)*sidelengthsquared(2) ...
                - (sidelengthsquared(1) + sidelengthsquared(2) -extradistsquared)^2  );
            triarea2 =  1/4*sqrt( 4*sidelengthsquared(3)*sidelengthsquared(4) ...
                - (sidelengthsquared(3) + sidelengthsquared(3) -extradistsquared)^2  );
            pistonarea = triarea1 + triarea2;

           n2 = Sinputdata.pistongaussorder^2;
           [x,w] = lgwt(Sinputdata.pistongaussorder,0,1);
           
           x = x(end:-1:1);
           xy = [repmat(x,Sinputdata.pistongaussorder,1) reshape(x(:,ones(1,Sinputdata.pistongaussorder)).',n2,1)];
           ww = prod([repmat(w,Sinputdata.pistongaussorder,1) reshape(w(:,ones(1,Sinputdata.pistongaussorder)).',n2,1)],2);   
            
           c1 = corners_onepiston(1,:);
           c2 = corners_onepiston(2,:);
           c3 = corners_onepiston(3,:);
           c4 = corners_onepiston(4,:);
           xvec1 = c2 - c1;
           alen = norm(xvec1);
           xvec2 = c3 - c4;
           clen = norm(xvec2);
           yvec1 = c4 - c1;
           dlen = norm(yvec1);
           yvec2 = c3 - c2;
           blen = norm(yvec2);
           crossvec1 = c3 - c1;
           plen = norm(crossvec1);
           crossvec2 = c4 - c2;
           qlen = norm(crossvec2);         

           startpointsx = c1(ones(n2,1),:);
           startpointsx = startpointsx + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);
    
           endpointsx = c4(ones(n2,1),:);
           endpointsx = endpointsx + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
           
           startpointsy = c1(ones(n2,1),:);
           startpointsy = startpointsy + yvec1(ones(n2,1),:).*xy(:,[2 2 2]);
    
           endpointsy = c2(ones(n2,1),:);
           endpointsy = endpointsy + yvec2(ones(n2,1),:).*xy(:,[2 2 2]);

           startpointsx = c1(ones(n2,1),:);
           startpointsx = startpointsx + xvec1(ones(n2,1),:).*xy(:,[1 1 1]);
    
           endpointsx = c4(ones(n2,1),:);
           endpointsx = endpointsx + xvec2(ones(n2,1),:).*xy(:,[1 1 1]);
           
           startpointsy = c1(ones(n2,1),:);
           startpointsy = startpointsy + yvec1(ones(n2,1),:).*xy(:,[2 2 2]);
    
           endpointsy = c2(ones(n2,1),:);
           endpointsy = endpointsy + yvec2(ones(n2,1),:).*xy(:,[2 2 2]);
    
           Sinputdata.pistongausscoordinates{ii} = startpointsx + (endpointsx - startpointsx).*xy(:,[2 2 2]);
           Sinputdata.pistongaussweights{ii} = ww;
        else
            % Check if the polygon is regular; then we can compute its area
            if std(sidelengthsquared)/mean(sidelengthsquared) > 1e-9
                error('ERROR: piston areas can only be calculated for 3- and 4-sided pistons, and regular polygons.')
            else
                cornerradius = mean( EDcalcdist(corners_onepiston,pistonmidpointcoordinates(ii,:)) );
                pistonarea = ncornersperpiston(ii)/2*cornerradius^2*sin(2*pi/ncornersperpiston(ii));

                % We should not use cornerradius; it would be a bit more
                % accurate to use the theoretical circle radius
                circlecoordinates = GEOcirclepoints(cornerradius,Sinputdata.pistongaussorder);
                
                Sinputdata.pistongausscoordinates{ii} = circlecoordinates + pistonmidpointcoordinates(ii,:);
                Sinputdata.pistongaussweights{ii} = 1/size(circlecoordinates,1)*ones(size(circlecoordinates,1),1);
                disp('WARNING! Regular polygon piston has only been properly implemented in the plane with nvec = [0 0 1]')
            end
        end
        pistonareavec(ii) = pistonarea;
    end
    Sinputdata.coordinates = pistonmidpointcoordinates;
    Sinputdata.ncornersperpiston = ncornersperpiston;
    Sinputdata.pistonareas = pistonareavec;

    inputorigin.Sinputdata.coordinates = 'derived';
    inputorigin.Sinputdata.pistoncornercoordinates = 'user';
    inputorigin.Sinputdata.pistoncornernumbers = 'user';
    inputorigin.Sinputdata.pistonplanes = 'user';
    if pistongaussordergiven
        inputorigin.Sinputdata.pistongaussorder = 'user';
    else
        inputorigin.Sinputdata.pistongaussorder = 'default';
    end
    for ii = 1:length(pistonderivedfields)
        inputorigin.Sinputdata.(pistonderivedfields{ii}) = 'derived';
    end

    if ~isfield(Sinputdata,'pistonplanes')
        error('ERROR: When Sinputdata.sourcetype is polygonpiston, then the field .pistonplanes must be specified.')
    else
        if size(Sinputdata.pistonplanes(:),1) ~= npistons
            error('ERROR: the fields .pistoncornernumbers and .pistonplanes must have the same number of rows.')
        end
    end

end

if nsources == 1
    if ~(doaddsourcesgiven && Sinputdata.doaddsources == 1)
        inputorigin.Sinputdata.doaddsources = 'derived';
    end
    Sinputdata.doaddsources = 1;
end    
if Sinputdata.doallSRcombinations == 0 && nsources ~= nreceivers
    disp(['   nsources = ',int2str(nsources)])
    disp(['   nreceivers = ',int2str(nreceivers)])
    error('ERROR: doallSRcombinations was set to 0, but the number of sources was not the same as the number of receivers'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the field .sourceamplitudes of the struct Sinputdata, and expand it
% to size [nsources,nfrequencies] (FD cases) or [nsources,1] (TD case).
% Allowed input sizes: scalar, one row, one column, or the full size.

if ismember(controlparameters.generalmethod,{'tf','tf_ESIEBEM'})
    nampcolumns = nfrequencies;
else
    nampcolumns = 1;
end
if sourceamplitudesgiven
    inputorigin.Sinputdata.sourceamplitudes = 'user';
else
    Sinputdata.sourceamplitudes = 1;
    inputorigin.Sinputdata.sourceamplitudes = 'default';
end
if ~isnumeric(Sinputdata.sourceamplitudes) || isempty(Sinputdata.sourceamplitudes) || ~ismatrix(Sinputdata.sourceamplitudes)
    error('EDinputcheck:badparam','Invalid value for Sinputdata.sourceamplitudes')
end
[n1,n2] = size(Sinputdata.sourceamplitudes);
if n1 ~= 1 && n1 ~= nsources
    error('EDinputcheck:badparam',['Sinputdata.sourceamplitudes must have one row, or one row per source, ',...
        'but it had the size ',int2str(n1),' by ',int2str(n2)])
end
if n2 ~= 1 && n2 ~= nampcolumns
    error('EDinputcheck:badparam',['Sinputdata.sourceamplitudes must have one column, or (for FD calculations) one column ',...
        'per frequency, but it had the size ',int2str(n1),' by ',int2str(n2)])
end
if n1 == 1
    Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(ones(nsources,1),:);
end
if n2 == 1
    Sinputdata.sourceamplitudes = Sinputdata.sourceamplitudes(:,ones(1,nampcolumns));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct filehandlingparameters
%
% If a geometry file was given, then .outputdirectory and .filestem get
% default values from the file name. Otherwise they must be given.

if ~isstruct(filehandlingparameters)
    error('ERROR: the struct filehandlingparameters was not specified')
end

outputdirectoryderived = false;
filestemderived = false;
if usesgeofile
    [geofilepath,geofilestem] = fileparts(geoinputdata.geoinputfile);
    if isempty(geofilepath)
        geofilepath = pwd;
    end
    if ~edisgiven(filehandlingparameters,'outputdirectory')
        filehandlingparameters.outputdirectory = [geofilepath,filesep,'results'];
        outputdirectoryderived = true;
    end
    if ~edisgiven(filehandlingparameters,'filestem')
        filehandlingparameters.filestem = geofilestem;
        filestemderived = true;
    end
end

paramspec = {
    'outputdirectory',          @(x) ischar(x) && ~isempty(x),              NaN, @(s) true, true; 
    'filestem',                 @(x) ischar(x) && ~isempty(x),              NaN, @(s) true, true; 
    'suppressresultrecycling',  @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) true, false; 
    'showtext',                 @(x) isnumeric(x) && x >= 0 && round(x)==x, 1, @(s) true, false; 
    'saveSRdatafiles',          @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'saveeddatafile',           @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'saveed2datafile',          @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'savesubmatrixdata',        @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'saveinteqsousigs',         @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) true, false; 
    'loadinteqsousigs',         @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) true, false; 
    'savepathsfile',            @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'savediff2result',          @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) true, false; 
    'savelogfile',              @(x) isnumeric(x) && (x == 0 | x == 1),     1, @(s) true, false; 
    'saveISEStree',             @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) strcmp(controlparameters.generalmethod,'ir'), false; 
    'savehodpaths',             @(x) isnumeric(x) && (x == 0 | x == 1),     0, @(s) strcmp(controlparameters.generalmethod,'ir'), false; 
    };

removedfields = {'savecadgeofile','savesetupfile','lineending'};
[filehandlingparameters,inputorigin.filehandlingparameters] = ...
    edapplyparamspec(filehandlingparameters,paramspec,'filehandlingparameters',{},removedfields);
if outputdirectoryderived
    inputorigin.filehandlingparameters.outputdirectory = 'derived';
end
if filestemderived
    inputorigin.filehandlingparameters.filestem = 'derived';
end

if exist(filehandlingparameters.outputdirectory,'dir') ~= 7
    mkdir(filehandlingparameters.outputdirectory)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,origin] = edapplyparamspec(s,paramspec,structname,extrafields,removedfields)
% Checks the fields of one input struct against a table, and sets default
% values.
%
% Input parameters:
%   s            The struct to check
%   paramspec    Cell array of size [nparameters,5], with the columns
%                   1. name          The field name
%                   2. validator     Function handle; takes the field value
%                                    and returns true if it is valid.
%                   3. default       The default value, or a function handle
%                                    which takes the struct and returns
%                                    the default value.
%                   4. applicability Function handle; takes the struct and
%                                    returns true if the field is relevant
%                                    for the current settings.
%                   5. required      true if the user must give the field
%                                    (when it is applicable).
%                The rows are processed in order, so a row can depend on
%                fields that were handled by earlier rows.
%   structname   The name of the struct, for the error messages
%   extrafields  Cell array with names of fields that are allowed in the
%                struct, but are not in the table (they are handled
%                elsewhere).
%   removedfields  Cell array with names of fields that existed in earlier
%                versions of the toolbox. They give a warning, and are
%                removed from the struct.
%
% Output parameters:
%   s            The struct, where all fields in the table exist. Fields
%                that are not applicable have the value NaN.
%   origin       A struct with one field for each row in the table, with
%                the value 'user', 'default' or 'notapplicable'.
%
% A field in the struct which is neither in the table, nor in extrafields
% or removedfields, gives an error. This catches misspelt field names.
%
% A field which has the value NaN is treated as not given. This makes it
% possible to send in a struct which has already been through this check
% (e.g., EDres.controlparameters from an earlier run).

isremoved = isfield(s,removedfields);
if any(isremoved)
    removednames = removedfields(isremoved);
    for k = 1:length(removednames)
        warning('EDinputcheck:removedparam', ...
            '%s.%s is no longer used by the EDtoolbox. It is ignored.',structname,removednames{k})
    end
    s = rmfield(s,removednames);
end

unknownnames = setdiff(fieldnames(s),[paramspec(:,1);extrafields(:)]);
if ~isempty(unknownnames)
    error('EDinputcheck:unknownparam', ...
        '%s.%s is not a known field. Please check the spelling. The known fields are: %s', ...
        structname,unknownnames{1},strjoin(sort([paramspec(:,1);extrafields(:)]).',', '))
end

origin = struct;

for k = 1:size(paramspec,1)
    name       = paramspec{k,1};
    validator  = paramspec{k,2};
    defaultval = paramspec{k,3};
    applicable = paramspec{k,4}(s);
    required   = paramspec{k,5};

    if applicable
        if ~edisgiven(s,name)
            if required
                error('EDinputcheck:missingparam', ...
                    '%s.%s is required and must be supplied by the user (no default exists)',structname,name)
            else
                if isa(defaultval,'function_handle')
                    defaultval = defaultval(s);
                end
                s.(name) = defaultval;
                origin.(name) = 'default';
            end
        else
            origin.(name) = 'user';
            % A validator that fails, or that does not return a single
            % true/false, counts as an invalid value.
            try
                isvalid = validator(s.(name));
                isvalid = islogical(isvalid) && isscalar(isvalid) && isvalid;
            catch
                isvalid = false;
            end
            if ~isvalid
                error('EDinputcheck:badparam','Invalid value for %s.%s',structname,name)
            end
        end
    else
        if isfield(s,name)
            val = s.(name);
            if ~(isnumeric(val) && isscalar(val) && isnan(val))
                warning('EDinputcheck:irrelevantparam', ...
                    '%s.%s is set but has no effect for the current settings (set to NaN)',structname,name)
            end
        end
        s.(name) = NaN;
        origin.(name) = 'notapplicable';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function given = edisgiven(s,name)
% Returns true if the field exists in the struct s and has a value other
% than a single NaN. NaN is the marker for "not applicable", so a field
% with that value counts as not given by the user.

given = isfield(s,name);
if given
    val = s.(name);
    given = ~(isnumeric(val) && isscalar(val) && isnan(val));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [corners,planecorners] = edreadgeometryfile(geoinputfile)
% Reads a geometry file and returns the two geometry matrices. The file
% format is chosen from the file extension. To add a format: add a case
% here, which calls the reader function for that format.

[~,~,fileext] = fileparts(char(geoinputfile));

switch lower(fileext)
    case '.cad'
        planedata = EDreadcad(char(geoinputfile),0);
    otherwise
        error('EDinputcheck:geometryformat', ...
            'The geometry file format ''%s'' is not supported. Supported formats: .cad',fileext)
end

corners = planedata.corners;
planecorners = double(planedata.planecorners);
