function [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,inputorigin] = ...
    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters)
% EDcheckinputstructs - Checks the input structs of EDmain_convex and sets the default values.
% Each of the six input structs is checked against a table (in the local
% function edparamspec), which for each field gives a validator, the
% default value, when the field is relevant, whether it is required, and
% a description. The descriptions can be read with EDinputdoc.
%
% A misspelt field name gives an error. A field which has been removed from
% the toolbox gives a warning and is ignored. A field which is not relevant
% for the settings gets the value NaN. If a geometry file is given, it is
% read here, and .corners and .planecorners are filled in.
%
% Each value gets an origin: 'user' (given by the user), 'default',
% 'derived' (computed from other input) or 'notapplicable'. A struct from
% an earlier run (e.g., EDres.controlparameters) can be sent in again: its
% field .inputorigin is taken out and used, so that a value that was a
% default, and has not been changed, stays 'default'.
%
% Input parameters:
%   geoinputdata, Sinputdata, Rinputdata, envdata, controlparameters,
%   filehandlingparameters
%                   The six input structs of EDmain_convex. For their
%                   fields, see EDinputdoc.
%   A special call, EDcheckinputstructs('paramspecs'), returns a struct
%   with the tables, without checking anything (used by EDinputdoc).
%
% Output parameters:
%   geoinputdata, Sinputdata, Rinputdata, envdata, controlparameters,
%   filehandlingparameters
%                   The six structs, where all fields have values. They do
%                   not contain the field .inputorigin.
%   inputorigin     A struct with one field per input struct, e.g.
%                   inputorigin.controlparameters.ngauss = 'default'.
%
% See also EDinputdoc, EDmain_convex, EDdatastructs
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters,inputorigin] = ...
%    EDcheckinputstructs(geoinputdata,Sinputdata,Rinputdata,envdata,controlparameters,filehandlingparameters);

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
% 11 Oct. 2026 An input struct can carry the field .inputorigin from an
% earlier run. Then 'default' and 'derived' are kept for values that have
% not been changed.
% 12 Oct. 2026 Introduced filehandlingparameters.savefullEDresfile
% 12 Oct. 2026 The defaults for controlparameters.ESIEsourcetermimplementation
% and .ESIEpropagationimplementation were changed from 'legacy' to
% 'accurate'.
% 12 Oct. 2026 controlparameters.difforder is obligatory: the unused
% default value 10 in the table was replaced by NaN, and the header text
% corrected.
% 14 Oct. 2026 filehandlingparameters.saveinteqsousigs and
% .loadinteqsousigs were removed (replaced by the hashed _qfinal file).
% New field .saveinteqsolution (default 1).
% 15 Oct. 2026 The tables were moved to the local function edparamspec, and
% got a sixth column with a description of each field (used by the new
% function EDinputdoc). New call EDcheckinputstructs('paramspecs'). The
% header was updated; the list of fields is now given by EDinputdoc. (By
% Claude, prompted by PS.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Documentation mode: EDcheckinputstructs('paramspecs') returns the tables
% (used by EDinputdoc), without checking anything.

if nargin == 1 && ischar(geoinputdata) && strcmp(geoinputdata,'paramspecs')
    docmode = struct('usesgeofile',true,'usesgeomatrices',true,...
        'geoinputdata',struct('planecorners',[]),'controlparameters',struct('generalmethod','tf'));
    specs = struct;
    structnames = {'geoinputdata','Sinputdata','Rinputdata','envdata','controlparameters','filehandlingparameters'};
    for ii = 1:numel(structnames)
        specs.(structnames{ii}) = edparamspec(structnames{ii},docmode);
    end
    specs.notes = ednotes;
    geoinputdata = specs;
    return
end

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
% If a struct carries an origin record from an earlier run (the field
% .inputorigin), take it out of the struct. It is used further down.

prevorigin = struct;
[geoinputdata,prevorigin.geoinputdata]                     = edtakeoutorigin(geoinputdata);
[Sinputdata,prevorigin.Sinputdata]                         = edtakeoutorigin(Sinputdata);
[Rinputdata,prevorigin.Rinputdata]                         = edtakeoutorigin(Rinputdata);
[envdata,prevorigin.envdata]                               = edtakeoutorigin(envdata);
[controlparameters,prevorigin.controlparameters]           = edtakeoutorigin(controlparameters);
[filehandlingparameters,prevorigin.filehandlingparameters] = edtakeoutorigin(filehandlingparameters);

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

% The tables for the six structs are defined in the local function
% edparamspec (at the end of this file), with the columns
%    name, validator, default, applicability, required, description
% and they are processed by the local function edapplyparamspec. The
% default can be a function handle, which then takes the struct as input.
% The applicability is a function handle which takes the struct as input.
% A field which is not applicable gets the value NaN. The description is
% used by EDinputdoc.

paramspec = edparamspec('controlparameters',struct);



inputorigin = struct;

removedfields = {'inteqnconsecutive','inteqrestimate','intequsemaxcheck','nedgepoints_visibility','specorder'};
[controlparameters,inputorigin.controlparameters] = ...
    edapplyparamspec(controlparameters,paramspec,'controlparameters',{},removedfields,prevorigin.controlparameters);

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
elseif controlparameters.difforder == 2 && strcmp(edprevorigin(prevorigin.controlparameters,'inteqsolutionmethod'),'derived')
    % A struct from an earlier run, where the method was enforced
    inputorigin.controlparameters.inteqsolutionmethod = 'derived';
end

nfrequencies = numel(controlparameters.frequencies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the struct envdata

if ~isstruct(envdata)
    error('ERROR: the struct envdata was not specified')
end

paramspec = edparamspec('envdata',struct);


[envdata,inputorigin.envdata] = edapplyparamspec(envdata,paramspec,'envdata',{},{},prevorigin.envdata);

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

paramspec = edparamspec('geoinputdata',struct('usesgeofile',usesgeofile,...
    'usesgeomatrices',usesgeomatrices,'geoinputdata',geoinputdata));


[geoinputdata,inputorigin.geoinputdata] = edapplyparamspec(geoinputdata,paramspec,'geoinputdata',{},{},prevorigin.geoinputdata);

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

paramspec = edparamspec('Rinputdata',struct);


[Rinputdata,inputorigin.Rinputdata] = edapplyparamspec(Rinputdata,paramspec,'Rinputdata',{},{},prevorigin.Rinputdata);

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
paramspec = edparamspec('Sinputdata',struct);


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
    [pistoninputfields,pistonderivedfields,{'sourceamplitudes'}],{},prevorigin.Sinputdata);

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
    if pistongaussordergiven && ~(strcmp(edprevorigin(prevorigin.Sinputdata,'pistongaussorder'),'default') && isequal(Sinputdata.pistongaussorder,3))
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
    if ~(doaddsourcesgiven && Sinputdata.doaddsources == 1) || strcmp(edprevorigin(prevorigin.Sinputdata,'doaddsources'),'derived')
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
if strcmp(edprevorigin(prevorigin.Sinputdata,'sourceamplitudes'),'default') && all(Sinputdata.sourceamplitudes(:) == 1)
    inputorigin.Sinputdata.sourceamplitudes = 'default';
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
    elseif strcmp(edprevorigin(prevorigin.filehandlingparameters,'outputdirectory'),'derived')
        outputdirectoryderived = isequal(filehandlingparameters.outputdirectory,[geofilepath,filesep,'results']);
    end
    if ~edisgiven(filehandlingparameters,'filestem')
        filehandlingparameters.filestem = geofilestem;
        filestemderived = true;
    elseif strcmp(edprevorigin(prevorigin.filehandlingparameters,'filestem'),'derived')
        filestemderived = isequal(filehandlingparameters.filestem,geofilestem);
    end
end

paramspec = edparamspec('filehandlingparameters',struct('controlparameters',controlparameters));


removedfields = {'savecadgeofile','savesetupfile','lineending','saveEDresfile','saveinteqsousigs','loadinteqsousigs'};
[filehandlingparameters,inputorigin.filehandlingparameters] = ...
    edapplyparamspec(filehandlingparameters,paramspec,'filehandlingparameters',{},removedfields,prevorigin.filehandlingparameters);
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

function [s,origin] = edapplyparamspec(s,paramspec,structname,extrafields,removedfields,prevorigin)
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
%   prevorigin   The origin record from an earlier run, if the struct
%                carried one (see edtakeoutorigin); otherwise an empty
%                struct. A field which was 'default' in the earlier run,
%                and still has the default value, keeps the origin
%                'default'.
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
            if strcmp(edprevorigin(prevorigin,name),'default') && ~required
                % Given, but it was a default value in an earlier run. If
                % the value still is the default value, it stays 'default'.
                try
                    currentdefault = defaultval;
                    if isa(currentdefault,'function_handle')
                        currentdefault = currentdefault(s);
                    end
                    if isequaln(s.(name),currentdefault)
                        origin.(name) = 'default';
                    end
                catch
                end
            end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,prevorigin] = edtakeoutorigin(s)
% If the struct s has the field .inputorigin (attached by EDmain_convex in
% an earlier run), that field is removed from s and returned. Otherwise an
% empty struct is returned.

prevorigin = struct;
if isstruct(s) && isfield(s,'inputorigin')
    if isstruct(s.inputorigin) && isscalar(s.inputorigin)
        prevorigin = s.inputorigin;
    end
    s = rmfield(s,'inputorigin');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function originvalue = edprevorigin(prevorigin,name)
% Returns the origin of the field "name" in an earlier run, or '' if there
% is no such record.

originvalue = '';
if isfield(prevorigin,name) && ischar(prevorigin.(name))
    originvalue = prevorigin.(name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function paramspec = edparamspec(structname,context)
% Returns the table for one input struct, with the columns
%    name, validator, default, applicability, required, description
% context is a struct with values that some of the handles need:
%    .usesgeofile, .usesgeomatrices, .geoinputdata (for geoinputdata)
%    .controlparameters (for filehandlingparameters)
% The rows are processed in order, so a row can depend on fields that were
% handled by earlier rows.

switch structname
    case 'controlparameters'
    paramspec = {
        'generalmethod', @(x) ismember(x,{'ir','tf','tf_ESIEBEM'}), 'tf', @(cp) true, true, ...
            'The calculation method: ''tf'' = transfer functions (ESIE), ''ir'' = impulse responses, ''tf_ESIEBEM'' = transfer functions via surface receivers (ESIEBEM).';
        'frequencies', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && all(x(:) >= 0), NaN, @(cp) ismember(cp.generalmethod,{'tf','tf_ESIEBEM'}), true, ...
            'The frequencies, in Hz. Required for ''tf'' and ''tf_ESIEBEM''.';
        'fs', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 44100, @(cp) strcmp(cp.generalmethod,'ir'), false, ...
            'The sampling frequency, in Hz. Only for ''ir''.';
        'Rstart', @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 0, 0, @(cp) true, false, ...
            'A distance, in m, that gives the time zero (IRs) or the phase reference (TFs): a path of this length has no delay.';
        'directsound', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(cp) true, false, ...
            '1: the direct sound is included. 0: it is not.';
        'skipfirstorder', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(cp) true, false, ...
            '1: the direct sound, the specular reflections and the first-order diffraction are not computed.';
        'difforder', @(x) isnumeric(x) && isscalar(x) && x>=0 && x==round(x), NaN, @(cp) true, true, ...
            'The highest diffraction order. Required. For the GMRES methods, any value larger than 2 gives the full solution of the integral equation.';
        'diff1method', @(x) ismember(x,{'GL','accurate'}), 'accurate', @(cp) cp.difforder>=1, false, ...
            'The method for first-order diffraction: ''accurate'' (adaptive integration, with an analytic treatment near the apex point) or ''GL'' (Gauss-Legendre quadrature, faster).';
        'diff1GLorder', @(x) isnumeric(x) && isreal(x) && isscalar(x) && x > 0 && x==round(x), 16, @(cp) strcmp(cp.diff1method,'GL'), false, ...
            'The number of Gauss-Legendre points per edge, for diff1method = ''GL''.';
        'inteqsolutionmethod', @(x) ismember(x,{'GMRES_fast','GMRES_large','iterative_fixeddifforder'}), 'GMRES_fast', @(cp) cp.difforder>=2, false, ...
            'How the integral equation is solved: ''GMRES_fast'' (with an explicit sparse matrix), ''GMRES_large'' (matrix-free, for large problems) or ''iterative_fixeddifforder'' (difforder-2 terms of the Neumann series). For difforder = 2, the iterative method is always used.';
        'inteqtargettol', @(x) isnumeric(x) && isscalar(x) && x>=0, 1e-10, @(cp) strcmp(cp.inteqsolutionmethod,'GMRES_fast')||strcmp(cp.inteqsolutionmethod,'GMRES_large'), false, ...
            'The relative tolerance for GMRES.';
        'inteqgmresmaxiter', @(x) isnumeric(x) && isscalar(x) && x>=0 && round(x)==x, 100, @(cp) strcmp(cp.inteqsolutionmethod,'GMRES_fast')||strcmp(cp.inteqsolutionmethod,'GMRES_large'), false, ...
            'The maximum number of GMRES iterations.';
        'ESIEsourcetermimplementation', @(x) ismember(x,{'legacy','accurate'}), 'accurate', @(cp) cp.difforder>=2, false, ...
            'The source term of the integral equation: ''accurate'' uses singularity subtraction when a source is close to the plane of a face; ''legacy'' is the earlier method.';
        'ESIEpropagationimplementation', @(x) ismember(x,{'legacy','accurate'}), 'accurate', @(cp) cp.difforder>=2, false, ...
            'The propagation from the edge sources to the receivers: ''accurate'' uses singularity subtraction near the zone boundaries; ''legacy'' uses plain Gauss-Legendre quadrature.';
        'skipreceiverpropagation', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 0, @(cp) cp.difforder>=2, false, ...
            '1: the integral equation is solved, but the propagation to the receivers is skipped (the result is zero).';
        'ngauss', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 16, @(cp) true, false, ...
            'The number of Gauss-Legendre points on the longest edge; shorter edges get proportionally fewer (at least 2). A negative value gives an odd number on every edge. The log file reports the points per wavelength.';
        'discretizationtype', @(x) isnumeric(x) && (x == 0 | x == 1 | x == 2), 2, @(cp) true, false, ...
            'The edge discretization: 2 = Gauss-Legendre (the only one that is supported).';
        'surfacegaussorder', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 5, @(cp) strcmp(cp.generalmethod,'tf_ESIEBEM'), false, ...
            'The quadrature order on the largest surface, for the surface receivers. Only for ''tf_ESIEBEM''.';
        'HODelemsize', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && all(x(:) >= 0), @(cp) 2*2.^(-[0:cp.difforder-1]), @(cp) strcmp(cp.generalmethod,'ir'), false, ...
            'For each diffraction order N, the edge element length is (cair/fs)/HODelemsize(N). Only for ''ir''. Default: 2*2.^(-[0:difforder-1]).';
        'savealldifforders', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(cp) strcmp(cp.generalmethod,'ir'), false, ...
            '1: the higher-order diffraction IRs are kept separate for each order. Only for ''ir''.';
        'saveindividualfirstdiff', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(cp) strcmp(cp.generalmethod,'ir'), false, ...
            '1: the first-order diffraction IR of each edge is kept separate. Only for ''ir''.';
        };
    case 'envdata'
    paramspec = {
        'cair', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 344, @(s) true, false, ...
            'The speed of sound, in m/s.';
        'rhoair', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && isscalar(x) && x > 0, 1.21, @(s) true, false, ...
            'The density of air, in kg/m3.';
        };
    case 'geoinputdata'
    paramspec = {
        'freefieldcase', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 0, @(g) true, false, ...
            '1: no scattering object; only the direct sound is computed.';
        'geoinputfile', @(x) (ischar(x) || (isstring(x) && isscalar(x))) && exist(char(x),'file')==2, NaN, @(g) context.usesgeofile, true, ...
            'The name of a geometry file (.cad). Either this, or .corners and .planecorners, must be given.';
        'corners', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3, NaN, @(g) context.usesgeomatrices, true, ...
            '[ncorners,3] the corner coordinates. If a geometry file is given, corners is filled in from the file.';
        'planecorners', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,1)>=2, NaN, @(g) context.usesgeomatrices, true, ...
            '[nplanes,nmaxcorners] the corner numbers of each plane; the order of the corners defines the front side of the plane according to the RH rule (thumb in normal vector direction). If a geometry file is given, planecorners is filled in from the file.';
        'planerefltypes', @(x) isnumeric(x) && isreal(x) && numel(x)==size(context.geoinputdata.planecorners,1), @(g) ones(size(g.planecorners,1),1), @(g) context.usesgeomatrices, false, ...
            '[nplanes,1] for each plane: +1 rigid, 0 totally absorbing, -1 pressure release. Only with .corners/.planecorners. Default: all rigid.';
        'firstcornertoskip', @(x) isnumeric(x) && isreal(x) && isscalar(x) && x>0, 1e6, @(g) true, false, ...
            'Edges with at least one corner with this number, or higher, are switched off.';
        'listofcornerstoskip', @(x) isnumeric(x) && isreal(x) && all(x(:) >= 0), [], @(g) true, false, ...
            'Edges with at least one corner in this list are switched off.';
        'listofedgestoskip', @(x) isnumeric(x) && isreal(x) && size(x,2)<=2 && all(x(:)>0), [], @(g) true, false, ...
            'Edge numbers (one column), or corner pairs (two columns), of edges that are switched off.';
        'planeseesplanestrategy', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 0, @(g) true, false, ...
            '1: the plane-to-plane visibility is checked only between the plane midpoints (enough for some geometries). 0: a full check.';
        'allownonconvex', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 0, @(g) true, false, ...
            '1: a non-convex model is accepted (EDmain_convex computes correct results only for convex models).';
        };
    case 'Rinputdata'
    paramspec = {
        'coordinates', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3, NaN, @(s) true, true, ...
            '[nreceivers,3] the receiver coordinates. Required.';
        'nedgesubs', @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 2 && round(x) == x, 2, @(s) true, false, ...
            'The number of points per edge that are used in the visibility tests (2 = the two end points).';
        };
    case 'Sinputdata'
        % NB! 'sourcetype' must come before 'coordinates', since the
        % applicability of 'coordinates' depends on it.
    paramspec = {
        'sourcetype', @(x) ismember(x,{'monopole','polygonpiston'}), 'monopole', @(s) true, false, ...
            '''monopole'' or ''polygonpiston'' (see the piston fields below).';
        'coordinates', @(x) isnumeric(x) && isreal(x) && ~isempty(x) && size(x,2)==3, NaN, @(s) strcmp(s.sourcetype,'monopole'), true, ...
            '[nsources,3] the source coordinates. Required for monopoles; for pistons, the center points are computed.';
        'doaddsources', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 0, @(s) true, false, ...
            '1: the contributions of all the sources are added. 0: each source gives a separate result. For a single source, it is always 1.';
        'doallSRcombinations', @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1), 1, @(s) true, false, ...
            '1: all source-receiver combinations are computed. 0: only source N to receiver N (then there must be as many sources as receivers).';
        'nedgesubs', @(x) isnumeric(x) && isreal(x) && isscalar(x) && x >= 2 && round(x) == x, 2, @(s) true, false, ...
            'The number of points per edge that are used in the visibility tests (2 = the two end points).';
        };
    case 'filehandlingparameters'
    paramspec = {
        'outputdirectory', @(x) ischar(x) && ~isempty(x), NaN, @(s) true, true, ...
            'The directory for all the result files. Required, unless a geometry file is given: then the default is the folder ''results'' next to it.';
        'filestem', @(x) ischar(x) && ~isempty(x), NaN, @(s) true, true, ...
            'The first part of the names of all the result files. Required, unless a geometry file is given: then the default is its name.';
        'suppressresultrecycling', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(s) true, false, ...
            '1: saved stage files are never reused, so everything is computed.';
        'showtext', @(x) isnumeric(x) && x >= 0 && round(x)==x, 1, @(s) true, false, ...
            '0: no text on the screen. 1: one line per stage. 2 or more: more details.';
        'saveSRdatafiles', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the stage files _Sdata and _Rdata are saved.';
        'saveeddatafile', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the stage file _eddata is saved.';
        'saveed2datafile', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the stage file _ed2data is saved.';
        'savesubmatrixdata', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the stage file _submatrixdata is saved.';
        'saveinteqsolution', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the solution of the integral equation is saved in the stage file _qfinal (if it is smaller than 4 GB), so that it can be reused for other receivers.';
        'savepathsfile', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the stage file _paths is saved.';
        'savediff2result', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(s) true, false, ...
            '1: the integral equation result is also computed without the second-order diffraction (extraoutputdata.tfinteqdiff_nodiff2 in the _tfinteq file).';
        'savelogfile', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: a log file is written.';
        'savefullEDresfile', @(x) isnumeric(x) && (x == 0 | x == 1), 1, @(s) true, false, ...
            '1: the file _EDres contains the result matrices. 0: they are left out (EDgetresults loads them from the stage files).';
        'saveISEStree', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(s) strcmp(context.controlparameters.generalmethod,'ir'), false, ...
            'Not used by EDmain_convex.';
        'savehodpaths', @(x) isnumeric(x) && (x == 0 | x == 1), 0, @(s) strcmp(context.controlparameters.generalmethod,'ir'), false, ...
            '1: the stage file _hodpaths is saved. Only for ''ir''.';
        };
    otherwise
        error('EDcheckinputstructs:unknownstruct','Unknown struct name: %s',structname)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notes = ednotes
% Descriptions of the input fields that are not in the tables, because they
% are checked separately (used by EDinputdoc).

notes = struct;
notes.Sinputdata = {
    'sourceamplitudes', 'Default: 1. A scalar, a row ([1,nfrequencies]), a column ([nsources,1]) or the full matrix [nsources,nfrequencies]; it is expanded to the full size. For ''ir'', one value per source.';
    'pistoncornercoordinates', '[npistoncorners,3] the corner coordinates of all the pistons. Required for polygon pistons.';
    'pistoncornernumbers', '[npistons,nmaxcorners] the corner numbers of each piston (0 for unused positions). Required for polygon pistons.';
    'pistonplanes', '[npistons,1] the number of the plane that each piston is on. Required for polygon pistons.';
    'pistongaussorder', 'Default: 3. The quadrature order on each piston.';
    };
notes.Sinputdata_text = {'For polygon pistons, the fields .coordinates, .pistongausscoordinates, .pistongaussweights, .ncornersperpiston and .pistonareas are computed (see help EDdatastructs).'};
notes.geoinputdata_text = {'Exactly one geometry source must be given (unless .freefieldcase = 1): .geoinputfile, or .corners and .planecorners.'};
notes.filehandlingparameters_text = {'filehandlingparameters can be left out if a geometry file is given.'};
