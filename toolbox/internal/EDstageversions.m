function stageinfo = EDstageversions(stagename)
% EDstageversions keeps the version numbers of the EDtoolbox functions that
% compute results, and, for each stage function (a function that saves a
% stage file, such as _eddata or _tf), the list of subroutines it uses.
%
% Each stage function includes the output of EDstageversions(itsname) in
% its EDinputdatahash. So, if the version of the stage function, or of
% any subroutine it uses, is increased, previously saved files of that
% stage are not reused, and neither are the files of all later stages
% (because of the chained hashes).
%
% HOW TO USE
%   * When a change in a function changes its results: increase its
%     number in the "versions" table below by one.
%   * When a stage function (or one of its subroutines) starts to call a
%     new function which affects the results: add that function to the
%     stage's list in the "stagesubroutines" table, and to "versions".
%   * Functions that do not affect any results (text output, file
%     handling, hashing) go in the "ignored" list instead.
%   EDverify checks that the lists agree with the actual calls in the code.
%
% Input parameters:
%   stagename   (optional) Name of a stage function, e.g. 'EDedgeo'.
%
% Output parameters:
%   stageinfo   If stagename is given: a struct with the fields
%                 .stage        the stage name
%                 .version      the version of the stage function
%                 .subroutines  a struct with one field per subroutine,
%                               with its version number
%               If no input parameter is given: a struct with the fields
%                 .versions          all version numbers (struct)
%                 .stagesubroutines  the subroutine list of each stage
%                                    (struct of cell arrays)
%                 .ignored           cell array with the ignored functions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% stageinfo = EDstageversions(stagename);

% 14 Oct. 2026 First version (by Claude, prompted by PS). The lists were
% made from the output of EDshowcalltree('EDmain_convex').
% 14 Oct. 2026 EDmakeHODirs, EDwedge1st_fd and EDpropagateESIEBEM to
% version 2 (bug fixes).
% 15 Oct. 2026 EDwedge1st_fd to version 3 (near-singular cases with 'GL').
% 15 Oct. 2026 EDcalcpropagatematrix to version 2 (bug fix).

%--------------------------------------------------------------------------
% Version numbers

versions = struct;

% Stage functions
versions.EDedgeo                      = 1;
versions.EDSorRgeo                    = 1;
versions.EDfindconvexGApaths          = 1;
versions.EDmakefirstordertfs          = 1;
versions.EDmakefirstorderirs          = 1;
versions.EDed2geo                     = 1;
versions.EDinteg_submatrixstructure   = 1;
versions.EDintegralequation_convex_tf = 1;
versions.EDfindHODpaths               = 1;
versions.EDmakeHODirs                 = 2;  % 2: source amplitudes applied also for doaddsources = 0
versions.EDpropagateESIEBEM           = 2;  % 2: complete hash; several sources give an error

% Subroutines
versions.EDbetaoverm_fd               = 1;
versions.EDbetaoverml                 = 1;
versions.EDbetaoverml_fd              = 1;
versions.EDbuildHsparse               = 1;
versions.EDcalccosfi                  = 1;
versions.EDcalcdist                   = 1;
versions.EDcalcedgeinteqmatrixsub2_mex = 1;  % Compiled file: increase by
                                             % hand when it is recompiled
                                             % from changed source code.
versions.EDcalcpropagatematrix        = 2;  % 2: bug fix, edge data for edges with the same number of nodes
versions.EDcheckobstr_edgetoedge      = 1;
versions.EDcheckobstr_pointtoedge     = 1;
versions.EDcheckobstrpaths            = 1;
versions.EDchkISvisible               = 1;
versions.EDchkISvisiblex              = 1;
versions.EDcompress3                  = 1;
versions.EDcompress7                  = 1;
versions.EDcoordtrans1                = 1;
versions.EDcoordtrans2                = 1;
versions.EDcreindexmatrix             = 1;
versions.EDcross                      = 1;
versions.EDdistelements               = 1;
versions.EDfindis                     = 1;
versions.EDgetedgepoints              = 1;
versions.EDgetGLrule                  = 1;
versions.EDinfrontofplane             = 1;
versions.EDinteg_souterm              = 1;
versions.EDintegratebetaoverm         = 1;
versions.EDmatvecH                    = 1;
versions.EDpistonedgeintegrand        = 1;
versions.EDpoinpla                    = 1;
versions.EDpoinplax                   = 1;
versions.EDquadgk                     = 1;
versions.EDquadstep                   = 1;
versions.EDwedge1st_fd                = 3;  % 2: 'GL' falls back to 'accurate' for singular cases; 3: also near-singular
versions.EDwedge1st_ir                = 1;
versions.EDwedge2nd                   = 1;
versions.EDwedgeN                     = 1;
versions.lgwt                         = 1;  % From Matlab Central

%--------------------------------------------------------------------------
% For each stage function: the subroutines it uses, directly or via other
% subroutines.

stagesubroutines = struct;
stagesubroutines.EDedgeo = {'EDcheckobstrpaths','EDchkISvisiblex','EDpoinplax',...
    'EDcoordtrans1','EDinfrontofplane'};
stagesubroutines.EDSorRgeo = {'EDcheckobstr_pointtoedge','EDpoinpla','EDcompress3',...
    'EDcoordtrans1','EDgetedgepoints','EDinfrontofplane'};
stagesubroutines.EDfindconvexGApaths = {'EDchkISvisible','EDpoinpla','EDfindis'};
stagesubroutines.EDmakefirstordertfs = {'EDpistonedgeintegrand','EDwedge1st_fd',...
    'EDbetaoverml_fd','EDquadgk','EDgetGLrule','EDdistelements','lgwt'};
stagesubroutines.EDmakefirstorderirs = {'EDcoordtrans2','EDcross','EDwedge1st_ir',...
    'EDbetaoverml','EDquadstep'};
stagesubroutines.EDed2geo = {'EDcheckobstr_edgetoedge','EDpoinpla','EDcompress7',...
    'EDcoordtrans2','EDcross','EDgetedgepoints','EDinfrontofplane'};
stagesubroutines.EDinteg_submatrixstructure = {'EDgetGLrule','EDdistelements','lgwt'};
stagesubroutines.EDintegralequation_convex_tf = {'EDinteg_souterm','EDintegratebetaoverm',...
    'EDbetaoverm_fd','EDcalcpropagatematrix','EDbuildHsparse','EDmatvecH',...
    'EDcalcedgeinteqmatrixsub2_mex','EDcoordtrans1','EDdistelements','lgwt'};
stagesubroutines.EDfindHODpaths = {};
stagesubroutines.EDmakeHODirs = {'EDwedge2nd','EDwedgeN','EDcreindexmatrix','EDcoordtrans1'};
stagesubroutines.EDpropagateESIEBEM = {'EDcalccosfi','EDcalcdist'};

%--------------------------------------------------------------------------
% Functions which are called by the stage functions but do not affect any
% results.

ignored = {'DataHash','EDrecycleresultfiles','EDparenthash','EDstageversions',...
    'EDmessage','EDpostfunctext','EDgetversion','EDinputsummary','EDinteqsummary',...
    'EDinputhash','EDassembleresults'};

%--------------------------------------------------------------------------

if nargin == 0
    stageinfo = struct;
    stageinfo.versions = versions;
    stageinfo.stagesubroutines = stagesubroutines;
    stageinfo.ignored = ignored;
    return
end

if ~isfield(stagesubroutines,stagename)
    error('EDstageversions:unknownstage','%s is not a known stage function',stagename)
end

stageinfo = struct;
stageinfo.stage = stagename;
stageinfo.version = versions.(stagename);
stageinfo.subroutines = struct;
sublist = stagesubroutines.(stagename);
for ii = 1:numel(sublist)
    if ~isfield(versions,sublist{ii})
        error('EDstageversions:noversion','The subroutine %s (used by %s) has no version number',sublist{ii},stagename)
    end
    stageinfo.subroutines.(sublist{ii}) = versions.(sublist{ii});
end
