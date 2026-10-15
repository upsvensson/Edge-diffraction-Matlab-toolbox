function[tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
    EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,...
	Sdata,Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters)
% EDintegralequation_convex_tf - Computes second- and higher-order diffraction TFs with the edge source integral equation.
% For each frequency, the integral equation q = q0 + H*q for the edge
% source amplitudes is set up (q0 = the source term, from EDinteg_souterm)
% and solved, either with GMRES (explicit sparse matrix, or matrix-free)
% or with a fixed number of terms of the Neumann series. The solution is
% then propagated to the receivers (EDcalcpropagatematrix).
%
% Two stage files are saved, each with its own hash:
%   <filestem>_qfinal.mat   the solution, Qfinal (and Q_firstterm), for all
%                           frequencies. It does not depend on the
%                           receivers, so it is reused when only the
%                           receivers or the propagation settings change.
%                           Saved if filehandlingparameters.saveinteqsolution
%                           = 1 and it is smaller than 4 GB.
%   <filestem>_tfinteq.mat  the result at the receivers.
%
% Input parameters:
%   Hsubmatrixdata, planedata, edgedata, edgetoedgedata, Sdata, Rdata
%                   structs, see help EDdatastructs (planedata is not
%                   used)
%   envdata         struct; the field .cair is used.
%   controlparameters  struct; the fields .frequencies, .difforder,
%                   .inteqsolutionmethod ('GMRES_fast', 'GMRES_large' or
%                   'iterative_fixeddifforder'), .inteqtargettol,
%                   .inteqgmresmaxiter, .ESIEsourcetermimplementation,
%                   .ESIEpropagationimplementation, .skipreceiverpropagation,
%                   .ngauss, .discretizationtype and .Rstart are used. For
%                   the GMRES methods, any difforder > 2 gives the full
%                   solution.
%   EDversionnumber  The version number of the EDtoolbox (stored in the
%                   stage file).
%   filehandlingparameters
%                   struct; the fields .outputdirectory, .filestem, .showtext,
%                   .suppressresultrecycling, .saveinteqsolution and
%                   .savediff2result are used.
%
% Output parameters:
%   tfinteqdiff     The second- and higher-order diffraction, of size
%                   [nfrequencies,nreceivers] if the sources are added
%                   (Sdata.doaddsources = 1), otherwise
%                   [nfrequencies,nreceivers,nsources].
%   timingdata      [1,4] the times, for the first frequency, for building
%                   the H-matrix, the source term, solving the integral
%                   equation, and the propagation to the receivers.
%   extraoutputdata  struct with the fields
%     .tfinteqdiff_nodiff2  the same as tfinteqdiff but without the second-
%                   order diffraction (only if savediff2result = 1)
%     .inteqdiagnostics  see help EDdatastructs
%   elapsedtimehodtf
%                   The time used by this function, in seconds. If an
%                   existing stage file was reused, a second value gives
%                   the time it took to compute that file.
%   existingfilename  The name of the stage file that was reused, or empty
%                   if the data were computed.
%
% See also EDdatastructs, EDstageversions
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
%    EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,...
%    Sdata,Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters);

% 31 March 2015 Introduced detailed timing, also as output parameter.
% 8 April 2015 Substantial speeding up by saving Hsubdata instead of Hsub
%              and using the accumarray function.
% 28 June 2015 Added the option of saving the edge source signals.
% 29 Jan. 2016 Introduced the input parameters calcinteq_souterms,calcinteq_edgeterms,calcinteq_propagate
% 30 Jan 2016 Changed the inteqsettings to a struct
% 29 April 2016 Fixed a small bug (niterations did not use the struct
%              correctly in a showtext case). Also corrected a bug with
%              edge-to-edge transfer between non-parallel edges.
% 2 March 2017 Introduced a swtich between mac and non-mac computers, since
%              a mex-file exists only for mac.
%              Also, corrected a use of uint8 such that more than 255 edge
%              points are possible.
% 10 March 2017 Removed the automatic saving of ivproblematic; now it will
%               happen only if ivproblematicmatrix contains non-zero values.
% 10 March 2017 Introduced doaddsources = 2 for regular source/receiver
%               matrices which exploits symmetries.
% 25 April 2017 A small change in the for loop over the Hsub matrices. A
%               list is constructed so that the iterative loop doesnt have
%               to check if the Hsub exists or not.
% 30 May 2017   Changed back what was changed on 25 April, because it was
%               wrong.
% 26 Oct. 2017 Fixed a bug that occurred for skewed edges: the r-values
%              were distributed linearly before.
% 4 Nov. 2017  Implemented the partly visible edge in the source term.
% 27 Nov. 2017 Copied from ESIE2toolbox and trimmed down + changed to the
%               use of cells for the Hsubs.
% 28 Nov. 2017 Cleaned up code a bit. Added timingdata.
% 29 Nov. 2017 Cleaned up: removed the old timingdata code. Reduced the
%              no. of displayed frequencies.
% 13 Dec. 2017 Added the sourceamplitudes to the Sinputdata struct
% 19 Jan 2018 Renamed the ESIE2...mex to ED...mex
% 14 Oct. 2026 The header was updated to the standard structure (by Claude,
% prompted by PS).
% 15 Oct. 2026 The list of called functions was replaced by a See also line.

% 10 Dec. 2017 TODO: finish the timingdata for the CASE2 : addsources
% 19 Jan 2018 The TODO above is probably done
% 23 Jan 2018 Changed to "~strcmp(computer, 'MACI64')" after hint from Jan
%             Slechta.
% 26 Jan 2018 Implemented the handling of doallSRcombinations = 0.
% 31 Jan 2018 Had not implemented sourceamplitudes for case 2.
% 9 Feb 2018 Introduced EDinputdatahash
% 28 Feb 2018 The saving of sousig files had an extra "results"
% subdirectory.
% 13 Apr 2018 Fixed a small bug on l143: thinplaneboostvec got the format
% uint8 but needs to be double. Found by Antoine.
% 3 June 2020 Fixed a bug: folder names with spaces can be handled now
% 28 Sep. 2023 Implemented version 2 of this function while maintaining
% compatibility with the old "version 1". v2 moves the check if an existing
% file can be reused inside this function. Also updated load and save to
% the function call form, which avoids problems with spaces in file names.
% 3 Oct. 2023 Introduced the EDsettingshash as input parameter.
% 6 Oct. 2023 Fixed a small error: the timing data wasn't returned when an
% existing file was reused.
% 27 Oct. 2023 Changed the input parameter list substantially.
% 30 Oct. 2023 Fine-tuned the EDinputdatahash
% 1 Nov. 2023 Fixed a small error: doaddsources and sourceamplitudes and 
% doallSRcombinations were missing the Sdata. in one location.
% 20 Nov. 2023 Implemented the piston source
% 24 Nov. 2023 Stores many gaussvalues, up to the max. needed
% 29 Nov. 2023 Adapted to the name change of the field pistongausspoints to
% pistongaussorder
% 11 Feb. 2025 Pre-calculate one higher order of quadrature points, just in
% case an odd number of quadrature points is enforced. Also, handle the
% hidden feature that the number of gauss points can be negative.
% 10 Sep. 2026 Implemented the GMRES solution approach and also the
% tracking of term decay in the iterative series. Exports an
% inteqdiagnostics struct for analysis of the calculation performance.
% Also introduced the possibility to skip the propagation step.
% 20 Sep. 2026 FOrwarded ESIEpropagationimplementation to
% EDcalcpropagatematrix.
% 6 Oct. 2026 Implemented the propdiagnostics info from
% EDcalcpropagatemtrix
% Code improvement: changed nbig to nQ, changed bigmatrix... to Qvec...
% 10 Oct. 2026 Pruning, deletions only (by Claude, prompted by PS): removed
% the 'iterative_auto' option, the diagnostics fields that were never
% filled or only served auto-stopping (normQfinal, r_hat, esttail etc),
% the copy of ESIEpropagationimplementation in the diagnostics, and dead
% code (old accumarray blocks, unused variables).
% 10 Oct. 2026 Diagnostics rewrite (by Claude, prompted by PS): the GMRES
% diagnostics are stored per frequency (and per source, if the sources are
% handled separately) instead of being overwritten; gmres_time1/2 renamed
% to assemblytime/solvetime; normQ is only allocated for the iterative
% method, with size [difforder-2,nfrequencies,ncols]; solutionmethod, nQ
% and nnzH are now set.
% 11 Oct. 2026 Bug fixes (by Claude, prompted by PS): elapsedtimehodtf is
% now the total time (t00 was reused for the section timers); propdiagnostics
% is defined also when skipreceiverpropagation = 1; for several separate
% sources (CASE 2): timingdata(2:4) are set, the sparse matrix for
% GMRES_fast is built once per frequency instead of once per source, and
% ISOU is defined also when loadinteqsousigs = 1.
% 11 Oct. 2026 Hash fix (by Claude, prompted by PS): the solution method,
% GMRES tolerance and maxiter, the two ESIE...implementation settings,
% skipreceiverpropagation, the source type/amplitudes/doaddsources/
% doallSRcombinations (and piston geometry), cair and savediff2result are
% now part of the EDinputdatahash.
% 12 Oct. 2026 Chained hash (by Claude, prompted by PS): the geometry and
% S/R coordinates in the EDinputdatahash were replaced by the hashes of
% Hsubmatrixdata, edgedata, edgetoedgedata, Sdata and Rdata, and a local
% stageversion replaces the overall EDversionnumber. For the GMRES
% methods, difforder is not hashed (any value > 2 gives the full
% solution). load/save of the _tfinteq file without eval.
% 14 Oct. 2026 The stageversion is taken from EDstageversions, and it
% also contains the versions of the subroutines that this function uses.
% 14 Oct. 2026 Solve/propagate split (by Claude, prompted by PS): the
% solution of the integral equation, Qfinal (and Q_firstterm), for all
% frequencies is saved in a stage file of its own, <filestem>_qfinal.mat,
% with its own hash (everything except the receivers and the propagation
% settings). If only the receivers, Rstart, ESIEpropagationimplementation,
% skipreceiverpropagation, doallSRcombinations or savediff2result are
% changed, Qfinal is taken from that file, and only the propagation to the
% receivers is computed. The file is saved if
% filehandlingparameters.saveinteqsolution = 1 (default) and it is smaller
% than 4 GB. This replaces the per-frequency _sousigs files and the
% parameters saveinteqsousigs/loadinteqsousigs (which had no hash check).
% When the _tfinteq file is reused, a matching _qfinal file is copied to
% this filestem too.
% 15 Oct. 2026 (by Claude, prompted by PS) The propagation to the
% receivers is done after the frequency loop: EDcalcpropagatematrix is
% called once per receiver (with frequency = 0, which gives the
% frequency-independent part and the path lengths), instead of once per
% frequency, receiver and source. Only exp(-1i*k*(L-Rstart)) is computed
% per frequency. The same results up to rounding errors, but much faster
% for many frequencies or sources.
% 15 Oct. 2026 (by Claude, prompted by PS) The H-matrix is built once,
% before the frequency loop, without the factor exp(-1i*k*m), which
% depends only on the column: H(k) = B*diag(dcol(k)). Per frequency, only
% dcol = exp(-1i*k*mcol) is computed, and H*x = B*(dcol.*x). For
% GMRES_fast, the sparse matrix is also built only once. The same results
% up to rounding errors. New field inteqdiagnostics.maxcolumndistancediff
% (a check, should be 0).
% 15 Oct. 2026 (by Claude, prompted by PS) The propagation to the
% receivers can run on parallel workers (parfor over the receivers), with
% at most filehandlingparameters.parallelworkers workers (default 0 = no
% parallel computing). The results are the same.
% If the Parallel Computing Toolbox is not available, a note is given and
% the loop runs serially. Qfirsttermstore is sent to the workers only if
% savediff2result = 1.

t00_total = clock;
t00 = clock;
showtext = filehandlingparameters.showtext;

if ~isfield(controlparameters,'skipreceiverpropagation')
    controlparameters.skipreceiverpropagation = 0;
end

% The version of this function, and of the subroutines it uses, are
% defined in EDstageversions. They are part of the hash, so that saved
% stage files are not reused after a change in the computations.
stageversion = EDstageversions('EDintegralequation_convex_tf');

% Two hashes:
%   EDinputdatahash_solve   for the solution of the integral equation,
%                           Qfinal (the _qfinal file). It does not depend
%                           on the receivers.
%   EDinputdatahash         for the final result, at the receivers (the
%                           _tfinteq file) = the solve hash + the receivers
%                           and the propagation settings.
%
% For the GMRES methods, the full solution is computed for any
% difforder > 2, so the value of difforder does not affect the result.
if strncmp(controlparameters.inteqsolutionmethod,'GMRES',5)
    hashdifforder = 'full';
else
    hashdifforder = controlparameters.difforder;
end
solvehashstruct = struct('stage','EDintegralequation_convex_tf_solve','stageversion',stageversion,...
    'submatrixhash',EDparenthash(Hsubmatrixdata),'edgedatahash',EDparenthash(edgedata),...
    'ed2datahash',EDparenthash(edgetoedgedata),'Sdatahash',EDparenthash(Sdata),...
    'difforder',hashdifforder,'frequencies',controlparameters.frequencies,...
    'ngauss',controlparameters.ngauss,'discretizationtype',controlparameters.discretizationtype,...
    'inteqsolutionmethod',controlparameters.inteqsolutionmethod,...
    'inteqtargettol',controlparameters.inteqtargettol,...
    'inteqgmresmaxiter',controlparameters.inteqgmresmaxiter,...
    'ESIEsourcetermimplementation',controlparameters.ESIEsourcetermimplementation,...
    'sourceamplitudes',Sdata.sourceamplitudes,'doaddsources',Sdata.doaddsources,...
    'cair',envdata.cair);
EDinputdatahash_solve = DataHash(solvehashstruct);

EDinputdatastruct = struct('stage','EDintegralequation_convex_tf','stageversion',stageversion,...
    'solvehash',EDinputdatahash_solve,'Rdatahash',EDparenthash(Rdata),...
    'Rstart',controlparameters.Rstart,...
    'ESIEpropagationimplementation',controlparameters.ESIEpropagationimplementation,...
    'skipreceiverpropagation',controlparameters.skipreceiverpropagation,...
    'doallSRcombinations',Sdata.doallSRcombinations,...
    'savediff2result',filehandlingparameters.savediff2result);
EDinputdatahash = DataHash(EDinputdatastruct);

%---------------------------------------------------------------
% Sort out the file business: can an existing file be used?
% Then copy the existing file to a new copy. 

if filehandlingparameters.suppressresultrecycling == 1
    foundmatch = 0;
    existingfilename = '';
else
    [foundmatch,existingfilename] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_tfinteq',EDinputdatahash);
end

EDinputdatahash_tfinteq = EDinputdatahash;

desiredname = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'];

if foundmatch == 1
    filedata = load(existingfilename,'tfinteqdiff','extraoutputdata','timingdata','elapsedtimehodtf');
    tfinteqdiff = filedata.tfinteqdiff;
    extraoutputdata = filedata.extraoutputdata;
    timingdata = filedata.timingdata;
    if ~strcmp(existingfilename,desiredname)
        copyfile(existingfilename,desiredname);
    end
    % If a matching _qfinal file exists (possibly with another filestem),
    % it is also copied, so that this filestem has the complete set of
    % stage files.
    desiredname_qfinal = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_qfinal.mat'];
    [qfinalfound,existingqfinalfile] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_qfinal',EDinputdatahash_solve);
    if qfinalfound == 1 && ~strcmp(existingqfinalfile,desiredname_qfinal)
        copyfile(existingqfinalfile,desiredname_qfinal);
    end
    elapsedtimehodtf_new = etime(clock,t00_total);
    elapsedtimehodtf = [elapsedtimehodtf_new filedata.elapsedtimehodtf];
    return
end

extraoutputdata = struct('tfinteqdiff_nodiff2',[]);

nfrequencies = length(controlparameters.frequencies);

timingdata = zeros(1,4);

propdiagnostics = [];

inteqdiagnostics = edinteqdiagnosticstemplate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect special cases: check if the object is thin
    
nsources = size(Sdata.coordinates,1);
nreceivers = size(Rdata.coordinates,1);

% 11 Feb. 2025 by PS Make sure that one higher quadrature order than before
% has its weights and nodes pre-calculated.
maxngauss = abs(controlparameters.ngauss);
if controlparameters.ngauss < 0
    maxngauss = maxngauss +2;
end
gaussvectors = cell(maxngauss,1);
for ii = 2:maxngauss
    [gaussvec1,gaussvec2] = EDdistelements(ii,controlparameters.discretizationtype);
    gaussvectors{ii} = [gaussvec1 gaussvec2];
end

ivproblematicmatrix = zeros(nreceivers,100);

listofactiveHsubs = zeros(size(Hsubmatrixdata.listofsubmatrices,1),1);
Hsubcounter = 1;

alltripletsthin = prod(double(Hsubmatrixdata.isthinplanetriplet));
notripletsthin = prod(1-double(Hsubmatrixdata.isthinplanetriplet));
 nQ = Hsubmatrixdata.Qvecendnums(end);
nsubmatrices = size(Hsubmatrixdata.listofsubmatrices,1);

refto = edgetoedgedata.reftoshortlistE;

if isempty(Hsubmatrixdata.reftoshortlist)
    symmetrycompression = 0;
else    
    symmetrycompression = 1;
end

if symmetrycompression < 1
   reftoshortlist = 1:size(Hsubmatrixdata.listofsubmatrices,1);
end

Hsubconstant = -1/4/pi/2;

if alltripletsthin == 1
    thinplaneboostvec = 2*ones(nsubmatrices,1);
else
    if notripletsthin == 1
        thinplaneboostvec = ones(nsubmatrices,1);
    else
        thinplaneboostvec = double(Hsubmatrixdata.isthinplanetriplet(1:nsubmatrices)+1);
    end
end

Hsubdatalists = cell(nsubmatrices,1);
ivuselists = cell(nsubmatrices,1);

isthinhole = 0;

% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% We must have an outer loop over frequency, since all matrices change with
% frequency. The for-loop stretches to the end of this function.

if showtext >= 1
    freqstep = round(nfrequencies/10);
    if freqstep == 0
       freqstep = 1; 
    end
end

inteqdiagnostics.solutionmethod = controlparameters.inteqsolutionmethod;
inteqdiagnostics.nQ = double(nQ);
% The diagnostics get one column per source if the sources are handled
% separately (CASE 2 below), otherwise a single column.
if Sdata.doaddsources == 1 || nsources == 1
    ndiagcols = 1;
else
    ndiagcols = nsources;
end
if strncmp(controlparameters.inteqsolutionmethod,'GMRES',5)
    inteqdiagnostics.gmresflag   = NaN(nfrequencies,ndiagcols);
    inteqdiagnostics.gmresrelres = NaN(nfrequencies,ndiagcols);
    inteqdiagnostics.gmresiter   = NaN(nfrequencies,ndiagcols);
    inteqdiagnostics.solvetime   = NaN(nfrequencies,ndiagcols);
    if strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast')
        inteqdiagnostics.assemblytime = NaN(nfrequencies,1);
    end
else
    inteqdiagnostics.normQ = zeros(max(controlparameters.difforder-2,0),nfrequencies,ndiagcols);
end

%---------------------------------------------------------------
% Can an existing solution of the integral equation (the _qfinal file) be
% used? Then the H-matrix is not built and the integral equation is not
% solved; only the propagation to the receivers is computed.

desiredname_qfinal = [filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_qfinal.mat'];
if filehandlingparameters.suppressresultrecycling == 1
    qfinalreused = 0;
    existingqfinalfile = '';
else
    [qfinalreused,existingqfinalfile] = EDrecycleresultfiles(filehandlingparameters.outputdirectory,'_qfinal',EDinputdatahash_solve);
end
if qfinalreused == 1
    qfiledata = load(existingqfinalfile,'Qfinalstore','Qfirsttermstore','inteqdiagnostics_solve','timingdata_solve');
    Qfinalstore = qfiledata.Qfinalstore;
    Qfirsttermstore = qfiledata.Qfirsttermstore;
    inteqdiagnostics = qfiledata.inteqdiagnostics_solve;
    timingdata_solve = qfiledata.timingdata_solve;
    clear qfiledata
    if ~strcmp(existingqfinalfile,desiredname_qfinal)
        copyfile(existingqfinalfile,desiredname_qfinal);
    end
    if showtext >= 1
        [~,qfilestem] = fileparts(existingqfinalfile);
        disp(['      Using the integral equation solution in ',qfilestem,'.mat'])
    end
else
    Qfinalstore = cell(nfrequencies,ndiagcols);
    Qfirsttermstore = cell(nfrequencies,ndiagcols);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1, before the frequency loop: the frequency-independent part of the
% H-matrix.
%
% Each element of an H-submatrix is
%   (2-acrossface_out)*Hsubconstant*ny*beta*exp(-1i*k*m)/m*dz*thinplaneboost
% where m is the distance from the node on the sending edge (edge1) to the
% node on the edge in between (edge2). Only exp(-1i*k*m) depends on the
% frequency, and m depends only on the column of H (a node pair of the
% sending edge pair). Therefore
%   H(k) = B*diag(dcol(k)),  dcol(k) = exp(-1i*k*mcol)
% where B is real and frequency-independent. B is built once, here (in
% Hsubdatalists, and for GMRES_fast also as the sparse matrix Bsp), and
% for each frequency only dcol is computed: H*x = B*(dcol.*x).
% Before 15 Oct. 2026, the H-matrix was built for every frequency, with the
% same result up to rounding errors.

% PART 1: Build the H-matrix

if showtext >= 2 && controlparameters.difforder > 2
    disp(' ')
    disp(['      Building the H matrix of size (',int2str( nQ),',',int2str( nQ),') with ',int2str(Hsubmatrixdata.submatrixcounter),' submatrices'])
end

t00 = clock;

% In the for-loop, nn will be the Hsub number
% The ii value will refer to the row-number in the edgepairlist for the 
% receiving edge-pair.
% The jj value will refer to the row-number in the edgepairlist for the 
% sending edge-pair.

% mcol: the distance m (from the node on the sending edge to the node
% on the edge in between) for each column of H. maxmcoldiff checks
% that all submatrices give the same value for a column.
mcol = zeros(nQ,1);
mcolset = false(nQ,1);
maxmcoldiff = 0;

n1prev = 0;
n2prev = 0;
n3prev = 0;
shortlistprev = 0;

% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% If difforder > 2, create the big H-matrix, in the form of separate sub-matrices
% The H-matrix is independent of sources and receivers.

if controlparameters.difforder > 2 && qfinalreused == 0
    for nn = 1:nsubmatrices

        edge1 = Hsubmatrixdata.edgetripletlist(nn,3);
        n1 = Hsubmatrixdata.nedgeelems(edge1);

        edge2 = Hsubmatrixdata.edgetripletlist(nn,2);                   
        n2 = Hsubmatrixdata.nedgeelems(edge2);
           len2 = edgedata.edgelengthvec(edge2);

        edge3 = Hsubmatrixdata.edgetripletlist(nn,1);
        n3 = Hsubmatrixdata.nedgeelems(edge3);
%             disp(['triplet no. ',int2str(nn),': n1,n2,n3 = ',int2str(n1),',',int2str(n2),',',int2str(n3)])
        shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);

%             ind_start = 1 + (Hsubmatrixdata.Qvecstartnums(ii)-1) + (Hsubmatrixdata.Qvecstartnums(jj)-1)* nQ;

       if showtext >= 3
           disp(' ')
            disp(['         H-submatrix ',int2str(nn),': From edge ',int2str(edge1),' via edge ',int2str(edge2),' to edge ',int2str(edge3)])
            disp(['         of size ',int2str(n3*n2),' by ',int2str(n2*n1),'. Thinplaneboost = ',int2str(thinplaneboostvec(nn))])
            if symmetrycompression ==1
                disp(['         stored as shortlist number ',int2str(shortlistnumber)])
            end
       end

       % Build vertical [n3,n2] matrices and expand them horizontally
       % NB! n1 refers to edge 1, which is the source edge
       %     n3 refers to edge 3, which is the target edge
       % edgetriplets written in order [edge3 edge2 edge 1]

       if n1 ~= n1prev || n2 ~= n2prev || n3 ~= n3prev

%                if computer ~= 'MACI64'
           if ~strcmp(computer, 'MACI64')
               if n3 <= 255
                   n3vertmat = uint8(1:n3);
               else
                   n3vertmat = uint16(1:n3);                   
               end
               n3vertmat = reshape(repmat(n3vertmat,n2,1),n2*n3,1);           
               n3vertmat = n3vertmat(:,ones(1,n1*n2));

               if n2 <= 255
                   n2vertmat = uint8(1:n2).';
               else
                  n2vertmat = uint16(1:n2).'; 
               end
               n2vertmat = repmat(n2vertmat,n3,1);
               n2vertmat = n2vertmat(:,ones(1,n1*n2));

           % Build horizontal n2 and n1 matrices and expand them vertically
                if n1 <= 255                    
                   n1hormat = uint8(1:n1);
                else
                    n1hormat = uint16(1:n1);
                end
               n1hormat = repmat(n1hormat,1,n2);
                   n1hormat = n1hormat(ones(n3*n2,1),:);

                if n2 <= 255
                   n2hormat = uint8(1:n2);
                else
                   n2hormat = uint16(1:n2); 
                end
                   n2hormat = reshape(repmat(n2hormat,n1,1),1,n2*n1);           
               n2hormat = n2hormat(ones(n3*n2,1),:);

           % Now, we should use only the combinations where n2vertmat ==
           % n2hormat! Otherwise we get too many combinations: there should
           % be n1*n2*n3 unique combinations, but the Hsub-matrix has size
           % [n3*n2,n2*n1], which has n2 too many elements!
           %
           % Once we have the iv, which identifies the wanted matrix
           % elements, we can throw away n2hormat because it is identical
           % to n2vertmat in the wanted matrix elements!

               ivuse = find(n2hormat==n2vertmat);
               n3vertmat = n3vertmat(ivuse);
               n2vertmat = n2vertmat(ivuse);
               n1hormat  = n1hormat(ivuse);

           else
%                     [ivuse,n1hormat,n2vertmat,n3vertmat] = ESIE2calcedgeinteqmatrixsub2_mex(n1,n2,n3);       
                [ivuse,n1hormat,n2vertmat,n3vertmat] = EDcalcedgeinteqmatrixsub2_mex(n1,n2,n3);       

           end

           ivuse = uint64(ivuse);

       end

       % Build vectors that will be referred to by the n1,n2,n3
       % matrices

       if n1 ~= n1prev
            % n1vec = full(Hsubmatrixdata.quadraturematrix_pos(n1,1:n1));
            % n1vec = n1vec(:);
            n1vec = full(Hsubmatrixdata.quadraturematrix_pos(1:n1,n1));
            % weightvec = full(Hsubmatrixdata.quadraturematrix_weights(n1,1:n1));
            % weightvec = weightvec(:);
            weightvec = full(Hsubmatrixdata.quadraturematrix_weights(1:n1,n1));
            n1prev = n1;
       end       
       dzvec = weightvec*edgedata.edgelengthvec(edge1);

        if n2 ~= n2prev
            % n2vec = full(Hsubmatrixdata.quadraturematrix_pos(n2,1:n2));
            % n2vec = n2vec(:);
            n2vec = full(Hsubmatrixdata.quadraturematrix_pos(1:n2,n2));
            n2prev = n2;
        end
        if n3 ~= n3prev        
            % n3vec = full(Hsubmatrixdata.quadraturematrix_pos(n3,1:n3));
            % n3vec = n3vec(:);
            n3vec = full(Hsubmatrixdata.quadraturematrix_pos(1:n3,n3));
            n3prev = n3;
        end

        % Bug found and fixed 22 Sep. 2014: lines below did not always detect if another
        % edge was skewed relative to the first one, because only the
        % thetae1sho was used. Necessary to check thetae1sho and thetae2sho.

       ny = pi/(2*pi-edgedata.closwedangvec(edge2));
       thetain  = edgetoedgedata.thetae1sho(refto(edge1,edge2));
       thetain_other_end_of_edge  = edgetoedgedata.thetae2sho(refto(edge1,edge2));
    %    if abs(thetain) < 1e-10 || abs(thetain-(2*pi-closwedangvec(edge2))) < 1e-10   % Old erroneous line
       if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_other_end_of_edge) < 1e-10
           acrossface_in = 1;       
       else
           acrossface_in = 0;
       end
       thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
       thetaout_other_end_of_edge = edgetoedgedata.thetae2sho(refto(edge3,edge2));
    %   if abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-closwedangvec(edge2))) < 1e-10    % Old erroneous line
       if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetaout-thetaout_other_end_of_edge) < 1e-10
           acrossface_out = 1;
       else
           acrossface_out = 0;
       end

       if showtext >= 3
           disp(['         Acrossface-out = ',int2str(acrossface_out),' and Acrossface-in = ',int2str(acrossface_in), ' (not used)'])
       end

       ze2_re2 = n2vec*len2;
       ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2))   );
%        re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

       if shortlistnumber ~= shortlistprev      

           n1vec = n1vec(:).';

           ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );

           % Bug fixed 26 Oct. 2017: the line below is wrong because the
           % r-values are not distributed linearly.
%            re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );


           if acrossface_in*acrossface_out == 0 % This case is a preparation for non-convex objects.
                % Reconvert the cylindrical coordinates of
                % edge1-re2 into cartesian! Then we can use the EDcoordtrans.
                % We define our own cartesian coord syst such that the reference
                % edge (edge2) has its starting point in [0 0 0], the z-axis along its
                % own edge, and the x-axis is along (onto) the reference plane.

                xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
                ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);

                [re1_re2,thetae1_re2,~] = EDcoordtrans1([xe1_re2.' ye1_re2.' ze1_re2.'],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
                re1_re2 = re1_re2(:).';
                thetae1_re2 = thetae1_re2(:).';

                xedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
                xedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
                yedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
                yedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
                xe3_re2 = xedge3start + n3vec.*(xedge3end-xedge3start);
                ye3_re2 = yedge3start + n3vec.*(yedge3end-yedge3start);

                [re3_re2,thetae3_re2,~] = EDcoordtrans1([xe3_re2 ye3_re2 ze3_re2],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
                re3_re2 = re3_re2(:);
                thetae3_re2 = thetae3_re2(:);

           else

                re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );
                re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

           end

           ldist = (ze3_re2(n3vertmat) - ze2_re2(n2vertmat)).^2   + re3_re2(n3vertmat).^2;
           ldist = sqrt(ldist);

           mdist = (ze1_re2(n1hormat).' - ze2_re2(n2vertmat)).^2   + (re1_re2(n1hormat).').^2;
           mdist = sqrt(mdist);

           % The factor exp(-1i*k*mdist) is not included here: it is
           % applied per frequency, as dcol (see above).
           dz_overm = dzvec(n1hormat)./mdist;

           ch = ((ze2_re2(n2vertmat) - ze1_re2(n1hormat).').*(ze2_re2(n2vertmat) - ze3_re2(n3vertmat)) + mdist.*ldist)./re1_re2(n1hormat).'./re3_re2(n3vertmat);
           ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
           ch = ( ch + 1./ch)/2;

           % PS 28 March 2011: The 4 beta terms are identical!       
           if acrossface_in == 1 && acrossface_out == 1
               beta = 4*( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))));
           else
               beta = ( sin(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                      ( sin(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                      ( sin(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                      ( sin(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10));
           end

           if ismember(edge1,edgedata.offedges)==1 || ismember(edge2,edgedata.offedges)==1 || ismember(edge3,edgedata.offedges)==1
               beta = beta*0;
           else
    %           disp(['Edge1: ',int2str(edge1),' edge2: ',int2str(edge2),' edge3: ',int2str(edge3)])
           end

            Hsubdatalists{nn} = (2-acrossface_out)*Hsubconstant*ny*beta.*dz_overm*thinplaneboostvec(nn);   
            ivuselists{nn} = ivuse;
            listofactiveHsubs(Hsubcounter) = nn;
            Hsubcounter = Hsubcounter + 1;

       end

        % The distance mdist for each column of this submatrix, stored
        % at the global column numbers. A column is a node pair of the
        % sending edge pair, so mdist depends only on the column. For a
        % submatrix that reuses the data of the previous one (symmetry
        % compression), ivuse and mdist are those of that submatrix.
        iiout = Hsubmatrixdata.listofsubmatrices(nn,2);
        jjin  = Hsubmatrixdata.listofsubmatrices(nn,3);
        nrowsHsub = double(Hsubmatrixdata.Qvecendnums(iiout)) - double(Hsubmatrixdata.Qvecstartnums(iiout)) + 1;
        ivusedouble = double(ivuse);
        rowHsub = rem(ivusedouble-1,nrowsHsub) + 1;
        globalcol = double(Hsubmatrixdata.Qvecstartnums(jjin)) - 1 + (ivusedouble-rowHsub)/nrowsHsub + 1;
        alreadyset = mcolset(globalcol);
        if any(alreadyset)
            maxmcoldiff = max([maxmcoldiff; abs(mcol(globalcol(alreadyset)) - mdist(alreadyset))]);
        end
        mcol(globalcol) = mdist;
        mcolset(globalcol) = true;

        shortlistprev = shortlistnumber;

    end
end

ivactive = listofactiveHsubs~=0;
listofactiveHsubs = listofactiveHsubs(ivactive);

if maxmcoldiff > 1e-9*max([mcol;1])
    warning('EDintegralequation_convex_tf:mcol',...
        'The column distances of the H-matrix differ by %g between submatrices (should be 0)',maxmcoldiff)
end
inteqdiagnostics.maxcolumndistancediff = maxmcoldiff;

% For GMRES_fast, the sparse matrix B is built once; per frequency,
% H*x = B*(dcol.*x).
if strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast') == 1 && ...
        controlparameters.difforder > 2 && qfinalreused == 0
    t00_gmres = clock;
    Bsp = EDbuildHsparse(Hsubmatrixdata,Hsubdatalists,ivuselists);
    inteqdiagnostics.assemblytime(1,1) = etime(clock,t00_gmres);
    inteqdiagnostics.nnzH = nnz(Bsp);
end

timingdata(1) = etime(clock,t00);
t00 = clock;

% savetemp
% pause

for ifreq = 1:nfrequencies     
    frequency = controlparameters.frequencies(ifreq);
    if showtext >= 1
        if round(ifreq/freqstep)*freqstep == ifreq
            disp(['      Frequency no. ',int2str(ifreq),' (of ',int2str(nfrequencies),'): ',num2str(frequency),' Hz'])
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k = 2*pi*frequency/envdata.cair;

    % PART 1: the frequency-dependent factor of the H-matrix, which was
    % built (without this factor) before the frequency loop:
    % H(k) = B*diag(dcol).
    if controlparameters.difforder > 2 && qfinalreused == 0
        dcol = exp(-1i*k*mcol);
    end
    if ifreq == 1
        t00 = clock;
    end

    % End of PART 1: All the H submatrices have been built (for one freq.)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASE 1: A single source, or many sources' contributions should be
    % added to give the the edge source signals
    % (doaddsources = 1, or we have only one source)
    % PART 2: Create Q_firstterm 
    
    if Sdata.doaddsources == 1 || nsources == 1
        if ifreq==1
            P_receiver = zeros(nfrequencies,nreceivers);
        end
        
        if qfinalreused == 0
                        
            Q_firstterm = 0;
            for isou = 1:nsources
                ISOU = int2str(isou);

                if strcmp(Sdata.sourcetype,'monopole') == 1
                     Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                        controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                        Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
                        Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),filehandlingparameters.showtext);
                elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
                    Q_firstterm_sumpistongaussorder = 0;
                    for np = 1:size(Sdata.pistongausscoordinates{1},1)
                        Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                            controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                            Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),filehandlingparameters.showtext);
                        Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
                    end
                    Q_firstterm_addition = Q_firstterm_sumpistongaussorder;
                end
%                 Q_firstterm = Q_firstterm + Q_firstterm_addition;
                Q_firstterm = Q_firstterm + Q_firstterm_addition*Sdata.sourceamplitudes(isou,ifreq);
            end
        
            if ifreq == 1
                timingdata(2) = etime(clock,t00);
                t00 = clock;
            end
            
            % End of CASE 1, PART 2 
            % Q_firstterm has been computed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CASE 1,PART 3: Starting from Q_firstterm, iterate to find Qfinal.
            %
            % Two for-loops: iterations (iiter) and submatrices (nn)            
        
            Q = Q_firstterm;
            Qfinal = Q;   
            if strcmp(controlparameters.inteqsolutionmethod,'iterative_fixeddifforder') == 1
                if showtext >= 2
                    disp('      Starting to solve the int.eq. with the iterative method')
                end

                for iiter = 1:(controlparameters.difforder-2)
                    if showtext >= 3
                       disp(['         Iteration no. ',int2str(iiter),' of ',int2str(controlparameters.difforder-2)])
                    end
                    Q = EDmatvecH(dcol.*Q, Hsubmatrixdata, Hsubdatalists, ivuselists);
                    Qfinal = Qfinal + Q; 
                    inteqdiagnostics.normQ(iiter,ifreq,1) = norm(Q);
                end

            elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast') == 1 % Here we do the GMRES, fast version
               if showtext >= 2
                    disp('      Starting to solve the int.eq. with GMRES, fast version (for smaller problems)')
                end
                  A = @(x) x - Bsp*(dcol.*x);     % H*x = B*(dcol.*x)
                  t00_gmres = clock;    
                  [Qfinal,flag,relres,iter] = gmres(A,Q_firstterm,[],controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
                  inteqdiagnostics.solvetime(ifreq,1) = etime(clock,t00_gmres);
                  inteqdiagnostics.gmresiter(ifreq,1) = iter(2);
                  inteqdiagnostics.gmresflag(ifreq,1) = flag;
                  inteqdiagnostics.gmresrelres(ifreq,1) = relres;
                  
            elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_large') == 1 % Here we do the GMRES for large problems

               if showtext >= 2
                    disp('      Starting to solve the int.eq. with GMRES, the version for larger problems')
               end
               t00_gmres = clock;    
               afun = @(q) q - EDmatvecH(dcol.*q, Hsubmatrixdata, Hsubdatalists, ivuselists);
               [Qfinal, flag, relres, iter] = gmres(afun, Q_firstterm, [], controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
               inteqdiagnostics.solvetime(ifreq,1) = etime(clock,t00_gmres);
               inteqdiagnostics.gmresiter(ifreq,1) = iter(2);
               inteqdiagnostics.gmresflag(ifreq,1) = flag;
               inteqdiagnostics.gmresrelres(ifreq,1) = relres;

            end

            doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
            for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
               sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter)); 
               doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
            end        

            if ifreq == 1
                timingdata(3) = etime(clock,t00);
                t00 = clock;
            end

            % End of CASE 1, PART 3 
            % Qfinal has been computed. It is stored, for the _qfinal file.
            
            Qfinalstore{ifreq,1} = Qfinal;
            Qfirsttermstore{ifreq,1} = Q_firstterm;

        else    % qfinalreused = 1: the solution is taken from the _qfinal file
            Qfinal = Qfinalstore{ifreq,1};
            Q_firstterm = Qfirsttermstore{ifreq,1};
            doesQsegmenthavevalues = edQsegmentflags(Qfinal,Hsubmatrixdata);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CASE 1, PART 4: the propagation of Qfinal to the receivers is
        % done after the frequency loop (see PART 4 below).

        % End of CASE 1, PART 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASE 2: doaddsources = 0  and we have several sources
    % We might also have the case that doallSRcombinations = 0.
    % 
    % PART 2: Create Q_firstterm for each source
    % There will be an outer loop over sources
         
    elseif Sdata.doaddsources == 0   % This means we have several sources, and they should not be added
        if ifreq==1
            P_receiver = zeros(nfrequencies,nreceivers,nsources);
        end

        % For GMRES_fast: the sparse matrix B was built before the
        % frequency loop; H*x = B*(dcol.*x).
        if strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast') == 1 && qfinalreused == 0
            A = @(x) x - Bsp*(dcol.*x);
        end

        for isou = 1:nsources
            ISOU = int2str(isou);
            
            if qfinalreused == 0
                
                ISOU = int2str(isou); 
    
                if showtext >= 2
                    disp(['Source no. ',ISOU,' of ',int2str(nsources)])
                end
                
                if strcmp(Sdata.sourcetype,'monopole') == 1
                    Q_firstterm = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                        controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                        Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),...
                        Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),filehandlingparameters.showtext);
                elseif strcmp(Sdata.sourcetype,'polygonpiston') == 1
                    Q_firstterm_sumpistongaussorder = 0;
                    for np = 1:size(Sdata.pistongausscoordinates{1},1)
                            Q_firstterm_addition = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...                    
                            controlparameters,Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),...
                            Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),filehandlingparameters.showtext);
                        Q_firstterm_sumpistongaussorder = Q_firstterm_sumpistongaussorder + Q_firstterm_addition*squeeze(Sdata.pistongaussweights{isou}(np,1));
                    end
                    Q_firstterm = Q_firstterm_sumpistongaussorder;

                end

                Q_firstterm = Q_firstterm*Sdata.sourceamplitudes(isou,ifreq);

                if ifreq == 1 && isou == 1
                    timingdata(2) = etime(clock,t00);
                    t00 = clock;
                end
                              
                % End of CASE 2, PART 2 
                % Q_firstterm has been computed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CASE 2,PART 3: Starting from Q_firstterm, iterate to find Qfinal.
                %
                % Two for-loops: iterations (iiter) and submatrices (nn)
                  

                 Q = Q_firstterm;
                Qfinal = Q;   
                
                if strcmp(controlparameters.inteqsolutionmethod,'iterative_fixeddifforder') == 1
                   if showtext >= 2
                        disp('      Starting to solve the int.eq. with the iterative method')
                    end

                    for iiter = 1:(controlparameters.difforder-2)
                       if showtext >= 3
                           disp(['         Iteration no. ',int2str(iiter),' of ',int2str(controlparameters.difforder-2)])
                       end
                       Q = EDmatvecH(dcol.*Q, Hsubmatrixdata, Hsubdatalists, ivuselists);
                        Qfinal = Qfinal + Q; 

                        inteqdiagnostics.normQ(iiter,ifreq,isou) = norm(Q);
                        
                    end % for iiter = ....
 
                elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast') == 1 
                   if showtext >= 2
                        disp('      Starting to solve the int.eq. with GMRES, fast version (for smaller problems)')
                   end
                     t00_gmres = clock;
                     [Qfinal,flag,relres,iter] = gmres(A,Q_firstterm,[],controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
                     inteqdiagnostics.solvetime(ifreq,isou) = etime(clock,t00_gmres);                    
                     inteqdiagnostics.gmresiter(ifreq,isou) = iter(2);
                     inteqdiagnostics.gmresflag(ifreq,isou) = flag;
                     inteqdiagnostics.gmresrelres(ifreq,isou) = relres;
                    
                elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_large') == 1 
                   if showtext >= 2
                        disp('      Starting to solve the int.eq. with GMRES, the version for larger problems')
                    end
                    t00_gmres = clock;
                    afun = @(q) q - EDmatvecH(dcol.*q, Hsubmatrixdata, Hsubdatalists, ivuselists);
                    [Qfinal, flag, relres, iter] = gmres(afun, Q_firstterm, [], controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
                    inteqdiagnostics.solvetime(ifreq,isou) = etime(clock,t00_gmres);                     
                    inteqdiagnostics.gmresiter(ifreq,isou) = iter(2);
                    inteqdiagnostics.gmresflag(ifreq,isou) = flag;
                    inteqdiagnostics.gmresrelres(ifreq,isou) = relres;

                end

                % savetemp
                % pause
                % 

                doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
                for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
                   sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter)); 
                   doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
                end

                if ifreq == 1 && isou == 1
                    timingdata(3) = timingdata(3) + etime(clock,t00);
                    t00 = clock;
                end

                % End of CASE 2, PART 3 
                % Qfinal has been computed. It is stored, for the _qfinal
                % file.

                Qfinalstore{ifreq,isou} = Qfinal;
                Qfirsttermstore{ifreq,isou} = Q_firstterm;
            
            else    % qfinalreused = 1: the solution is taken from the _qfinal file
                Qfinal = Qfinalstore{ifreq,isou};
                Q_firstterm = Qfirsttermstore{ifreq,isou};
                doesQsegmenthavevalues = edQsegmentflags(Qfinal,Hsubmatrixdata);
            end  

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CASE 2, PART 4: the propagation of Qfinal to the receivers is
            % done after the frequency loop (see PART 4 below).
        end   % for isou = 1:nsources

    end  % if Sdata.doaddsources == 1 || nsources == 1
end  % for ifreq = 1:nfrequencies     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4: Propagate Qfinal to the receivers, for all frequencies (and all
% sources, if they are handled separately).
%
% The propagation row vector depends on the frequency only through the
% factor exp(-1i*k*(L-Rstart)), where L = l + m is the path length via the
% two edge nodes to the receiver. Everything else (the beta terms, the
% singularity treatment, 1/(l*m) and the quadrature weights) depends only
% on the geometry. EDcalcpropagatematrix is therefore called once per
% receiver, with frequency = 0, which gives the frequency-independent part
% Fgeom and the path lengths. For each frequency, only the exponential
% factor is computed.

t00_prop = clock;
% The receivers are independent of each other, so this loop can run on
% several parallel workers: filehandlingparameters.parallelworkers is the
% largest number of workers (0 = an ordinary serial loop; parfor with 0
% workers runs in this Matlab session). Every receiver gets its results in
% cell arrays, which are collected after the loop.
nworkers = 0;
if isfield(filehandlingparameters,'parallelworkers') && ~isempty(filehandlingparameters.parallelworkers)
    nworkers = filehandlingparameters.parallelworkers;
end
% Without the Parallel Computing Toolbox, parfor runs serially anyway; then
% a note is given (on the screen and in the log file).
parallelnote = '';
if nworkers > 0
    try
        hasparallel = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    catch
        hasparallel = false;
    end
    if ~hasparallel
        parallelnote = sprintf(['parallelworkers = %g, but the Parallel Computing Toolbox is not ',...
            'available: the receivers are computed serially.'],nworkers);
        nworkers = 0;
        if showtext >= 1
            disp(['      ',parallelnote])
        end
    end
end
inteqdiagnostics.parallelworkers = nworkers;
inteqdiagnostics.parallelnote = parallelnote;
if controlparameters.skipreceiverpropagation == 0
    if Sdata.doaddsources == 1 || nsources == 1
        nsolutions = 1;
    else
        nsolutions = nsources;
    end
    % The edge pairs where Qfinal has values, for some frequency and source
    doesQsegmenthavevalues = zeros(size(Hsubmatrixdata.edgepairlist,1),1);
    for ifreq = 1:nfrequencies
        for isol = 1:nsolutions
            doesQsegmenthavevalues = doesQsegmenthavevalues | edQsegmentflags(Qfinalstore{ifreq,isol},Hsubmatrixdata);
        end
    end
    doesQsegmenthavevalues = double(doesQsegmenthavevalues);
    kvec = 2*pi*controlparameters.frequencies/envdata.cair;
    savediff2 = filehandlingparameters.savediff2result;
    doallSRcombinations = Sdata.doallSRcombinations;
    Rstartvalue = controlparameters.Rstart;
    propimplementation = controlparameters.ESIEpropagationimplementation;

    % Qfirsttermstore is only needed for savediff2result = 1; otherwise an
    % empty cell array is used in the loop, so that it is not copied to the
    % parallel workers.
    if savediff2 == 1
        Qfirsttermforloop = Qfirsttermstore;
    else
        Qfirsttermforloop = {};
    end

    Pcell = cell(nreceivers,1);
    P2cell = cell(nreceivers,1);
    ivproblematiccell = cell(nreceivers,1);
    propdiagnosticscell = cell(nreceivers,1);

    parfor (irec = 1:nreceivers, nworkers)
        if showtext >=2
            disp(['      Receiver number ',int2str(irec)])
        end
        % The solutions (sources) that are propagated to this receiver
        % (doallSRcombinations = 0: only source N to receiver N)
        if nsolutions == 1
            solutionlist = 1;
        elseif doallSRcombinations == 1
            solutionlist = 1:nsolutions;
        elseif irec <= nsolutions
            solutionlist = irec;
        else
            solutionlist = [];
        end
        Prec = zeros(nfrequencies,nsolutions);
        P2rec = zeros(nfrequencies,nsolutions);
        ivproblematicrec = [];
        propdiagnosticsrec = [];

        if ~isempty(solutionlist)
            [Fgeom,ivproblematicrec,propdiagnosticsrec,pathlengths] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,...
                Hsubmatrixdata,isthinhole,Rdata.vispartedgesfromr(:,irec),0,Rstartvalue,...
                Rdata.rRsho(Rdata.reftoshortlistR(:,irec)),Rdata.thetaRsho(Rdata.reftoshortlistR(:,irec)),...
                Rdata.zRsho(Rdata.reftoshortlistR(:,irec)),doesQsegmenthavevalues,showtext,...
                propimplementation);

            % Only the non-zero elements are used
            ivnonzero = find(Fgeom);
            Fgeomnonzero = Fgeom(ivnonzero);
            pathlengthsnonzero = pathlengths(ivnonzero) - Rstartvalue;

            for ifreq = 1:nfrequencies
                Fvec = Fgeomnonzero.*exp(-1i*kvec(ifreq)*pathlengthsnonzero);
                for isol = solutionlist
                    Qfinalone = Qfinalstore{ifreq,isol};
                    Prec(ifreq,isol) = Fvec*Qfinalone(ivnonzero);
                    if savediff2 == 1
                        Qfirsttermone = Qfirsttermforloop{ifreq,isol};
                        P2rec(ifreq,isol) = Fvec*(Qfinalone(ivnonzero)-Qfirsttermone(ivnonzero));
                    end
                end
            end
        end
        Pcell{irec} = Prec;
        P2cell{irec} = P2rec;
        ivproblematiccell{irec} = ivproblematicrec;
        propdiagnosticscell{irec} = propdiagnosticsrec;
    end

    % Collect the results of all the receivers
    for irec = 1:nreceivers
        P_receiver(:,irec,:) = reshape(Pcell{irec},nfrequencies,1,nsolutions);
        if savediff2 == 1
            extraoutputdata.tfinteqdiff_nodiff2(:,irec,:) = reshape(P2cell{irec},nfrequencies,1,nsolutions);
        end
        ivproblematic = ivproblematiccell{irec};
        if ~isempty(ivproblematic)
            existinglength = size(ivproblematicmatrix,2);
            newlength = length(ivproblematic);
            if newlength > existinglength
                ivproblematicmatrix = [ivproblematicmatrix zeros(nreceivers,newlength-existinglength)]; %#ok<AGROW>
            end
            ivproblematicmatrix(irec,1:newlength) = ivproblematic;
        end
        if ~isempty(propdiagnosticscell{irec})
            propdiagnostics = propdiagnosticscell{irec};
        end
    end
    clear Pcell P2cell ivproblematiccell propdiagnosticscell
end
% timingdata(4): the propagation time per frequency, as before
timingdata(4) = etime(clock,t00_prop)/nfrequencies;

%---------------------------------------------------------------
% Save the solution of the integral equation in the _qfinal file, so that
% it can be reused for other receivers. Not saved if it would be larger
% than maxqfinalbytes.

maxqfinalbytes = 4e9;
if qfinalreused == 1
    timingdata(1:3) = timingdata_solve(1:3);
    [~,qfilestem] = fileparts(existingqfinalfile);
    qfinalfileinfo = ['reused from ',qfilestem,'.mat'];
else
    timingdata_solve = timingdata(1:3);
    qfinalbytes = 0;
    for ii = 1:numel(Qfinalstore)
        qfinalbytes = qfinalbytes + 16*(numel(Qfinalstore{ii}) + numel(Qfirsttermstore{ii}));
    end
    [~,qfilestem] = fileparts(desiredname_qfinal);
    if filehandlingparameters.saveinteqsolution ~= 1
        qfinalfileinfo = 'computed, not saved (saveinteqsolution = 0)';
    elseif qfinalbytes > maxqfinalbytes
        qfinalfileinfo = sprintf('computed, not saved (%.1f GB, more than the limit of %.1f GB)',...
            qfinalbytes/1e9,maxqfinalbytes/1e9);
    else
        qdata = struct;
        qdata.Qfinalstore = Qfinalstore;
        qdata.Qfirsttermstore = Qfirsttermstore;
        qdata.inteqdiagnostics_solve = inteqdiagnostics;
        qdata.timingdata_solve = timingdata_solve;
        qdata.EDinputdatahash = EDinputdatahash_solve;
        qdata.stageversion = stageversion;
        qdata.EDversionnumber = EDversionnumber;
        if qfinalbytes > 1.9e9
            save(desiredname_qfinal,'-struct','qdata','-v7.3');
        else
            save(desiredname_qfinal,'-struct','qdata');
        end
        clear qdata
        qfinalfileinfo = ['computed and saved in ',qfilestem,'.mat'];
    end
end
clear Qfinalstore Qfirsttermstore
inteqdiagnostics.qfinalfile = qfinalfileinfo;

elapsedtimehodtf = etime(clock,t00_total);
tfinteqdiff = P_receiver;

EDinputdatahash = EDinputdatahash_tfinteq;
inteqdiagnostics.propdiagnostics = propdiagnostics;
extraoutputdata.inteqdiagnostics = inteqdiagnostics;

save(desiredname,'tfinteqdiff','extraoutputdata','EDinputdatahash','timingdata','elapsedtimehodtf','stageversion','EDversionnumber');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inteqdiagnosticstruct = edinteqdiagnosticstemplate

inteqdiagnosticstruct = struct(...
    'solutionmethod',    [], ...
    'nQ',                NaN,...
    'nnzH',              NaN, ...
    'gmresflag',         [], ...
    'gmresrelres',       [], ...
    'gmresiter',         [], ...
    'assemblytime',      [], ...
    'solvetime',         [], ...
    'normQ',             [], ...
    'propdiagnostics',   [], ...
    'qfinalfile',        '');
end

function doesQsegmenthavevalues = edQsegmentflags(Qfinal,Hsubmatrixdata)
% For each edge pair: 1 if any of its elements in Qfinal is non-zero.
nedgepairs = length(Hsubmatrixdata.edgepairlist(:,1));
doesQsegmenthavevalues = zeros(nedgepairs,1);
for iicounter = 1:nedgepairs
    sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter));
    doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
end
end
    
end
