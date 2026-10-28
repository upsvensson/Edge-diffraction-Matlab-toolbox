function[tfinteqdiff,timingdata,extraoutputdata,elapsedtimehodtf,existingfilename] = ...
    EDintegralequation_convex_tf(Hsubmatrixdata,planedata,edgedata,edgetoedgedata,...
	Sdata,Rdata,envdata,controlparameters,EDversionnumber,filehandlingparameters)
% EDintegralequation_convex_tf - Computes second- and higher-order diffraction TFs with the edge source integral equation.
% For each frequency, the integral equation q = q0 + H*q for the edge
% source amplitudes is set up (q0 = the source term, from EDinteg_souterm)
% and solved, either with GMRES (explicit sparse matrix, or matrix-free)
% or with a fixed number of terms of the Neumann series. The solution is
% then propagated to the receivers (EDpropagateQvec). Optionally, the
% solution is first refined to more nodes per edge (EDrefineQvec, with the
% Nystrom interpolation formula), and propagated on the refined nodes
% (controlparameters.ESIEpropagationrefine).
%
% The source term (EDinteg_soutermsources), the kernel of the H-matrix
% (EDHkernelvalues) and the propagation (EDpropagateQvec) are separate
% functions, which EDrefineQvec uses too.
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
% EDGE NUMBERING (the same in all the integral equation functions)
% Each kernel is one first-order diffraction step, from an incoming point
% via a diffracting edge to an outgoing point, and the diffracting edge is
% always edge 2 (position z2):
%   source term:      q0(z3 <- z2) = G(z3 <- z2 <- S)
%   integral eq.:     q(z3 <- z2)  = q0(z3 <- z2) + int H(z3 <- z2 <- z1) q(z2 <- z1) dz1
%   propagation:      p(R)         = int int F(R <- z2 <- z1) q(z2 <- z1) dz1 dz2
% q(zb <- za) is the edge source amplitude at za (on the edge that the q
% segment is "from"), radiating towards zb (on the "to" edge). The wedge
% index ny, the beta terms and the cylindrical coordinates (the suffix
% _re2, e.g. ze1_re2) are those of edge 2. S or R takes the place of
% edge 1 or edge 3 where the kernel starts or ends at a point.
% In Hsubmatrixdata, edgepairlist(ii,:) = [to-edge from-edge] of the q
% segment ii, and edgetripletlist(nn,:) = [edge3 edge2 edge1]. Inside a q
% segment, the node index of the from-edge runs fastest: element
% ifrom + (ito-1)*nfrom.
% In this function (the H-matrix): edge1 -> edge2 -> edge3. The rows of an
% H-submatrix belong to the q segment [edge3 edge2], and the columns to
% the q segment [edge2 edge1]; edgetripletlist(nn,:) = [edge3 edge2 edge1].
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
%                   .ESIEpropagationimplementation, .ESIEpropagationrefine,
%                   .skipreceiverpropagation, .ngauss, .discretizationtype
%                   and .Rstart are used. For
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
%   timingdata      [1,5] the times for building the H-matrix (once), the
%                   source term and the solution of the integral equation
%                   (both summed over all frequencies and sources), the
%                   propagation to the receivers (all frequencies), and the
%                   refinement for the propagation (0 if it is not used).
%                   If the solution was taken from a _qfinal file, the
%                   first three values come from that file.
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
% See also EDdatastructs, EDstageversions, EDpropagateQvec, EDrefineQvec,
% EDinteg_soutermsources, EDHkernelvalues
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Oct. 2026
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
% 16 Oct. 2026 (by Claude, prompted by PS) For the 'accurate' propagation,
% the receivers are handled in blocks by the new function
% EDcalcpropagatematrix_block (all the receivers of a block at once for
% each edge pair), with parfor over the blocks. inteqdiagnostics.
% propdiagnostics then has totals over all the receivers (before: the
% values of the last receiver). The 'legacy' propagation is unchanged.
% 18 Oct. 2026 (by Claude, prompted by PS) timingdata has totals: the
% source term and the solution summed over all frequencies and sources,
% and the propagation for all frequencies (before: the first frequency,
% and the propagation per frequency).
% 18 Oct. 2026 (by Claude, prompted by PS) The edge numbering (edge 2 is
% always the diffracting edge) is described in the header.
% 28 Oct. 2026 (by Claude, prompted by PS) Refactoring, without changes of
% the results: the propagation (PART 4) was moved to the new function
% EDpropagateQvec, the source term for one or several sources (PART 2) to
% EDinteg_soutermsources, and the kernel of the H-matrix (PART 1) to
% EDHkernelvalues, so that EDrefineQvec uses the same code. The solution
% (PART 3), which was written out twice (CASE 1 and CASE 2), is in the
% local function edsolveinteq.
% 28 Oct. 2026 (by Claude, prompted by PS) New option
% controlparameters.ESIEpropagationrefine (default 1 = not used): the
% solution is refined to more nodes per edge (EDrefineQvec) before the
% propagation. It is part of the _tfinteq hash only if it is used. New
% timingdata(5), the time for the refinement, and
% inteqdiagnostics.nQpropagated, the number of propagated unknowns (only
% when the refinement is used).
% 28 Oct. 2026 The refinement uses the parallel workers too
% (filehandlingparameters.parallelworkers).
% 20 Oct. 2026 The node and weight table (gaussvectors) also covers edges
% with more points than abs(ngauss) (controlparameters.ngaussmin).

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

% The optional refinement of the nodes for the propagation
% (controlparameters.ESIEpropagationrefine). It is hashed only if it is
% used, so that the hash stays the same as before for the default value.
dorefine = edisrefinementused(controlparameters);
EDinputdatastruct = struct('stage','EDintegralequation_convex_tf','stageversion',stageversion,...
    'solvehash',EDinputdatahash_solve,'Rdatahash',EDparenthash(Rdata),...
    'Rstart',controlparameters.Rstart,...
    'ESIEpropagationimplementation',controlparameters.ESIEpropagationimplementation,...
    'skipreceiverpropagation',controlparameters.skipreceiverpropagation,...
    'doallSRcombinations',Sdata.doallSRcombinations,...
    'savediff2result',filehandlingparameters.savediff2result);
if dorefine
    EDinputdatastruct.ESIEpropagationrefine = controlparameters.ESIEpropagationrefine;
end
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

timingdata = zeros(1,5);

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
% 20 Oct. 2026: with controlparameters.ngaussmin, short edges can have more
% points than abs(ngauss)
maxngauss = max([maxngauss; double(Hsubmatrixdata.nedgeelems(:))]);
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

       % The kernel values (without the factor exp(-1i*k*m)) are computed
       % by EDHkernelvalues, which EDrefineQvec uses too. With symmetry
       % compression, a submatrix with the same shortlist number as the
       % previous one reuses its data.
       if shortlistnumber ~= shortlistprev
           [Hsubdatalists{nn},mdist] = EDHkernelvalues(edgedata,edgetoedgedata,edge1,edge2,edge3,...
               n1vec,n2vec,n3vec,dzvec,n1hormat,n2vertmat,n3vertmat,thinplaneboostvec(nn));
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

% dcol is only defined if the H-matrix is used (difforder > 2), and Bsp
% only for GMRES_fast; they are given to edsolveinteq in any case.
dcol = [];
if ~exist('Bsp','var')
    Bsp = [];
end
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
    % The times for the source term and the solution are summed over all
    % frequencies (and sources): timingdata(2) and timingdata(3).
    t00 = clock;

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
            Q_firstterm = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                controlparameters,Sdata,frequency,ifreq,gaussvectors,1:nsources,1,filehandlingparameters.showtext);

            timingdata(2) = timingdata(2) + etime(clock,t00);
            t00 = clock;

            % End of CASE 1, PART 2
            % Q_firstterm has been computed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CASE 1,PART 3: Starting from Q_firstterm, find Qfinal (the
            % local function edsolveinteq).

            [Qfinal,solvediag] = edsolveinteq(Q_firstterm,dcol,Bsp,Hsubmatrixdata,Hsubdatalists,ivuselists,...
                controlparameters,showtext);
            inteqdiagnostics = edstoresolvediag(inteqdiagnostics,solvediag,ifreq,1);

            doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
            for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
               sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter)); 
               doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
            end        

            timingdata(3) = timingdata(3) + etime(clock,t00);
            t00 = clock;

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


        for isou = 1:nsources
            ISOU = int2str(isou);
            
            if qfinalreused == 0
                
                ISOU = int2str(isou); 
    
                if showtext >= 2
                    disp(['Source no. ',ISOU,' of ',int2str(nsources)])
                end
                
                Q_firstterm = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
                    controlparameters,Sdata,frequency,ifreq,gaussvectors,isou,0,filehandlingparameters.showtext);

                timingdata(2) = timingdata(2) + etime(clock,t00);
                t00 = clock;

                % End of CASE 2, PART 2
                % Q_firstterm has been computed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CASE 2,PART 3: Starting from Q_firstterm, find Qfinal (the
                % local function edsolveinteq).

                [Qfinal,solvediag] = edsolveinteq(Q_firstterm,dcol,Bsp,Hsubmatrixdata,Hsubdatalists,ivuselists,...
                    controlparameters,showtext);
                inteqdiagnostics = edstoresolvediag(inteqdiagnostics,solvediag,ifreq,isou);

                % savetemp
                % pause
                % 

                doesQsegmenthavevalues = zeros(length(Hsubmatrixdata.edgepairlist(:,1)),1);
                for iicounter = 1:length(Hsubmatrixdata.edgepairlist(:,1))
                   sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter)); 
                   doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
                end

                timingdata(3) = timingdata(3) + etime(clock,t00);
                t00 = clock;

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
% sources, if they are handled separately), with EDpropagateQvec.
%
% Optionally (controlparameters.ESIEpropagationrefine), the solution is
% first refined to more nodes, with the Nystrom interpolation formula
% (EDrefineQvec), and propagated on the refined nodes. The solution of the
% integral equation (and the _qfinal file) is not changed; only the
% quadrature of the propagation integral is finer.

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
    propsettings = struct('Rstart',controlparameters.Rstart,...
        'implementation',controlparameters.ESIEpropagationimplementation,...
        'doallSRcombinations',Sdata.doallSRcombinations,...
        'savediff2',filehandlingparameters.savediff2result,...
        'nworkers',nworkers,'showtext',showtext,'isthinhole',isthinhole);
    if dorefine
        % The refined discretization: the number of nodes on each edge is
        % multiplied by its factor from ESIEpropagationrefine.
        t00_refine = clock;
        nedgeelems_refined = edrefinednodes(double(Hsubmatrixdata.nedgeelems(:)),...
            edgedata.edgelengthvec(:),controlparameters.ESIEpropagationrefine);
        refineinput = struct('Hsubmatrixdata',Hsubmatrixdata,'edgedata',edgedata,...
            'edgetoedgedata',edgetoedgedata,'Sdata',Sdata,'envdata',envdata,...
            'controlparameters',controlparameters);
        refineinput.Qfinalstore = Qfinalstore;
        refined = EDrefineQvec(refineinput,struct('method','nystrom','nedgeelems',nedgeelems_refined,...
            'propagate',0,'showtext',double(showtext >= 2),'parallelworkers',nworkers));
        clear refineinput
        Hsubmatrixdata_prop = refined.Hsubmatrixdata;
        Qfinalstore_prop = refined.Qfinalstore;
        if isfield(refined,'Qfirsttermstore')
            Qfirsttermstore_prop = refined.Qfirsttermstore;
        else
            Qfirsttermstore_prop = Qfirsttermstore;
        end
        clear refined
        timingdata(5) = etime(clock,t00_refine);
        inteqdiagnostics.nQpropagated = double(Hsubmatrixdata_prop.Qvecendnums(end));
    else
        Hsubmatrixdata_prop = Hsubmatrixdata;
        Qfinalstore_prop = Qfinalstore;
        Qfirsttermstore_prop = Qfirsttermstore;
        if isfield(inteqdiagnostics,'nQpropagated')
            inteqdiagnostics = rmfield(inteqdiagnostics,'nQpropagated');
        end
    end

    [P_receiver,P_nodiff2,propinfo] = EDpropagateQvec(edgedata,edgetoedgedata,Hsubmatrixdata_prop,Rdata,...
        Qfinalstore_prop,Qfirsttermstore_prop,controlparameters.frequencies,envdata,propsettings);
    clear Hsubmatrixdata_prop Qfinalstore_prop Qfirsttermstore_prop
    if filehandlingparameters.savediff2result == 1
        extraoutputdata.tfinteqdiff_nodiff2 = P_nodiff2;
    end
    propdiagnostics = propinfo.propdiagnostics;
    if ~isempty(propinfo.propagationblocksize)
        inteqdiagnostics.propagationblocksize = propinfo.propagationblocksize;
    end
    clear propinfo P_nodiff2
end
% timingdata(4): the total propagation time, for all frequencies (without
% the refinement, which is in timingdata(5))
timingdata(4) = etime(clock,t00_prop) - timingdata(5);

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

%--------------------------------------------------------------------------
function [Qfinal,solvediag] = edsolveinteq(Q_firstterm,dcol,Bsp,Hsubmatrixdata,Hsubdatalists,ivuselists,...
    controlparameters,showtext)
% Solves the integral equation Q = Q_firstterm + H*Q for one frequency
% (and one source, or the added sources), with H(k)*x = B*(dcol.*x).
% (Before 28 Oct. 2026, this code was written out twice, for CASE 1 and
% CASE 2.) solvediag has the diagnostics for this solution: .normQ for
% the iterative method, and .solvetime, .gmresiter, .gmresflag,
% .gmresrelres for the GMRES methods.

solvediag = struct;
Q = Q_firstterm;
Qfinal = Q;
if strcmp(controlparameters.inteqsolutionmethod,'iterative_fixeddifforder') == 1
    if showtext >= 2
        disp('      Starting to solve the int.eq. with the iterative method')
    end
    solvediag.normQ = zeros(max(controlparameters.difforder-2,0),1);
    for iiter = 1:(controlparameters.difforder-2)
        if showtext >= 3
           disp(['         Iteration no. ',int2str(iiter),' of ',int2str(controlparameters.difforder-2)])
        end
        Q = EDmatvecH(dcol.*Q, Hsubmatrixdata, Hsubdatalists, ivuselists);
        Qfinal = Qfinal + Q;
        solvediag.normQ(iiter) = norm(Q);
    end

elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_fast') == 1 % Here we do the GMRES, fast version
    if showtext >= 2
        disp('      Starting to solve the int.eq. with GMRES, fast version (for smaller problems)')
    end
    A = @(x) x - Bsp*(dcol.*x);     % H*x = B*(dcol.*x)
    t00_gmres = clock;
    [Qfinal,flag,relres,iter] = gmres(A,Q_firstterm,[],controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
    solvediag.solvetime = etime(clock,t00_gmres);
    solvediag.gmresiter = iter(2);
    solvediag.gmresflag = flag;
    solvediag.gmresrelres = relres;

elseif strcmp(controlparameters.inteqsolutionmethod,'GMRES_large') == 1 % Here we do the GMRES for large problems
    if showtext >= 2
        disp('      Starting to solve the int.eq. with GMRES, the version for larger problems')
    end
    t00_gmres = clock;
    afun = @(q) q - EDmatvecH(dcol.*q, Hsubmatrixdata, Hsubdatalists, ivuselists);
    [Qfinal, flag, relres, iter] = gmres(afun, Q_firstterm, [], controlparameters.inteqtargettol,controlparameters.inteqgmresmaxiter);
    solvediag.solvetime = etime(clock,t00_gmres);
    solvediag.gmresiter = iter(2);
    solvediag.gmresflag = flag;
    solvediag.gmresrelres = relres;
end
end

%--------------------------------------------------------------------------
function inteqdiagnostics = edstoresolvediag(inteqdiagnostics,solvediag,ifreq,icol)
% Stores the diagnostics of one solution (from edsolveinteq) in column icol
% of the per-frequency diagnostics.

if isfield(solvediag,'normQ')
    for iiter = 1:numel(solvediag.normQ)
        inteqdiagnostics.normQ(iiter,ifreq,icol) = solvediag.normQ(iiter);
    end
end
if isfield(solvediag,'solvetime')
    inteqdiagnostics.solvetime(ifreq,icol) = solvediag.solvetime;
    inteqdiagnostics.gmresiter(ifreq,icol) = solvediag.gmresiter;
    inteqdiagnostics.gmresflag(ifreq,icol) = solvediag.gmresflag;
    inteqdiagnostics.gmresrelres(ifreq,icol) = solvediag.gmresrelres;
end
end

%--------------------------------------------------------------------------
function dorefine = edisrefinementused(controlparameters)
% 1 if controlparameters.ESIEpropagationrefine asks for a refinement (not
% the default value 1).

dorefine = 0;
if isfield(controlparameters,'ESIEpropagationrefine')
    r = controlparameters.ESIEpropagationrefine;
    if isnumeric(r) && ~isempty(r) && ~(isscalar(r) && (isnan(r) || r == 1))
        dorefine = 1;
    end
end
end

%--------------------------------------------------------------------------
function nnew = edrefinednodes(nold,edgelengths,refinespec)
% The refined number of nodes on each edge, for the propagation.
%   refinespec  a scalar factor (for all edges), or rows [maxlength factor]:
%               the edges not longer than maxlength get the factor (the
%               smallest maxlength that fits); longer edges are not refined.
% The new number is ceil(factor*nold).

if isscalar(refinespec)
    factors = refinespec*ones(size(nold));
else
    [maxlengths,iv] = sort(refinespec(:,1));
    tablefactors = refinespec(iv,2);
    factors = ones(size(nold));
    for ee = 1:numel(nold)
        irow = find(edgelengths(ee) <= maxlengths*(1+1e-9),1);
        if ~isempty(irow)
            factors(ee) = tablefactors(irow);
        end
    end
end
% (The small subtraction avoids that e.g. 1.1*10 = 11.000000000000002
% gives 12.)
nnew = max(ceil(factors.*nold - 1e-9),nold);
end

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
