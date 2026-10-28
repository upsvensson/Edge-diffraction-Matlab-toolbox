function refined = EDrefineQvec(inteqdata,refineopts)
% EDrefineQvec - Refines a solution of the edge source integral equation to more nodes.
% Takes a computed solution q of the integral equation (the vectors
% Qfinal, for all frequencies, from EDintegralequation_convex_tf) and
% gives the values on a finer (or coarser) Gauss-Legendre discretization
% of some or all edges, without solving the integral equation again.
% Optionally, the refined solution is propagated to the receivers.
% PROTOTYPE (19 Oct. 2026): for tests and convergence studies.
%
% Two methods:
%   'nystrom'  (default) The Nystrom interpolation formula,
%                q(z3 <- z2) = q0(z3 <- z2) + sum_z1 w1*H(z3 <- z2 <- z1)*q(z2 <- z1)
%              which holds at any point z3, not only at the nodes. The
%              sum uses the original nodes z1 and the original solution.
%              Refining the "to" variable (z3) takes one step; refining
%              the "from" variable (z2) takes a second step, which uses
%              the values from the first step on the right-hand side:
%                step 1: q(f' <- g) for all q segments with a refined
%                        to-edge f (f' = new nodes, g = old nodes)
%                step 2: q(t' <- f') for all q segments with a refined
%                        from-edge f, using step 1 on the right-hand side.
%              The result has the accuracy of the original solution
%              (the errors of the H-matrix quadrature over z1 remain),
%              but it contains no polynomial interpolation error.
%   'lagrange' Polynomial (Lagrange) interpolation of each q segment, in
%              both variables, from the old to the new Gauss-Legendre
%              nodes. Fast, and needs no kernels, but inaccurate where q
%              is not smooth (near the edge ends).
%
% Applications:
% - Convergence studies that separate the errors of the integral equation
%   and of the propagation: with a solution q_N (N nodes per edge),
%     p(N,N)    the normal result
%     p(N->M,M) q_N refined to M nodes, propagated on the M nodes
%     p(M,M)    the solution with M nodes
%   then p(N,N) - p(N->M,M) is (mostly) the error of the propagation
%   quadrature, and p(N->M,M) - p(M,M) the error of the integral equation.
% - More nodes on one or a few edges, e.g. for a receiver very close to
%   an edge, where the propagation needs a denser discretization.
% - A starting guess, or the prolongation step of a two-grid solver.
%
% Input parameters:
%   inteqdata       EITHER a struct with the fields
%                     .Hsubmatrixdata, .edgedata, .edgetoedgedata, .Sdata,
%                     .envdata, .controlparameters, .Qfinalstore
%                     and optionally .Rdata (needed for the propagation)
%                   OR a text string with the output directory and the
%                   filestem, e.g. '/Users/me/results/cube' (without a
%                   suffix). The structs are then loaded from the files
%                   <filestem>_submatrixdata, _eddata, _ed2data, _Sdata,
%                   _Rdata, _qfinal and _settings (.mat).
%   refineopts      (optional) struct with the fields
%     .method       'nystrom' (default), 'lagrange', or 'none' (no
%                   refinement: the solution is propagated on its own
%                   nodes, whatever .ngauss etc. say; e.g. with
%                   .propagationsplit)
%     .edges        the edges to refine (default: all edges)
%     .nedgeelems   the new number of nodes on each of .edges: a scalar
%                   (the same for all) or a vector with one value per edge
%     .ngauss       instead of .nedgeelems: the new number of nodes is
%                   given by the same rule as in EDinteg_submatrixstructure,
%                   ceil(ngauss*edgelength/maxedgelength), at least 2 (a
%                   negative value gives odd numbers). Default: 2 times the
%                   original ngauss, unless controlparameters.nedgeelemsuser
%                   was used (see .nfactor).
%     .nfactor      instead of .nedgeelems and .ngauss: each of .edges gets
%                   ceil(nfactor*n) nodes, where n is its number of nodes in
%                   the solution. This keeps any node distribution, e.g.
%                   from controlparameters.nedgeelemsuser on a graded
%                   geometry. Default: 2 if controlparameters.nedgeelemsuser
%                   was used and neither .nedgeelems nor .ngauss is given
%                   (the ngauss rule would then give the short pieces fewer
%                   nodes than they have).
%     .frequencyindices  the frequencies (indices into
%                   controlparameters.frequencies) to refine. Default: all.
%     .propagate    1 (default, if Rdata is available): propagate the
%                   refined solution to the receivers.
%     .propagateoriginal  1: also propagate the original solution, with
%                   the same code (default 0). Gives p(N,N) above.
%     .minnedgeelems  a minimum number of nodes on the refined edges,
%                   applied after .ngauss or .nedgeelems (default 0). With
%                   the toolbox rule, very short edges get only 2 nodes for
%                   any ngauss, so they are not refined without this.
%     .propagationsplit  the propagated pressure is also given in parts,
%                   by edge pair (in the order of
%                   Hsubmatrixdata.edgepairlist). The parts add up to the
%                   total. One of
%                   'sharecorner'  part 1: the two edges share a corner
%                                  (Hsubmatrixdata.edgepairsharecorner =
%                                  1), part 2: the others
%                   'sharecorner_near'  part 1: shared corner, part 2: the
%                                  shortest distance between the two edges
%                                  is below .neardistance (default 0.05
%                                  times the longest edge), part 3: the
%                                  others
%                   a logical vector: part 1 = true, part 2 = false
%                   a vector of part numbers 1,2,..., one per edge pair
%     .showtext     0 (default) or 1
%     .parallelworkers  the largest number of parallel workers for the
%                   kernel blocks of the Nystrom step and for the
%                   propagation (default 0 = serial; see
%                   filehandlingparameters.parallelworkers)
%     .forcerecompute  1: a test; the edges in .edges are recomputed with
%                   the Nystrom formula even if their number of nodes is
%                   unchanged. The result should equal the original
%                   solution, up to the tolerance of the solver.
%
% Output parameters:
%   refined         struct with the fields
%     .Hsubmatrixdata   the structure of the refined discretization (the
%                       same edge pairs, but new .nedgeelems,
%                       .Qvecstartnums, .Qvecendnums and quadrature
%                       matrices). Only the fields that EDinteg_souterm and
%                       EDcalcpropagatematrix_block use are valid.
%     .Qfinalstore      cell {nfrequencies,nsolutions} of the refined
%                       solution vectors
%     .Qfirsttermstore  ('nystrom' only) the source terms on the refined
%                       discretization
%     .nedgeelems_old, .nedgeelems_new   [nedges,1]
%     .frequencies      the frequencies that were refined
%     .tfinteqdiff      the refined solution at the receivers, as
%                       tfinteqdiff from EDintegralequation_convex_tf
%                       (if propagated)
%     .tfinteqdiff_original   the original solution at the receivers (if
%                       propagateoriginal = 1)
%     .tfinteqdiff_parts, .tfinteqdiff_original_parts
%                       (if propagationsplit is given) the parts,
%                       [nfrequencies,nreceivers,nsolutions,nparts]
%     .splitlabels      [nedgepairs,1] the part number of each edge pair
%     .splitnames       the names of the parts
%     .splitmask        [nedgepairs,1] the edge pairs of part 1
%     .splitsizes       [2,nparts] the number of edge pairs (row 1) and of
%                       unknowns in the refined solution (row 2) per part
%     .pairdistances    ('sharecorner_near') the shortest distance between
%                       the two edges of each edge pair
%     .nkernelblocks    ('nystrom', difforder > 2) [the number of kernel
%                       blocks, the number that were computed]; blocks
%                       with identical input data are computed once
%     .method, .refinededges, .timing ([refinement propagation] in s)
%                       .method is 'none' if no edge was refined (then the
%                       original solution is returned, e.g. to propagate
%                       it with the split)
%
% The kernel of the H-matrix (EDHkernelvalues), the source term
% (EDinteg_soutermsources) and the propagation (EDpropagateQvec) are the
% same functions that EDintegralequation_convex_tf uses. The propagation
% uses controlparameters.ESIEpropagationimplementation ('accurate' or
% 'legacy'); the split into parts (.propagationsplit) needs 'accurate'.
% EDintegralequation_convex_tf calls this function (with .propagate = 0)
% when controlparameters.ESIEpropagationrefine is used.
%
% Limitations:
% - Only Gauss-Legendre nodes (a higher number of nodes per edge), not
%   arbitrary or graded node distributions.
% - Symmetric submatrices are not reused.
% - For inteqsolutionmethod = 'iterative_fixeddifforder', Qfinal is a
%   truncated series, and the Nystrom step adds one more diffraction
%   order.
%
% See also EDintegralequation_convex_tf, EDpropagateQvec, EDHkernelvalues,
% EDinteg_soutermsources, EDinteg_submatrixstructure
%
% Peter Svensson (peter.svensson@ntnu.no) 28 Oct. 2026
%
% refined = EDrefineQvec(inteqdata,refineopts);

% 19 Oct. 2026 First version, a prototype (by Claude, prompted by PS).
% 20 Oct. 2026 New option .propagationsplit: the propagated pressure in two
% parts, e.g. from the edge pairs that share a corner and from the others.
% If no edge is refined, the original solution is returned directly.
% New options .minnedgeelems and propagationsplit = 'sharecorner_near'
% (three parts) or part numbers.
% 20 Oct. 2026 New method 'none': no refinement, whatever ngauss etc. say.
% 25 Oct. 2026 New option .nfactor: the number of nodes on each refined edge
% is scaled, which keeps the distribution of controlparameters.nedgeelemsuser.
% It is the default (nfactor = 2) when nedgeelemsuser was used. An error if
% an edge would get fewer nodes than before with .ngauss.
% 28 Oct. 2026 The local copies of the H-matrix kernel, the source term and
% the propagation were replaced by the functions that
% EDintegralequation_convex_tf uses: EDHkernelvalues, EDinteg_soutermsources
% and EDpropagateQvec. The results are the same as before, except that the
% propagation now follows controlparameters.ESIEpropagationimplementation
% (before: always 'accurate', the default). New option .parallelworkers.
% 28 Oct. 2026 Kernel blocks with identical input data (e.g. symmetric edge
% triplets) are computed only once, optionally on parallel workers
% (.parallelworkers), and shared. New output field .nkernelblocks =
% [number of blocks, number computed]. The results are unchanged.

t00 = clock;

if nargin < 2
    refineopts = struct;
end

%---------------------------------------------------------------
% The input data

if ischar(inteqdata) || isstring(inteqdata)
    inteqdata = edloadinteqdata(char(inteqdata));
end
Hsub = inteqdata.Hsubmatrixdata;
edgedata = inteqdata.edgedata;
edgetoedgedata = inteqdata.edgetoedgedata;
Sdata = inteqdata.Sdata;
envdata = inteqdata.envdata;
controlparameters = inteqdata.controlparameters;
Qfinalstore = inteqdata.Qfinalstore;
haveRdata = isfield(inteqdata,'Rdata') && ~isempty(inteqdata.Rdata);

method = edgetopt(refineopts,'method','nystrom');
if ~any(strcmpi(method,{'nystrom','lagrange','none'}))
    error('EDrefineQvec:method','refineopts.method must be ''nystrom'', ''lagrange'' or ''none''')
end
method = lower(method);
showtext = edgetopt(refineopts,'showtext',0);
nfreqall = numel(controlparameters.frequencies);
freqindices = edgetopt(refineopts,'frequencyindices',1:nfreqall);
dopropagate = edgetopt(refineopts,'propagate',double(haveRdata));
dopropagateoriginal = edgetopt(refineopts,'propagateoriginal',0);
if (dopropagate || dopropagateoriginal) && ~haveRdata
    error('EDrefineQvec:noRdata','Rdata is needed for the propagation')
end
if size(Qfinalstore,1) ~= nfreqall
    error('EDrefineQvec:Qfinalstore','Qfinalstore must have one row per frequency')
end
nsolutions = size(Qfinalstore,2);

%---------------------------------------------------------------
% The new number of nodes per edge

nedges = numel(Hsub.nedgeelems);
nold = double(Hsub.nedgeelems(:));
edgelist = edgetopt(refineopts,'edges',1:nedges);
edgelist = edgelist(:);
nnew = nold;
haveuserdistribution = isfield(controlparameters,'nedgeelemsuser') && ~isempty(controlparameters.nedgeelemsuser);
givennedgeelems = isfield(refineopts,'nedgeelems') && ~isempty(refineopts.nedgeelems);
givenngauss = isfield(refineopts,'ngauss') && ~isempty(refineopts.ngauss);
nfactor = edgetopt(refineopts,'nfactor',[]);
if isempty(nfactor) && haveuserdistribution && ~givennedgeelems && ~givenngauss
    nfactor = 2;
end
if ~isempty(nfactor) && (givennedgeelems || givenngauss)
    error('EDrefineQvec:nfactor','Give only one of refineopts.nfactor, .nedgeelems and .ngauss')
end
if ~isempty(nfactor)
    if ~isscalar(nfactor) || nfactor <= 0
        error('EDrefineQvec:nfactor','refineopts.nfactor must be a positive scalar')
    end
    nnew(edgelist) = max(ceil(nfactor*nold(edgelist)),2);
elseif givennedgeelems
    nspec = refineopts.nedgeelems(:);
    if numel(nspec) == 1
        nspec = nspec*ones(numel(edgelist),1);
    elseif numel(nspec) ~= numel(edgelist)
        error('EDrefineQvec:nedgeelems','refineopts.nedgeelems must be a scalar or have one value per edge in refineopts.edges')
    end
    nnew(edgelist) = nspec;
else
    newngauss = edgetopt(refineopts,'ngauss',2*controlparameters.ngauss);
    lenratio = edgedata.edgelengthvec(:)/max(edgedata.edgelengthvec);
    if newngauss > 0
        nrule = ceil(newngauss*lenratio);
    else
        nrule = ceil(abs(newngauss)*lenratio/2)*2+1;
    end
    nrule(nrule < 2) = 2;
    if any(nrule(edgelist) < nold(edgelist))
        if haveuserdistribution
            error('EDrefineQvec:ngauss',['With controlparameters.nedgeelemsuser, the ngauss rule would give ',...
                'some edges fewer nodes than in the solution. Use refineopts.nfactor instead.'])
        else
            warning('EDrefineQvec:fewernodes',['The ngauss rule gives some edges fewer nodes than in the ',...
                'solution (e.g. if the solution used nedgeelemsuser). Consider refineopts.nfactor.'])
        end
    end
    nnew(edgelist) = nrule(edgelist);
end
% A minimum number of nodes on the refined edges (e.g. for very short
% edges, which get 2 nodes from the rule above)
minn = edgetopt(refineopts,'minnedgeelems',0);
nnew(edgelist) = max(nnew(edgelist),minn);
if any(nnew < 2) || any(nnew ~= round(nnew))
    error('EDrefineQvec:nedgeelems','The number of nodes per edge must be an integer >= 2')
end
if strcmp(method,'none')
    % No refinement: the discretization of the solution itself (e.g. to
    % propagate it with propagationsplit)
    nnew = nold;
end
isrefined = (nnew ~= nold);
% Test option: recompute the listed edges with the Nystrom formula even if
% their number of nodes is unchanged (the result should then equal the
% original solution, up to the solver tolerance).
if edgetopt(refineopts,'forcerecompute',0)
    isrefined(edgelist) = true;
end

HsubN = edmakelayout(Hsub,nnew,1:size(Hsub.edgepairlist,1));
npairs = size(Hsub.edgepairlist,1);
toedge = Hsub.edgepairlist(:,1);
fromedge = Hsub.edgepairlist(:,2);

% Nodes and weights for the source term (EDinteg_souterm uses
% EDdistelements via gaussvectors, as in EDintegralequation_convex_tf)
maxn = max([nnew;nold]) + 2;
gaussvectors = cell(maxn,1);
for ii = 2:maxn
    [gaussvec1,gaussvec2] = EDdistelements(ii,controlparameters.discretizationtype);
    gaussvectors{ii} = [gaussvec1 gaussvec2];
end

refined = struct;
refined.method = method;
refined.refinededges = find(isrefined);
refined.nedgeelems_old = nold;
refined.nedgeelems_new = nnew;
refined.frequencies = controlparameters.frequencies(freqindices);
refined.Hsubmatrixdata = HsubN;
nfreq = numel(freqindices);
refined.Qfinalstore = cell(nfreq,nsolutions);
if strcmp(method,'nystrom')
    refined.Qfirsttermstore = cell(nfreq,nsolutions);
end

if showtext >= 1
    disp(['   EDrefineQvec (',method,'): ',int2str(sum(isrefined)),' edges refined, nQ = ',...
        int2str(Hsub.Qvecendnums(end)),' -> ',int2str(HsubN.Qvecendnums(end))])
end

%---------------------------------------------------------------
% LAGRANGE: interpolate each q segment in both variables

if ~any(isrefined)
    % Nothing to refine: the original solution (useful for propagating it
    % with the split below)
    method = 'none';
    refined.method = method;
    refined.Qfinalstore = Qfinalstore(freqindices,:);
end

if strcmp(method,'lagrange')
    Lmats = cell(nedges,1);
    for ee = 1:nedges
        if isrefined(ee)
            xold = full(Hsub.quadraturematrix_pos(1:nold(ee),nold(ee)));
            xnew = full(HsubN.quadraturematrix_pos(1:nnew(ee),nnew(ee)));
            Lmats{ee} = edlagrangematrix(xold,xnew);
        else
            Lmats{ee} = speye(nold(ee));
        end
    end
    for jf = 1:nfreq
        for isol = 1:nsolutions
            Qold = Qfinalstore{freqindices(jf),isol};
            Qnew = zeros(HsubN.Qvecendnums(end),1);
            for ii = 1:npairs
                Qseg = Qold(Hsub.Qvecstartnums(ii):Hsub.Qvecendnums(ii));
                % element ifrom + (ito-1)*nfrom: [nfrom,nto]
                Qseg = reshape(Qseg,nold(fromedge(ii)),nold(toedge(ii)));
                Qseg = Lmats{fromedge(ii)}*Qseg*Lmats{toedge(ii)}.';
                Qnew(HsubN.Qvecstartnums(ii):HsubN.Qvecendnums(ii)) = Qseg(:);
            end
            refined.Qfinalstore{jf,isol} = Qnew;
        end
    end
end

%---------------------------------------------------------------
% NYSTROM

if strcmp(method,'nystrom')
    usekernel = controlparameters.difforder > 2;
    if usekernel && strcmp(controlparameters.inteqsolutionmethod,'iterative_fixeddifforder')
        warning('EDrefineQvec:iterative',['With inteqsolutionmethod = ''iterative_fixeddifforder'', ',...
            'the Nystrom step adds one more diffraction order to the truncated series.'])
    end

    % The submatrices (edge triplets) that end in each q segment:
    % listofsubmatrices(nn,2) = the row segment [edge3 edge2],
    % listofsubmatrices(nn,3) = the column segment [edge2 edge1].
    rowseg = double(Hsub.listofsubmatrices(:,2));
    colseg = double(Hsub.listofsubmatrices(:,3));

    % Step 1: the segments with a refined to-edge: q(f' <- g), g old.
    step1segs = find(isrefined(toedge));
    % Step 2: the segments with a refined from-edge: q(t' <- f').
    step2segs = find(isrefined(fromedge));

    % The geometric (frequency-independent) part of the kernel blocks
    % that are needed, and the distances m for the frequency factor.
    % step 1: rows (f' new, g old), columns (g old, h old)
    % step 2: rows (t' new or old, f' new), columns (f' new, g old)
    % Edge triplets with identical input data (the same edge-to-edge
    % geometry, numbers of nodes, edge lengths etc., e.g. symmetric
    % triplets) have identical kernel blocks, which are computed only once
    % (optionally on parallel workers) and shared.
    kernel1 = cell(npairs,1);
    kernel2 = cell(npairs,1);
    if usekernel
        % All the (step, triplet) combinations that are needed
        nnlist1 = find(ismember(rowseg,step1segs));
        nnlist2 = find(ismember(rowseg,step2segs));
        tripletlist = [ones(numel(nnlist1),1) nnlist1(:); 2*ones(numel(nnlist2),1) nnlist2(:)];
        ntrip = size(tripletlist,1);
        refto = edgetoedgedata.reftoshortlistE;
        keys = zeros(ntrip,12);
        nvecs = zeros(ntrip,3);    % [n1 n2 n3]
        for it = 1:ntrip
            nn = tripletlist(it,2);
            e3 = Hsub.edgetripletlist(nn,1); e2 = Hsub.edgetripletlist(nn,2); e1 = Hsub.edgetripletlist(nn,3);
            if tripletlist(it,1) == 1
                nvecs(it,:) = [nold(e1) nold(e2) nnew(e3)];
            else
                nvecs(it,:) = [nold(e1) nnew(e2) nnew(e3)];
            end
            keys(it,:) = [refto(e1,e2) refto(e3,e2) nvecs(it,:) edgedata.edgelengthvec([e1 e2]).' ...
                edgedata.closwedangvec(e2) any(edgedata.offedges(:) == e1) any(edgedata.offedges(:) == e2) ...
                any(edgedata.offedges(:) == e3) double(Hsub.isthinplanetriplet(nn))+1];
        end
        [~,iunique,ikey] = unique(keys,'rows');
        Bstore = cell(numel(iunique),1);
        mstore = cell(numel(iunique),1);
        nworkers = edgetopt(refineopts,'parallelworkers',0);
        parfor (iu = 1:numel(iunique), nworkers)
            it = iunique(iu);
            nn = tripletlist(it,2);
            e1 = Hsub.edgetripletlist(nn,3);
            n1 = nvecs(it,1); n2 = nvecs(it,2); n3 = nvecs(it,3);
            n2vec = [];
            if tripletlist(it,1) == 1
                n2vec = edglnodes(Hsub,n2);
            else
                n2vec = edglnodes(HsubN,n2);
            end
            [Bstore{iu},mstore{iu}] = edHgeomblock(edgedata,edgetoedgedata,Hsub,nn,...
                edglnodes(HsubN,n3),n2vec,edglnodes(Hsub,n1),...
                edglweights(Hsub,n1)*edgedata.edgelengthvec(e1));
        end
        % The blocks of each segment (the arrays are shared, not copied)
        for it = 1:ntrip
            nn = tripletlist(it,2);
            K = struct('B',Bstore{ikey(it)},'mdist',mstore{ikey(it)},'colseg',colseg(nn),'n2',nvecs(it,2));
            ii = rowseg(nn);
            if tripletlist(it,1) == 1
                kernel1{ii}{end+1} = K;
            else
                kernel2{ii}{end+1} = K;
            end
        end
        clear Bstore mstore
        refined.nkernelblocks = [ntrip numel(iunique)];
    end

    % Layouts for the source term of step 1, for the segments [f g] where
    % both f and g are refined: the old number of nodes on g. One layout
    % per refined to-edge f.
    mixlayouts = {};
    mixedges = [];
    for ee = find(isrefined).'
        segs = step1segs(toedge(step1segs) == ee & isrefined(fromedge(step1segs)));
        if ~isempty(segs)
            nmix = nold;
            nmix(ee) = nnew(ee);
            mixlayouts{end+1} = struct('Hsub',edmakelayout(Hsub,nmix,segs),'segs',segs); %#ok<AGROW>
            mixedges(end+1) = ee; %#ok<AGROW>
        end
    end

    for jf = 1:nfreq
        ifreq = freqindices(jf);
        frequency = controlparameters.frequencies(ifreq);
        k = 2*pi*frequency/envdata.cair;
        for isol = 1:nsolutions
            Qold = Qfinalstore{ifreq,isol};

            % The source term on the refined discretization. For the
            % segments with an unrefined from-edge, it is also the source
            % term of step 1.
            Q0N = edsourceterm(envdata,edgedata,edgetoedgedata,HsubN,controlparameters,...
                Sdata,gaussvectors,ifreq,isol,nsolutions);
            refined.Qfirsttermstore{jf,isol} = Q0N;

            % Step 1
            Qstep1 = cell(npairs,1);
            for ii = step1segs(:).'
                if isrefined(fromedge(ii))
                    imix = find(mixedges == toedge(ii));
                    L = mixlayouts{imix};
                    if ~isfield(L,'Q0') || ~isequal(L.Q0key,[jf isol])
                        L.Q0 = edsourceterm(envdata,edgedata,edgetoedgedata,L.Hsub,controlparameters,...
                            Sdata,gaussvectors,ifreq,isol,nsolutions);
                        L.Q0key = [jf isol];
                        mixlayouts{imix} = L;
                    end
                    iloc = find(L.segs == ii);
                    q0seg = L.Q0(L.Hsub.Qvecstartnums(iloc):L.Hsub.Qvecendnums(iloc));
                else
                    q0seg = Q0N(HsubN.Qvecstartnums(ii):HsubN.Qvecendnums(ii));
                end
                qseg = q0seg;
                if usekernel
                    for jj = 1:numel(kernel1{ii})
                        K = kernel1{ii}{jj};
                        rhs = Qold(Hsub.Qvecstartnums(K.colseg):Hsub.Qvecendnums(K.colseg));
                        qseg = qseg + edapplykernel(K,rhs,k);
                    end
                end
                Qstep1{ii} = qseg;
            end

            % Step 2, and the assembly of the refined vector
            Qnew = zeros(HsubN.Qvecendnums(end),1);
            for ii = 1:npairs
                ivnew = HsubN.Qvecstartnums(ii):HsubN.Qvecendnums(ii);
                if isrefined(fromedge(ii))
                    qseg = Q0N(ivnew);
                    if usekernel
                        for jj = 1:numel(kernel2{ii})
                            K = kernel2{ii}{jj};
                            % the right-hand side q(f' <- g): from step 1
                            rhs = Qstep1{K.colseg};
                            qseg = qseg + edapplykernel(K,rhs,k);
                        end
                    end
                    Qnew(ivnew) = qseg;
                elseif isrefined(toedge(ii))
                    Qnew(ivnew) = Qstep1{ii};
                else
                    Qnew(ivnew) = Qold(Hsub.Qvecstartnums(ii):Hsub.Qvecendnums(ii));
                end
            end
            refined.Qfinalstore{jf,isol} = Qnew;
        end
    end
end

trefine = etime(clock,t00);
t01 = clock;

%---------------------------------------------------------------
% Propagation to the receivers

% The optional split of the propagated pressure into two parts, by edge
% pair (q segment): part 1 = the segments where splitmask is true, part 2
% = the others. The parts add up to the total.
splitopt = edgetopt(refineopts,'propagationsplit',[]);
splitlabels = [];
if ~isempty(splitopt)
    if ischar(splitopt) || isstring(splitopt)
        sharecorner = full(Hsub.edgepairsharecorner(:)) == 1;
        switch lower(char(splitopt))
            case 'sharecorner'
                splitlabels = 2 - double(sharecorner);
                refined.splitnames = {'shared corner','other'};
            case 'sharecorner_near'
                neardistance = edgetopt(refineopts,'neardistance',0.05*max(edgedata.edgelengthvec));
                pairdist = zeros(npairs,1);
                for ii = 1:npairs
                    ea = Hsub.edgepairlist(ii,1);
                    eb = Hsub.edgepairlist(ii,2);
                    pairdist(ii) = edsegmentdistance(edgedata.edgestartcoords(ea,:),edgedata.edgeendcoords(ea,:),...
                        edgedata.edgestartcoords(eb,:),edgedata.edgeendcoords(eb,:));
                end
                splitlabels = 3*ones(npairs,1);
                splitlabels(pairdist < neardistance) = 2;
                splitlabels(sharecorner) = 1;
                refined.splitnames = {'shared corner',sprintf('near (< %g m)',neardistance),'other'};
                refined.pairdistances = pairdist;
            otherwise
                error('EDrefineQvec:split','Unknown refineopts.propagationsplit: %s',char(splitopt))
        end
    elseif islogical(splitopt)
        splitlabels = 2 - double(splitopt(:));
        refined.splitnames = {'true','false'};
    else
        splitlabels = double(splitopt(:));
        refined.splitnames = arrayfun(@(x) sprintf('part %d',x),1:max(splitlabels),'UniformOutput',false);
    end
    if numel(splitlabels) ~= npairs || any(splitlabels < 1) || any(splitlabels ~= round(splitlabels))
        error('EDrefineQvec:split','refineopts.propagationsplit must give one part number (1,2,...) per edge pair')
    end
    nparts = max(splitlabels);
    refined.splitlabels = splitlabels;
    refined.splitmask = (splitlabels == 1);
    refined.splitsizes = zeros(2,nparts);
    for jp = 1:nparts
        refined.splitsizes(:,jp) = [sum(splitlabels == jp); edcountunknowns(HsubN,splitlabels == jp)];
    end
end

if dopropagate
    [refined.tfinteqdiff,parts] = edpropagate(edgedata,edgetoedgedata,HsubN,inteqdata.Rdata,...
        refined.Qfinalstore,controlparameters,envdata,Sdata,freqindices,showtext,splitlabels,...
        edgetopt(refineopts,'parallelworkers',0));
    if ~isempty(splitlabels)
        refined.tfinteqdiff_parts = parts;
    end
end
if dopropagateoriginal
    [refined.tfinteqdiff_original,parts] = edpropagate(edgedata,edgetoedgedata,Hsub,inteqdata.Rdata,...
        Qfinalstore(freqindices,:),controlparameters,envdata,Sdata,freqindices,showtext,splitlabels,...
        edgetopt(refineopts,'parallelworkers',0));
    if ~isempty(splitlabels)
        refined.tfinteqdiff_original_parts = parts;
    end
end
refined.timing = [trefine etime(clock,t01)];

if showtext >= 1
    disp(['   EDrefineQvec: refinement ',num2str(refined.timing(1),3),' s, propagation ',num2str(refined.timing(2),3),' s'])
end

%==========================================================================
% LOCAL FUNCTIONS
%==========================================================================

function val = edgetopt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name))
    val = opts.(name);
else
    val = default;
end

%--------------------------------------------------------------------------
function D = edloadinteqdata(stem)
% Loads the structs from the stage files <stem>_*.mat

D = struct;
f = load([stem,'_submatrixdata.mat'],'Hsubmatrixdata');   D.Hsubmatrixdata = f.Hsubmatrixdata;
f = load([stem,'_ed2data.mat'],'edgedata','edgetoedgedata');
D.edgedata = f.edgedata;
D.edgetoedgedata = f.edgetoedgedata;
f = load([stem,'_Sdata.mat'],'Sdata');                    D.Sdata = f.Sdata;
f = load([stem,'_qfinal.mat'],'Qfinalstore');             D.Qfinalstore = f.Qfinalstore;
f = load([stem,'_settings.mat'],'envdata','controlparameters');
D.envdata = f.envdata;
D.controlparameters = f.controlparameters;
if exist([stem,'_Rdata.mat'],'file') == 2
    f = load([stem,'_Rdata.mat'],'Rdata');                D.Rdata = f.Rdata;
end

%--------------------------------------------------------------------------
function x = edglnodes(Hsub,n)
x = full(Hsub.quadraturematrix_pos(1:n,n));

function w = edglweights(Hsub,n)
w = full(Hsub.quadraturematrix_weights(1:n,n));

%--------------------------------------------------------------------------
function HsubX = edmakelayout(Hsub,nedgeelems,pairlist)
% A copy of Hsub with another number of nodes per edge, and only the edge
% pairs in pairlist (in that order). The Qvec layout and the quadrature
% matrices are recomputed; the submatrix fields refer to the full list
% and are only valid if pairlist is the full list.

HsubX = Hsub;
HsubX.nedgeelems = cast(nedgeelems(:),class(Hsub.nedgeelems));
HsubX.edgepairlist = Hsub.edgepairlist(pairlist,:);
if isfield(Hsub,'isthinplaneedgepair') && numel(Hsub.isthinplaneedgepair) == size(Hsub.edgepairlist,1)
    HsubX.isthinplaneedgepair = Hsub.isthinplaneedgepair(pairlist);
end
nelem = double(nedgeelems(HsubX.edgepairlist(:,1))).*double(nedgeelems(HsubX.edgepairlist(:,2)));
HsubX.Qvecendnums = cumsum(nelem(:));
HsubX.Qvecstartnums = [1;HsubX.Qvecendnums(1:end-1)+1];
nmax = max([double(nedgeelems(:));size(Hsub.quadraturematrix_pos,1)]);
qpos = sparse(nmax,nmax);
qw = sparse(nmax,nmax);
[nprev1,nprev2] = size(Hsub.quadraturematrix_pos);
qpos(1:nprev1,1:nprev2) = Hsub.quadraturematrix_pos;
qw(1:nprev1,1:nprev2) = Hsub.quadraturematrix_weights;
for n = unique(double(nedgeelems(:))).'
    if n > nprev2 || nnz(qpos(:,n)) == 0
        [xn,wn] = EDgetGLrule(n);
        qpos(1:n,n) = xn(:);
        qw(1:n,n) = wn(:);
    end
end
HsubX.quadraturematrix_pos = qpos;
HsubX.quadraturematrix_weights = qw;

%--------------------------------------------------------------------------
function L = edlagrangematrix(xold,xnew)
% Lagrange interpolation matrix: values at xnew = L*values at xold
% (barycentric formula)

xold = xold(:);
xnew = xnew(:);
n = numel(xold);
wbary = ones(n,1);
for jj = 1:n
    wbary(jj) = 1/prod(xold(jj) - xold([1:jj-1 jj+1:n]));
end
D = xnew - xold.';                 % [nnew,nold]
exact = (abs(D) < 1e-14);
D(exact) = 1;
T = (wbary.')./D;
L = T./sum(T,2);
[ir,jc] = find(exact);
for kk = 1:numel(ir)
    L(ir(kk),:) = 0;
    L(ir(kk),jc(kk)) = 1;
end

%--------------------------------------------------------------------------
function Q0 = edsourceterm(envdata,edgedata,edgetoedgedata,HsubX,controlparameters,...
    Sdata,gaussvectors,ifreq,isol,nsolutions)
% The source term on the layout HsubX, as in EDintegralequation_convex_tf
% (with EDinteg_soutermsources): the sum over the sources if they are
% added (nsolutions = 1), otherwise source number isol.

frequency = controlparameters.frequencies(ifreq);
nsources = size(Sdata.coordinates,1);
if nsolutions == 1
    Q0 = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,HsubX,controlparameters,...
        Sdata,frequency,ifreq,gaussvectors,1:nsources,1,0);
else
    Q0 = EDinteg_soutermsources(envdata,edgedata,edgetoedgedata,HsubX,controlparameters,...
        Sdata,frequency,ifreq,gaussvectors,isol,0,0);
end
if isscalar(Q0)
    Q0 = zeros(HsubX.Qvecendnums(end),1);
end

%--------------------------------------------------------------------------
function qadd = edapplykernel(K,rhs,k)
% One kernel block times the right-hand side segment:
%   qadd(row) = sum_i1 B(row,i1)*exp(-1i*k*m(i2,i1))*rhs(i1 + (i2-1)*n1)
% where row = i2 + (i3-1)*n2.

n2 = K.n2;
n1 = size(K.B,2);
n3 = size(K.B,1)/n2;
Rmat = reshape(rhs,n1,n2).*exp(-1i*k*K.mdist.');     % [n1,n2]
qadd = zeros(n2,n3);
Bthree = reshape(K.B,n2,n3,n1);
for i1 = 1:n1
    qadd = qadd + Bthree(:,:,i1).*Rmat(i1,:).';
end
qadd = qadd(:);

%--------------------------------------------------------------------------
function [B,mdist] = edHgeomblock(edgedata,edgetoedgedata,Hsub,nn,n3vec,n2vec,n1vec,dzvec)
% The frequency-independent part of one H-submatrix (edge triplet nn =
% [edge3 edge2 edge1]), for the given normalized node positions on the
% three edges and the weights dzvec (in m) on edge 1, computed by
% EDHkernelvalues (the same kernel as in EDintegralequation_convex_tf) for
% the full grid of node combinations.
%   B      [n3*n2,n1]; row i2 + (i3-1)*n2, column i1. The element of H is
%          B*exp(-1i*k*mdist(i2,i1)).
%   mdist  [n2,n1] the distance from node i1 on edge 1 to node i2 on edge 2

edge3 = Hsub.edgetripletlist(nn,1);
edge2 = Hsub.edgetripletlist(nn,2);
edge1 = Hsub.edgetripletlist(nn,3);
n1 = numel(n1vec); n2 = numel(n2vec); n3 = numel(n3vec);
thinplaneboost = double(Hsub.isthinplanetriplet(nn)) + 1;

% Index lists for the elements of B, column by column: element e
% (0-based) has i2 = mod(e,n2), i3 = mod(floor(e/n2),n3),
% i1 = floor(e/(n2*n3)) (plus 1)
e = (0:n1*n2*n3-1).';
n2ind = mod(e,n2) + 1;
n3ind = mod(floor(e/n2),n3) + 1;
n1ind = floor(e/(n2*n3)) + 1;

[Bvalues,mvalues] = EDHkernelvalues(edgedata,edgetoedgedata,edge1,edge2,edge3,...
    n1vec,n2vec,n3vec,dzvec(:),n1ind,n2ind,n3ind,thinplaneboost);
B = reshape(Bvalues,n3*n2,n1);
mdist = reshape(mvalues,n3*n2,n1);
mdist = mdist(1:n2,:);

%--------------------------------------------------------------------------
function [tfinteqdiff,parts] = edpropagate(edgedata,edgetoedgedata,HsubX,Rdata,Qstore,...
    controlparameters,envdata,Sdata,freqindices,showtext,splitlabels,nworkers)
% Propagation of the solutions in Qstore to all receivers, with
% EDpropagateQvec (as in EDintegralequation_convex_tf). If splitlabels (a
% part number 1..nparts per edge pair) is given, parts gets the parts
% [nfreq,nreceivers,nsolutions,nparts], which add up to the total.

implementation = 'accurate';
if isfield(controlparameters,'ESIEpropagationimplementation') && ischar(controlparameters.ESIEpropagationimplementation)
    implementation = controlparameters.ESIEpropagationimplementation;
end
propsettings = struct('Rstart',controlparameters.Rstart,'implementation',implementation,...
    'doallSRcombinations',Sdata.doallSRcombinations,'savediff2',0,'nworkers',nworkers,...
    'showtext',showtext);
propsettings.splitlabels = splitlabels;
[tfinteqdiff,~,propinfo] = EDpropagateQvec(edgedata,edgetoedgedata,HsubX,Rdata,Qstore,{},...
    controlparameters.frequencies(freqindices),envdata,propsettings);
parts = [];
if isfield(propinfo,'parts')
    parts = propinfo.parts;
end

%--------------------------------------------------------------------------
function n = edcountunknowns(HsubX,mask)
% The number of unknowns in the edge pairs where mask is true
n = sum(double(HsubX.Qvecendnums(mask)) - double(HsubX.Qvecstartnums(mask)) + 1);

%--------------------------------------------------------------------------
function d = edsegmentdistance(p1,q1,p2,q2)
% The shortest distance between the line segments p1-q1 and p2-q2 in 3D
% (the closest points of the two segments; Ericson, Real-Time Collision
% Detection, 2005, section 5.1.9).

d1 = q1 - p1;
d2 = q2 - p2;
r = p1 - p2;
a = dot(d1,d1);
e = dot(d2,d2);
f = dot(d2,r);
if a <= eps && e <= eps
    d = norm(p1 - p2);
    return
end
if a <= eps
    s = 0;
    t = min(max(f/e,0),1);
else
    c = dot(d1,r);
    if e <= eps
        t = 0;
        s = min(max(-c/a,0),1);
    else
        b = dot(d1,d2);
        denom = a*e - b*b;
        if denom > eps*a*e
            s = min(max((b*f - c*e)/denom,0),1);
        else
            s = 0;
        end
        t = (b*s + f)/e;
        if t < 0
            t = 0;
            s = min(max(-c/a,0),1);
        elseif t > 1
            t = 1;
            s = min(max((b - c)/a,0),1);
        end
    end
end
d = norm((p1 + d1*s) - (p2 + d2*t));
