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
%     .method       'nystrom' (default) or 'lagrange'
%     .edges        the edges to refine (default: all edges)
%     .nedgeelems   the new number of nodes on each of .edges: a scalar
%                   (the same for all) or a vector with one value per edge
%     .ngauss       instead of .nedgeelems: the new number of nodes is
%                   given by the same rule as in EDinteg_submatrixstructure,
%                   ceil(ngauss*edgelength/maxedgelength), at least 2 (a
%                   negative value gives odd numbers). Default: 2 times the
%                   original ngauss.
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
%     .method, .refinededges, .timing ([refinement propagation] in s)
%                       .method is 'none' if no edge was refined (then the
%                       original solution is returned, e.g. to propagate
%                       it with the split)
%
% Limitations of this prototype:
% - Only Gauss-Legendre nodes (a higher number of nodes per edge), not
%   arbitrary or graded node distributions.
% - The H-matrix kernel is a copy of the code in
%   EDintegralequation_convex_tf (local function edHgeomblock below);
%   symmetric submatrices are not reused.
% - For inteqsolutionmethod = 'iterative_fixeddifforder', Qfinal is a
%   truncated series, and the Nystrom step adds one more diffraction
%   order.
% - The propagation always uses the 'accurate' implementation
%   (EDcalcpropagatematrix_block), serially.
%
% See also EDintegralequation_convex_tf, EDinteg_souterm,
% EDcalcpropagatematrix_block, EDinteg_submatrixstructure
%
% Peter Svensson (peter.svensson@ntnu.no) 20 Oct. 2026
%
% refined = EDrefineQvec(inteqdata,refineopts);

% 19 Oct. 2026 First version, a prototype (by Claude, prompted by PS).
% 20 Oct. 2026 New option .propagationsplit: the propagated pressure in two
% parts, e.g. from the edge pairs that share a corner and from the others.
% If no edge is refined, the original solution is returned directly.
% New options .minnedgeelems and propagationsplit = 'sharecorner_near'
% (three parts) or part numbers.

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
if ~any(strcmpi(method,{'nystrom','lagrange'}))
    error('EDrefineQvec:method','refineopts.method must be ''nystrom'' or ''lagrange''')
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
if isfield(refineopts,'nedgeelems') && ~isempty(refineopts.nedgeelems)
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
    nnew(edgelist) = nrule(edgelist);
end
% A minimum number of nodes on the refined edges (e.g. for very short
% edges, which get 2 nodes from the rule above)
minn = edgetopt(refineopts,'minnedgeelems',0);
nnew(edgelist) = max(nnew(edgelist),minn);
if any(nnew < 2) || any(nnew ~= round(nnew))
    error('EDrefineQvec:nedgeelems','The number of nodes per edge must be an integer >= 2')
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
    kernel1 = cell(npairs,1);
    kernel2 = cell(npairs,1);
    if usekernel
        for ii = step1segs(:).'
            nnlist = find(rowseg == ii);
            kernel1{ii} = cell(numel(nnlist),1);
            for jj = 1:numel(nnlist)
                nn = nnlist(jj);
                e3 = Hsub.edgetripletlist(nn,1); e2 = Hsub.edgetripletlist(nn,2); e1 = Hsub.edgetripletlist(nn,3);
                [B,mdist] = edHgeomblock(edgedata,edgetoedgedata,Hsub,nn,...
                    edglnodes(HsubN,nnew(e3)),edglnodes(Hsub,nold(e2)),edglnodes(Hsub,nold(e1)),...
                    edglweights(Hsub,nold(e1))*edgedata.edgelengthvec(e1));
                kernel1{ii}{jj} = struct('B',B,'mdist',mdist,'colseg',colseg(nn),'n2',nold(e2));
            end
        end
        for ii = step2segs(:).'
            nnlist = find(rowseg == ii);
            kernel2{ii} = cell(numel(nnlist),1);
            for jj = 1:numel(nnlist)
                nn = nnlist(jj);
                e3 = Hsub.edgetripletlist(nn,1); e2 = Hsub.edgetripletlist(nn,2); e1 = Hsub.edgetripletlist(nn,3);
                [B,mdist] = edHgeomblock(edgedata,edgetoedgedata,Hsub,nn,...
                    edglnodes(HsubN,nnew(e3)),edglnodes(HsubN,nnew(e2)),edglnodes(Hsub,nold(e1)),...
                    edglweights(Hsub,nold(e1))*edgedata.edgelengthvec(e1));
                kernel2{ii}{jj} = struct('B',B,'mdist',mdist,'colseg',colseg(nn),'n2',nnew(e2));
            end
        end
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
        refined.Qfinalstore,controlparameters,envdata,Sdata,freqindices,showtext,splitlabels);
    if ~isempty(splitlabels)
        refined.tfinteqdiff_parts = parts;
    end
end
if dopropagateoriginal
    [refined.tfinteqdiff_original,parts] = edpropagate(edgedata,edgetoedgedata,Hsub,inteqdata.Rdata,...
        Qfinalstore(freqindices,:),controlparameters,envdata,Sdata,freqindices,showtext,splitlabels);
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
% The source term on the layout HsubX, as in EDintegralequation_convex_tf:
% the sum over the sources if they are added (nsolutions = 1), otherwise
% source number isol.

frequency = controlparameters.frequencies(ifreq);
nsources = size(Sdata.coordinates,1);
if nsolutions == 1
    sourcelist = 1:nsources;
else
    sourcelist = isol;
end
Q0 = 0;
for isou = sourcelist
    if strcmp(Sdata.sourcetype,'monopole')
        Q0add = EDinteg_souterm(envdata,edgedata,edgetoedgedata,HsubX,controlparameters,...
            Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
            Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,...
            Sdata.rSsho(Sdata.reftoshortlistS(:,isou)),Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou)),...
            Sdata.zSsho(Sdata.reftoshortlistS(:,isou)),0);
    else   % polygonpiston
        Q0add = 0;
        for np = 1:size(Sdata.pistongausscoordinates{1},1)
            Q0one = EDinteg_souterm(envdata,edgedata,edgetoedgedata,HsubX,controlparameters,...
                Sdata.vispartedgesfroms(:,isou),Sdata.vispartedgesfroms_start(:,isou),...
                Sdata.vispartedgesfroms_end(:,isou),frequency,gaussvectors,...
                Sdata.rSsho(Sdata.reftoshortlistS(:,isou,np)),Sdata.thetaSsho(Sdata.reftoshortlistS(:,isou,np)),...
                Sdata.zSsho(Sdata.reftoshortlistS(:,isou,np)),0);
            Q0add = Q0add + Q0one*Sdata.pistongaussweights{isou}(np,1);
        end
    end
    Q0 = Q0 + Q0add*Sdata.sourceamplitudes(isou,ifreq);
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
% three edges and the weights dzvec (in m) on edge 1.
%   B      [n3*n2,n1]; row i2 + (i3-1)*n2, column i1. The element of H is
%          B*exp(-1i*k*mdist(i2,i1)).
%   mdist  [n2,n1] the distance from node i1 on edge 1 to node i2 on edge 2
% A COPY of the kernel in EDintegralequation_convex_tf (PART 1), for
% arbitrary node sets.

Hsubconstant = -1/4/pi/2;
refto = edgetoedgedata.reftoshortlistE;

edge3 = Hsub.edgetripletlist(nn,1);
edge2 = Hsub.edgetripletlist(nn,2);
edge1 = Hsub.edgetripletlist(nn,3);
len2 = edgedata.edgelengthvec(edge2);
n1 = numel(n1vec); n2 = numel(n2vec); n3 = numel(n3vec);
n1vec = n1vec(:).';
n2vec = n2vec(:);
n3vec = n3vec(:);
dzvec = dzvec(:).';

thinplaneboost = double(Hsub.isthinplanetriplet(nn)) + 1;

ny = pi/(2*pi-edgedata.closwedangvec(edge2));
thetain  = edgetoedgedata.thetae1sho(refto(edge1,edge2));
thetain_other_end_of_edge  = edgetoedgedata.thetae2sho(refto(edge1,edge2));
acrossface_in = double((abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) ...
    && abs(thetain-thetain_other_end_of_edge) < 1e-10);
thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
thetaout_other_end_of_edge = edgetoedgedata.thetae2sho(refto(edge3,edge2));
acrossface_out = double((abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) ...
    && abs(thetaout-thetaout_other_end_of_edge) < 1e-10);

ze2_re2 = n2vec*len2;                                        % [n2,1]
ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2)) );
ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2)) );

if acrossface_in*acrossface_out == 0
    xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
    xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
    yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
    yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
    xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
    ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);
    [re1_re2,thetae1_re2,~] = EDcoordtrans1([xe1_re2.' ye1_re2.' ze1_re2.'],[0 0 0;0 0 len2],[0 1 0]);
    re1_re2 = re1_re2(:).';
    thetae1_re2 = thetae1_re2(:).';

    xedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
    xedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
    yedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
    yedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
    xe3_re2 = xedge3start + n3vec.*(xedge3end-xedge3start);
    ye3_re2 = yedge3start + n3vec.*(yedge3end-yedge3start);
    [re3_re2,thetae3_re2,~] = EDcoordtrans1([xe3_re2 ye3_re2 ze3_re2],[0 0 0;0 0 len2],[0 1 0]);
    re3_re2 = re3_re2(:);
    thetae3_re2 = thetae3_re2(:);
else
    re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2)) );
    re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2)) );
end

% Index vectors for the rows: row = i2 + (i3-1)*n2
i2row = repmat((1:n2).',n3,1);
i3row = reshape(repmat(1:n3,n2,1),n2*n3,1);

mdist = sqrt( (ze1_re2 - ze2_re2).^2 + re1_re2.^2 );          % [n2,n1]
ldist = sqrt( (ze3_re2(i3row) - ze2_re2(i2row)).^2 + re3_re2(i3row).^2 );   % [n3*n2,1]

M = mdist(i2row,:);                                          % [n3*n2,n1]
dz_overm = dzvec./M;
ch = ((ze2_re2(i2row) - ze1_re2).*(ze2_re2(i2row) - ze3_re2(i3row)) + M.*ldist)./re1_re2./re3_re2(i3row);
ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
ch = ( ch + 1./ch)/2;

if acrossface_in == 1 && acrossface_out == 1
    beta = 4*( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))));
else
    t1 = thetae1_re2;            % [1,n1]
    t3 = thetae3_re2(i3row);     % [n3*n2,1]
    beta = ( sin(ny*(pi + t1 + t3))./(ch - cos(ny*(pi + t1 + t3))+eps*10)) + ...
           ( sin(ny*(pi + t1 - t3))./(ch - cos(ny*(pi + t1 - t3))+eps*10)) + ...
           ( sin(ny*(pi - t1 + t3))./(ch - cos(ny*(pi - t1 + t3))+eps*10)) + ...
           ( sin(ny*(pi - t1 - t3))./(ch - cos(ny*(pi - t1 - t3))+eps*10));
end

if ismember(edge1,edgedata.offedges) || ismember(edge2,edgedata.offedges) || ismember(edge3,edgedata.offedges)
    beta = beta*0;
end

B = (2-acrossface_out)*Hsubconstant*ny*beta.*dz_overm*thinplaneboost;

%--------------------------------------------------------------------------
function [tfinteqdiff,parts] = edpropagate(edgedata,edgetoedgedata,HsubX,Rdata,Qstore,...
    controlparameters,envdata,Sdata,freqindices,showtext,splitlabels)
% Propagation of the solutions in Qstore to all receivers, as in PART 4
% of EDintegralequation_convex_tf ('accurate', in blocks, serially).
% If splitlabels (a part number 1..nparts per edge pair) is given, parts
% gets the parts [nfreq,nreceivers,nsolutions,nparts]. The propagation
% matrix is computed once; each column belongs to one edge pair, so the
% parts add up to the total exactly.

nfreq = numel(freqindices);
nsolutions = size(Qstore,2);
nreceivers = size(Rdata.coordinates,1);
nQ = double(HsubX.Qvecendnums(end));
kvec = 2*pi*controlparameters.frequencies(freqindices)/envdata.cair;

npairs = size(HsubX.edgepairlist,1);
flags = false(npairs,1);
for jf = 1:nfreq
    for isol = 1:nsolutions
        for ii = 1:npairs
            flags(ii) = flags(ii) || any(Qstore{jf,isol}(HsubX.Qvecstartnums(ii):HsubX.Qvecendnums(ii)));
        end
    end
end
flags = double(flags);

tfinteqdiff = zeros(nfreq,nreceivers,nsolutions);
dosplit = ~isempty(splitlabels);
parts = [];
if dosplit
    nparts = max(splitlabels);
    parts = zeros(nfreq,nreceivers,nsolutions,nparts);
    % the part number of each element of the Q vector
    partofQ = zeros(nQ,1);
    for ii = 1:npairs
        partofQ(HsubX.Qvecstartnums(ii):HsubX.Qvecendnums(ii)) = splitlabels(ii);
    end
end
blocksize = max(1,min(nreceivers,floor(4e6/nQ)));
for r0 = 1:blocksize:nreceivers
    recs = r0:min(r0+blocksize-1,nreceivers);
    nb = numel(recs);
    if showtext >= 1
        disp(['      Propagation: receivers ',int2str(recs(1)),' to ',int2str(recs(end))])
    end
    [Fgeom,pathlengths] = EDcalcpropagatematrix_block(edgedata,edgetoedgedata,HsubX,...
        Rdata.vispartedgesfromr(:,recs),...
        reshape(Rdata.rRsho(Rdata.reftoshortlistR(:,recs)),[],nb),...
        reshape(Rdata.thetaRsho(Rdata.reftoshortlistR(:,recs)),[],nb),...
        reshape(Rdata.zRsho(Rdata.reftoshortlistR(:,recs)),[],nb),flags);
    ivnonzero = find(any(Fgeom,1));
    Fgeom = Fgeom(:,ivnonzero);
    pathlengths = pathlengths(:,ivnonzero) - controlparameters.Rstart;
    for jf = 1:nfreq
        Fk = Fgeom.*exp(-1i*kvec(jf)*pathlengths);
        if nsolutions == 1 || Sdata.doallSRcombinations == 1
            for isol = 1:nsolutions
                Qone = Qstore{jf,isol}(ivnonzero);
                tfinteqdiff(jf,recs,isol) = Fk*Qone;
                if dosplit
                    pnz = partofQ(ivnonzero);
                    for jp = 1:nparts
                        inp = (pnz == jp);
                        parts(jf,recs,isol,jp) = Fk(:,inp)*Qone(inp);
                    end
                end
            end
        else
            for jr = 1:nb
                isol = recs(jr);
                if isol <= nsolutions
                    Qone = Qstore{jf,isol}(ivnonzero);
                    tfinteqdiff(jf,isol,isol) = Fk(jr,:)*Qone;
                    if dosplit
                        pnz = partofQ(ivnonzero);
                        for jp = 1:nparts
                            inp = (pnz == jp);
                            parts(jf,isol,isol,jp) = Fk(jr,inp)*Qone(inp);
                        end
                    end
                end
            end
        end
    end
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
