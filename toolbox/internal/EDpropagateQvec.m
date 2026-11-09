function [P_receiver,P_nodiff2,propinfo] = EDpropagateQvec(edgedata,edgetoedgedata,Hsubmatrixdata,Rdata,...
    Qfinalstore,Qfirsttermstore,frequencies,envdata,propsettings)
% EDpropagateQvec - Propagates solutions of the edge source integral equation to the receivers.
% The edge source amplitudes q (Qfinal, one vector per frequency and
% solution) are propagated to the receivers,
%   p(R) = int int F(R <- z2 <- z1) q(z2 <- z1) dz1 dz2,
% with the quadrature nodes and weights of the layout Hsubmatrixdata. The
% propagation row vector depends on the frequency only through the factor
% exp(-1i*k*(L-Rstart)), where L = l + m is the path length via the two
% edge nodes to the receiver. Everything else (the beta terms, the
% singularity treatment, 1/(l*m) and the quadrature weights) depends only
% on the geometry, so the geometric part is computed once per receiver and
% only the exponential factor per frequency.
%
% For the 'accurate' and 'exactpole' propagation, the receivers are handled
% in blocks by EDcalcpropagatematrix_block, which processes all the
% receivers of a block at once for each edge pair; 'exactpole' uses its
% exact-pole treatment of the near-singular rows (polemethod = 'exact'). The
% 'legacy' propagation uses EDcalcpropagatematrix for one receiver at a
% time, and so do thin holes (isthinhole = 1), where 'exactpole' is replaced
% by 'accurate' (noted in propinfo.implementationnote). Both loops can run
% on parallel workers (parfor).
%
% Used by EDintegralequation_convex_tf (PART 4) and EDrefineQvec. The
% layout can be that of the solution, or a refined one (EDrefineQvec).
%
% Input parameters:
%   edgedata, edgetoedgedata, Hsubmatrixdata, Rdata
%                   structs, see help EDdatastructs. Hsubmatrixdata gives
%                   the layout of the Q vectors (.edgepairlist,
%                   .nedgeelems, .Qvecstartnums, .Qvecendnums, and the
%                   quadrature matrices).
%   Qfinalstore     cell {nfrequencies,nsolutions} with the Q vectors
%   Qfirsttermstore cell of the same size with the source terms; only used
%                   if propsettings.savediff2 = 1 (otherwise it can be {})
%   frequencies     [1,nfrequencies] the frequencies of the rows of
%                   Qfinalstore, in Hz
%   envdata         struct; the field .cair is used
%   propsettings    struct with the fields
%     .Rstart           the phase reference distance (controlparameters.Rstart)
%     .implementation   'accurate', 'exactpole' or 'legacy'
%                       (controlparameters.ESIEpropagationimplementation)
%     .doallSRcombinations   (Sdata.doallSRcombinations)
%     .savediff2        1: also the result without the second-order
%                       diffraction, propagated Qfinal - Q_firstterm
%                       (filehandlingparameters.savediff2result)
%     .nworkers         the largest number of parallel workers (0 = serial)
%     .showtext         (optional, default 0)
%     .isthinhole       (optional, default 0)
%     .splitlabels      (optional, default []) one part number (1,2,...)
%                       per edge pair (row of Hsubmatrixdata.edgepairlist).
%                       The propagated pressure is then also given in
%                       parts, which add up to the total. Only for the
%                       'accurate' and 'exactpole' propagation.
%
% Output parameters:
%   P_receiver      [nfrequencies,nreceivers,nsolutions] the propagated
%                   pressure. With doallSRcombinations = 0, only source N
%                   to receiver N is computed.
%   P_nodiff2       the same, without the second-order diffraction, if
%                   savediff2 = 1; otherwise [].
%   propinfo        struct with the fields
%     .propdiagnostics    the totals of the singularity corrections (see
%                         EDcalcpropagatematrix_block); for 'legacy', the
%                         values of the last receiver. For 'exactpole' also
%                         .nrows_exactpole and .nrows_fallback.
%     .propagationblocksize  the number of receivers per block ('accurate'),
%                         or [] ('legacy')
%     .ivproblematicmatrix   ('legacy') the near-singular entries per receiver
%     .doesQsegmenthavevalues  [nedgepairs,1] the edge pairs where some Q
%                         vector has non-zero values
%     .parts              (if splitlabels is given) [nfrequencies,
%                         nreceivers,nsolutions,nparts]
%     .implementationnote  '' or a note that 'exactpole' was replaced by
%                         'accurate' (thin holes)
%
% See also EDintegralequation_convex_tf, EDrefineQvec,
% EDcalcpropagatematrix_block, EDcalcpropagatematrix
%
% Peter Svensson (peter.svensson@ntnu.no) 9 Nov. 2026
%
% [P_receiver,P_nodiff2,propinfo] = EDpropagateQvec(edgedata,edgetoedgedata,Hsubmatrixdata,Rdata,...
%    Qfinalstore,Qfirsttermstore,frequencies,envdata,propsettings);

% 28 Oct. 2026 First version (by Claude, prompted by PS): PART 4 of
% EDintegralequation_convex_tf, moved here unchanged, so that
% EDrefineQvec can use the same code. New option: the split into parts by
% edge pair (from EDrefineQvec).
% 9 Nov. 2026 (by Claude, prompted by PS) New implementation 'exactpole':
% the block propagation with polemethod = 'exact'. The totals of the new
% propdiagnostics fields are collected. For thin holes, 'exactpole' is
% replaced by 'accurate'.

showtext = edgetsetting(propsettings,'showtext',0);
isthinhole = edgetsetting(propsettings,'isthinhole',0);
splitlabels = edgetsetting(propsettings,'splitlabels',[]);
nworkers = propsettings.nworkers;
savediff2 = propsettings.savediff2;
doallSRcombinations = propsettings.doallSRcombinations;
Rstartvalue = propsettings.Rstart;
propimplementation = propsettings.implementation;
if ~any(strcmpi(propimplementation,{'accurate','exactpole','legacy'}))
    error('EDpropagateQvec:implementation','propsettings.implementation must be ''accurate'', ''exactpole'' or ''legacy''')
end
implementationnote = '';
if strcmpi(propimplementation,'exactpole') && isthinhole ~= 0
    propimplementation = 'accurate';
    implementationnote = 'ESIEpropagationimplementation = ''exactpole'' is not implemented for thin holes; ''accurate'' was used.';
end
if strcmpi(propimplementation,'exactpole')
    polemethod = 'exact';
else
    polemethod = 'taylor';
end

[nfrequencies,nsolutions] = size(Qfinalstore);
nreceivers = size(Rdata.coordinates,1);
nQ = Hsubmatrixdata.Qvecendnums(end);

P_receiver = zeros(nfrequencies,nreceivers,nsolutions);
P_nodiff2 = [];
propdiagnostics = [];
ivproblematicmatrix = zeros(nreceivers,100);

% The edge pairs where Qfinal has values, for some frequency and source
doesQsegmenthavevalues = zeros(size(Hsubmatrixdata.edgepairlist,1),1);
for ifreq = 1:nfrequencies
    for isol = 1:nsolutions
        doesQsegmenthavevalues = doesQsegmenthavevalues | edQsegmentflags(Qfinalstore{ifreq,isol},Hsubmatrixdata);
    end
end
doesQsegmenthavevalues = double(doesQsegmenthavevalues);
kvec = 2*pi*frequencies/envdata.cair;

% Qfirsttermstore is only needed for savediff2result = 1; otherwise an
% empty cell array is used in the loop, so that it is not copied to the
% parallel workers.
if savediff2 == 1
    Qfirsttermforloop = Qfirsttermstore;
else
    Qfirsttermforloop = {};
end

% The optional split into parts, by edge pair: the part number of each
% element of the Q vector
dosplit = ~isempty(splitlabels);
parts = [];
if dosplit
    npairs = size(Hsubmatrixdata.edgepairlist,1);
    splitlabels = double(splitlabels(:));
    if numel(splitlabels) ~= npairs || any(splitlabels < 1) || any(splitlabels ~= round(splitlabels))
        error('EDpropagateQvec:split','propsettings.splitlabels must give one part number (1,2,...) per edge pair')
    end
    nparts = max(splitlabels);
    partofQ = zeros(double(nQ),1);
    for ii = 1:npairs
        partofQ(Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) = splitlabels(ii);
    end
else
    nparts = 0;
    partofQ = [];
end

% For the 'accurate' propagation (the default), the receivers are
% handled in blocks by EDcalcpropagatematrix_block, which processes all
% the receivers of a block at once for each edge pair. The block size
% is chosen so that the arrays [blocksize,nQ] stay moderate (about
% 4e6 elements); with parallel workers, smaller blocks are used so that
% there are enough blocks to share. The 'legacy' propagation uses
% EDcalcpropagatematrix for one receiver at a time.
useblocks = any(strcmpi(propimplementation,{'accurate','exactpole'})) && isthinhole == 0;
if dosplit && ~useblocks
    error('EDpropagateQvec:split','The split into parts is only implemented for the ''accurate'' and ''exactpole'' propagation')
end
propagationblocksize = [];
if useblocks
    blocksize = max(1,min(nreceivers,floor(4e6/double(nQ))));
    if nworkers > 0
        blocksize = min(blocksize,max(1,ceil(nreceivers/40)));
    end
    blockstarts = 1:blocksize:nreceivers;
    nblocks = numel(blockstarts);
    propagationblocksize = blocksize;
    Pcell = cell(nblocks,1);
    P2cell = cell(nblocks,1);
    Partscell = cell(nblocks,1);
    propdiagnosticscell = cell(nblocks,1);

    parfor (ib = 1:nblocks, nworkers)
        recs = blockstarts(ib):min(blockstarts(ib)+blocksize-1,nreceivers);
        nb = numel(recs);
        if showtext >= 2
            disp(['      Receivers ',int2str(recs(1)),' to ',int2str(recs(end))])
        end
        [Fgeom,pathlengthsblock,propdiagnosticscell{ib}] = EDcalcpropagatematrix_block(edgedata,edgetoedgedata,...
            Hsubmatrixdata,Rdata.vispartedgesfromr(:,recs),...
            reshape(Rdata.rRsho(Rdata.reftoshortlistR(:,recs)),[],nb),...
            reshape(Rdata.thetaRsho(Rdata.reftoshortlistR(:,recs)),[],nb),...
            reshape(Rdata.zRsho(Rdata.reftoshortlistR(:,recs)),[],nb),doesQsegmenthavevalues,polemethod);
        % Only the columns with non-zero elements are used
        ivnonzero = find(any(Fgeom,1));
        Fgeom = Fgeom(:,ivnonzero);
        pathlengthsblock = pathlengthsblock(:,ivnonzero) - Rstartvalue;
        Pblock = zeros(nfrequencies,nb,nsolutions);
        P2block = zeros(nfrequencies,nb,nsolutions);
        if dosplit
            Partsblock = zeros(nfrequencies,nb,nsolutions,nparts);
            pnz = partofQ(ivnonzero);
        else
            Partsblock = [];
            pnz = [];
        end
        for ifreq = 1:nfrequencies
            Fk = Fgeom.*exp(-1i*kvec(ifreq)*pathlengthsblock);    % [nb x nnonzero]
            if nsolutions == 1 || doallSRcombinations == 1
                for isol = 1:nsolutions
                    Qfinalone = Qfinalstore{ifreq,isol};
                    Pblock(ifreq,:,isol) = Fk*Qfinalone(ivnonzero);
                    if savediff2 == 1
                        Qfirsttermone = Qfirsttermforloop{ifreq,isol};
                        P2block(ifreq,:,isol) = Fk*(Qfinalone(ivnonzero)-Qfirsttermone(ivnonzero));
                    end
                    if dosplit
                        Qnz = Qfinalone(ivnonzero);
                        for jp = 1:nparts
                            inp = (pnz == jp);
                            Partsblock(ifreq,:,isol,jp) = Fk(:,inp)*Qnz(inp);
                        end
                    end
                end
            else
                % doallSRcombinations = 0: only source N to receiver N
                for jr = 1:nb
                    jsol = recs(jr);   % (not isol: parfor does not allow assigning a for-loop variable)
                    if jsol <= nsolutions
                        Qfinalone = Qfinalstore{ifreq,jsol};
                        Pblock(ifreq,jr,jsol) = Fk(jr,:)*Qfinalone(ivnonzero);
                        if savediff2 == 1
                            Qfirsttermone = Qfirsttermforloop{ifreq,jsol};
                            P2block(ifreq,jr,jsol) = Fk(jr,:)*(Qfinalone(ivnonzero)-Qfirsttermone(ivnonzero));
                        end
                        if dosplit
                            Qnz = Qfinalone(ivnonzero);
                            for jp = 1:nparts
                                inp = (pnz == jp);
                                Partsblock(ifreq,jr,jsol,jp) = Fk(jr,inp)*Qnz(inp);
                            end
                        end
                    end
                end
            end
        end
        Pcell{ib} = Pblock;
        P2cell{ib} = P2block;
        Partscell{ib} = Partsblock;
    end

    % Collect the results of all the blocks
    propdiagnostics = struct('npairs',0,'npairs_eligible',0,'npairs_corrected',0, ...
        'nrows',0,'nrows_corrected',0,'correctedpairs',zeros(0,4));
    if strcmp(polemethod,'exact')
        propdiagnostics.nrows_exactpole = 0;
        propdiagnostics.nrows_fallback = 0;
    end
    if savediff2 == 1
        P_nodiff2 = zeros(size(P_receiver));
    end
    if dosplit
        parts = zeros(nfrequencies,nreceivers,nsolutions,nparts);
    end
    for ib = 1:nblocks
        recs = blockstarts(ib):min(blockstarts(ib)+blocksize-1,nreceivers);
        P_receiver(:,recs,:) = Pcell{ib};
        if savediff2 == 1
            P_nodiff2(:,recs,:) = P2cell{ib};
        end
        if dosplit
            parts(:,recs,:,:) = Partscell{ib};
        end
        pd = propdiagnosticscell{ib};
        propdiagnostics.npairs = propdiagnostics.npairs + pd.npairs;
        propdiagnostics.npairs_eligible = propdiagnostics.npairs_eligible + pd.npairs_eligible;
        propdiagnostics.npairs_corrected = propdiagnostics.npairs_corrected + pd.npairs_corrected;
        propdiagnostics.nrows = propdiagnostics.nrows + pd.nrows;
        propdiagnostics.nrows_corrected = propdiagnostics.nrows_corrected + pd.nrows_corrected;
        propdiagnostics.correctedpairs = [propdiagnostics.correctedpairs; pd.correctedpairs];
        if strcmp(polemethod,'exact')
            propdiagnostics.nrows_exactpole = propdiagnostics.nrows_exactpole + pd.nrows_exactpole;
            propdiagnostics.nrows_fallback = propdiagnostics.nrows_fallback + pd.nrows_fallback;
        end
    end
    clear Pcell P2cell Partscell propdiagnosticscell
else
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
            P_nodiff2(:,irec,:) = reshape(P2cell{irec},nfrequencies,1,nsolutions);
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

propinfo = struct;
propinfo.propdiagnostics = propdiagnostics;
propinfo.propagationblocksize = propagationblocksize;
propinfo.ivproblematicmatrix = ivproblematicmatrix;
propinfo.doesQsegmenthavevalues = doesQsegmenthavevalues;
propinfo.implementationnote = implementationnote;
if dosplit
    propinfo.parts = parts;
end

end

%--------------------------------------------------------------------------
function val = edgetsetting(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
else
    val = default;
end
end

%--------------------------------------------------------------------------
function doesQsegmenthavevalues = edQsegmentflags(Qfinal,Hsubmatrixdata)
% For each edge pair: 1 if any of its elements in Qfinal is non-zero.
nedgepairs = length(Hsubmatrixdata.edgepairlist(:,1));
doesQsegmenthavevalues = zeros(nedgepairs,1);
for iicounter = 1:nedgepairs
    sectionofQvec = (Hsubmatrixdata.Qvecstartnums(iicounter):Hsubmatrixdata.Qvecendnums(iicounter));
    doesQsegmenthavevalues(iicounter) = any(Qfinal(sectionofQvec));
end
end
