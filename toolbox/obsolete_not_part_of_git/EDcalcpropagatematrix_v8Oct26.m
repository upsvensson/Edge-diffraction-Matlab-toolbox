function [Fmatrix,ivproblematic,propdiagnostics] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
isthinhole,vispartedgesfromr,frequency,Rstart,rRvec,thetaRvec,zRvec,doesQsegmenthavevalues,...
showtext,ESIEpropagationimplementation)
% EDcalcpropagatematrix - propagates the sound pressure at the edges to a receiver
% position.
%
% Input parameters:
%   envdata,edgedata,edgetoedgedata,Hsubmatrixdata   Structs
%   isthinhole
%   vispartedgesfromr
%   frequency
%   Rstart
%   rRvec, thetaRvec, zRvec
%   doesQsegmenthavevalues
%   showtext (optional) 0 -> No text displayed on screen. Default: 0
%   ESIEpropagationimplementation (optional) 'accurate' (default) or 'legacy'.
%       'legacy'   : original, byte-identical code path. Near the
%                    propagation singularity (where the far-edge GL node's
%                    apex point za(z1) falls between two z2 GL nodes and
%                    the Lorentzian half-width B1(z1) is much narrower
%                    than the local z2 node spacing) the plain fixed-node
%                    GL sum can badly under/over-resolve that node's
%                    contribution. This path only *detects* that via
%                    ivproblematic; it does not correct it.
%       'accurate' : same result away from the singularity, but for
%                    flagged (i,k) far-edge-node/beta-term combinations,
%                    replaces the plain GL weight on the z2 nodes with a
%                    modified weight computed via Neumann's formula for
%                    Legendre Q-functions, so that the near-singular
%                    Lorentzian part of that beta term is integrated
%                    (semi-)analytically against the still-only-known-at-
%                    fixed-nodes source amplitude, instead of naively
%                    summed. NOTE: this name intentionally does not say
%                    "fast" -- unlike the EDwedge1st_fd toggle, the point
%                    of this mode is correctness near the propagation
%                    singularity, not speed; it does a little more work
%                    than 'legacy', only for the flagged entries.
%
% Output parameters:
%   Fmatrix
%   ivproblematic
%   propdiagnostics
%
% Uses EDcoordtrans1.
%
% Peter Svensson (peter.svensson@ntnu.no) 8 Oct. 2026
%
% [Fmatrix,ivproblematic,propdiagnostics] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
%     isthinhole,vispartedgesfromr,frequency,Rstart,...
%     rRvec,thetaRvec,zRvec,doesQsegmenthavevalues,showtext,ESIEpropagationimplementation);

% 20121025 Added the isthinhole parameter
% 20130131 Cleaned up
% 20130313 Reintroduced thetae2sho to handle non-convex geometries
% 20140207 Removed the factor 1/2 for alongplane radiation - moved to
%           calcedgeinteqmatrix!
% 20141201 Added the input parameter edgerelatedcoordsysmatrices
% 20141202 Replaced the input parameter edgerelatedcoordsysmatrices with
%          rRvec,thetaRvec,zRvec
% 20160128 Added the output parameter ivproblematic: a list of entries
%          where the integrand is close to singular.
% 20160129 Added the input parameter doesQsegmenthavevalues which is 0 if one
%          entire segment (corresponding to one edge-to-edge example) of
%          the edge source signal is zero. The Fmatrix for that specific
%          segment doesn't need to be computed then.
% 29 April 2016 Fixed a bug, handling the edge-to-edge propagation between
%          non-parallel edges.
% 24 Oct. 2017 Made a small correction for skewed edge pairs
% 17 Nov. 2017 Fixed a bug; could not handle more than 255 points per edge.
% 27 Nov. 2017 Copied from ESIE2toolbox
% 28 Nov. 2017 Introduced the non-global input parameter showtext
% 13 Apr 2018 Fixed a small bug; thinplaneboostvec got the format uint8 by
% mistake. Found by Antoine.
% 2 Nov 2023 Changed that showtext must be 3 instead of 2 to show all the
% annoying "Building...."
% 11 Feb 2025 Made a small step towards handling singularities, by
% explicitly calculating the absnyfivec.
% 27 Aug. 2026 Code by Claude, prompted by PS
% Added the 'implementation' input ('legacy'/'accurate') to
%          fix (not just flag) the propagation singularity for the
%          isthinhole==0 case, via a per-far-edge-node apex point za(z1)
%          and Lorentzian half-width B1(z1) (derived from the same small-
%          z_rel expansion as EDwedge1st_fd's near-apex approximation,
%          applied per virtual source z1 instead of once), and a
%          modified-GL-weight (Neumann Q-function) product-integration
%          correction on the existing, fixed z2 GL nodes. See the
%          'accurate' branch below and its local functions for details.
%          The isthinhole~=0 branch is NOT covered by 'accurate' yet
%          (falls back to the legacy formula unconditionally) -- flagged,
%          not fixed, since that singularity structure hasn't been worked
%          out.
% 20 Sep. 2026 Changed optional input parameter name to
% "ESIEpropagationimplementation".
% 25 Sep. 2026 Adapted to the change in orientation of the sparse
% quadraturematrix data.
% 8 Oct. 2026 Modified the new "accurate" implementation to be used more
% often than the first version. Set "accurate" to the default.
% 8 Oct 2026 A rewriting by Claude to speed things up.

if nargin < 13
    showtext = 0;
    ESIEpropagationimplementation = 'accurate';    
end
if nargin < 14 || isempty(ESIEpropagationimplementation)
    ESIEpropagationimplementation = 'accurate';
end

% ---------------------------------------------------------------
% TEMPORARY DIAGNOSTIC -- a persistent call counter that (as long as this
% function is called exactly once per receiver per frequency, matching
% how EDintegralequation_convex_tf drives your 629-receiver/1-frequency
% sweep) equals the receiver number, without needing to thread a receiver
% index through the whole call chain. Used below to gate a detailed
% corrected-vs-uncorrected delta printout for specific receivers of
% interest (edit DEBUG_TARGET_RECEIVERS as needed).
% IMPORTANT: persistent variables survive across calls within the same
% MATLAB session, so run "clear EDcalcpropagatematrix" (or "clear
% functions"/"clear all") before each fresh sweep, or this will keep
% counting up from wherever the previous run left off.
% Set DEBUG_TARGET_RECEIVERS = [] to silence this once you're done.
% ---------------------------------------------------------------
% persistent EDcalcpropagatematrix_callcount
% if isempty(EDcalcpropagatematrix_callcount)
%     EDcalcpropagatematrix_callcount = 0;
% end
% EDcalcpropagatematrix_callcount = EDcalcpropagatematrix_callcount + 1;
% receivernumber_fordebug = EDcalcpropagatematrix_callcount;
% DEBUG_TARGET_RECEIVERS = [211 577];
% debugtarget_fordebug = ismember(receivernumber_fordebug,DEBUG_TARGET_RECEIVERS);

% ----------------------------------------------------------------------------------------------
% Create the F-matrix by first creating a matrix showing which
% submatrices should be included.

 nQ = Hsubmatrixdata.Qvecendnums(end);
Fmatrix = zeros(1, nQ);
refto = edgetoedgedata.reftoshortlistE;

k = 2*pi*frequency/envdata.cair;

if showtext >= 3
    disp(' ')
    disp(['      Building the F matrix of size (1,',int2str( nQ),')'])
end

n1previous = 0;
n2previous = 0;
edge2previous = 0;
edge1previous = 0;

alledgepairsthin = prod(double(Hsubmatrixdata.isthinplaneedgepair));
noedgepairsthin  = prod(1-double(Hsubmatrixdata.isthinplaneedgepair));

ivrelevantedgepairs = find( vispartedgesfromr(Hsubmatrixdata.edgepairlist(:,1)) > 0 );

ivproblematic = [];
propdiagnostics = struct('npairs',0,'npairs_eligible',0,'npairs_corrected',0, ...
                  'nrows',0,'nrows_corrected',0,'correctedpairs',zeros(0,4));

for iicounter = 1:length(ivrelevantedgepairs)

    ii = ivrelevantedgepairs(iicounter);

    if doesQsegmenthavevalues(ii) == 1

        edge2 = Hsubmatrixdata.edgepairlist(ii,1);
        edge1 = Hsubmatrixdata.edgepairlist(ii,2);

        n1 = Hsubmatrixdata.nedgeelems(edge1);
        n2 = Hsubmatrixdata.nedgeelems(edge2);
    %     len1 = edgelengthvec(edge1);
        len2 = edgedata.edgelengthvec(edge2);

        if alledgepairsthin == 1
            thinplaneboost = 2;
        else
            if noedgepairsthin == 1
                thinplaneboost = 1;
            else
                thinplaneboost = double(Hsubmatrixdata.isthinplaneedgepair(ii)+1);
            end
        end

        ivaddition = [];

    %     if vispartedgesfromr(edge2) > 0
            if showtext >= 3
                disp(['      From edge ',int2str(edge1),' via edge ',int2str(edge2),' to R with thinplaneboost = ',int2str(thinplaneboost)])
                disp(['         Startcol in F-matrix   = (',int2str(Hsubmatrixdata.Qvecstartnums(ii)),')'])
                disp(['         Endcol in F-matrix   = (',int2str(Hsubmatrixdata.Qvecendnums(ii)),')'])
            end

            ny = pi/(2*pi-edgedata.closwedangvec(edge2));

            if n2 ~= n2previous
%                n2vec = Hsubmatrixdata.quadraturematrix_pos(n2,1:n2);
                n2vec = Hsubmatrixdata.quadraturematrix_pos(1:n2,n2);
 %               n2vec = n2vec(:);
%                weightvec2 = Hsubmatrixdata.quadraturematrix_weights(n2,1:n2);
                weightvec2 = Hsubmatrixdata.quadraturematrix_weights(1:n2,n2);
%                weightvec2 = weightvec2(:);

                if edge2 ~= edge2previous
                    dzvec2 = weightvec2*edgedata.edgelengthvec(edge2);
                    ze2_re2 = (0 + n2vec*len2);
                    edge2previous = edge2;
                end

            end

            if n1 ~= n1previous
%                n1vec = Hsubmatrixdata.quadraturematrix_pos(n1,1:n1);
               n1vec = Hsubmatrixdata.quadraturematrix_pos(1:n1,n1);
%                n1vec = n1vec(:);
%                weightvec1 = Hsubmatrixdata.quadraturematrix_weights(n1,1:n1);
                weightvec1 = Hsubmatrixdata.quadraturematrix_weights(1:n1,n1);
%                weightvec1 = weightvec1(:);
                if edge1 ~= edge1previous
                    dzvec1 = weightvec1*edgedata.edgelengthvec(edge1);
                    edge1previous = edge1;
                end

            end

            ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );
            re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );

%             edgecoords = [edgedata.edgestartcoords(edge2,:);edgedata.edgeendcoords(edge2,:)];
            rR = rRvec(edge2);
            thetaout = thetaRvec(edge2);
            zR = zRvec(edge2);

            thetain = edgetoedgedata.thetae1sho(refto(edge1,edge2));
            thetain_end = edgetoedgedata.thetae2sho(refto(edge1,edge2));
            if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_end) < 1e-10
                acrossface_in = 1;
            else
                acrossface_in = 0;

                % Reconvert the cylindrical coordinates of
                % edge1-re2 into cartesian! Then we can use the ESIE2coordtrans.
                % We define our own cartesian coord syst such that the reference
                % edge (edge2) has its starting point in [0 0 0], the z-axis along its
                % own edge, and the x-axis is along (onto) the reference plane.

                xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
                yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
                xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
                ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);
%                 thetae1_re2 = atan2(ye1_re2,xe1_re2);
%                 thetae1_re2 = thetae1_re2 + 2*pi*(thetae1_re2<0);

                [re1_re2,thetae1_re2,ze1_re2] = EDcoordtrans1([xe1_re2 ye1_re2 ze1_re2],[0 0 0;0 0 len2],[0 1 0]);

            end

            % Build vertical [n2,n1] matrices and expand them horizontally
            if n2 ~= n2previous || n1 ~= n1previous
                if n2 <= 255
                    n2vertmat = uint8(1:n2);
                else
                    n2vertmat = uint16(1:n2);
                end
                n2vertmat = reshape(repmat(n2vertmat,n1,1),n1*n2,1);
                n2previous = n2;

                if n1 <= 255
                    n1vertmat = uint8(1:n1).';
                else
                     n1vertmat = uint16(1:n1).';
                end
                n1vertmat = repmat(n1vertmat,n2,1);
                n1previous = n1;
            end

            % New 11 Feb. 2025 by PS: create absnyfivec for each
            % edge-to-edge pair, as a preparation for handling
            % singularities in the propagation integral.
            fivec = pi + thetain*[1 1 -1 -1] + thetaout*[1 -1 1 -1];
            absnyfivec = abs(ny*fivec);

            ldist = (zR - ze2_re2(n2vertmat)).^2   + rR.^2;
            ldist = sqrt(ldist);

            mdist = (ze1_re2(n1vertmat) - ze2_re2(n2vertmat)).^2   + re1_re2(n1vertmat).^2;
            mdist = sqrt(mdist);

            dz_expjkmloverml = exp(-1i*k*(ldist+mdist-Rstart))./ldist./mdist.*dzvec2(n2vertmat).*dzvec1(n1vertmat);

            if strcmpi(ESIEpropagationimplementation,'legacy')

                ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
                ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
                ch = ( ch + 1./ch)/2;

                % 28 March 2011: thetain is always zero (or = theta_wedge) so terms
                % 1 and 3 will always be identical; and terms 2 and 4.

                if isthinhole == 0
                    if acrossface_in == 1
                        beta = 2* ( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))) + ...
                             sin(ny*(pi + thetain - thetaout  ))./(ch - cos(ny*(pi + thetain - thetaout   ))));
                        % The lines below were added after Claude added the
                        % 'accuracy' mode.
                        if abs(thetain - thetaout - pi) < 0.01 || abs(thetain - thetaout + pi) < 0.01
                            % disp([num2str(thetain),', ',num2str(thetaout)])
                            rhoangle = rR./re1_re2(n1vertmat);
                            R02 = (re1_re2(n1vertmat) + rR).^2 + (zR - ze1_re2(n1vertmat)).^2;
                            B1beta = 4*sin(ny*(pi + thetain + thetaout   )/2).^2*rhoangle.^2.*R02/ny^2/(1+rhoangle).^4;
                            B0beta = sin(ny*(pi + thetain + thetaout   ))*2*rhoangle.^2.*R02/ny^2./(1+rhoangle).^4;
                            zapex = (zR*re1_re2(n1vertmat) + ze1_re2(n1vertmat)*rR)./(re1_re2(n1vertmat)+rR);
                            betaapprox = B0beta./((ze2_re2(n2vertmat)-zapex).^2 + B1beta);
                            % savetemp
                            % pause
                        end

                    else
                        beta =  sin(ny*(pi + thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) + thetaout   ))) + ...
                                sin(ny*(pi + thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) - thetaout   ))) + ...
                                sin(ny*(pi - thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) + thetaout   ))) + ...
                                sin(ny*(pi - thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) - thetaout   )));
                    end

                    % SM - 2 last condictions added
                    if cos(ny*(pi + thetain + thetaout   )) > 0.99 || cos(ny*(pi + thetain - thetaout   )) > 0.99 %|| cos(ny*(pi - thetain + thetaout   )) > 0.99 || cos(ny*(pi - thetain - thetaout   )) > 0.99
                        ivaddition = find(ch<1.01);
    %                    filename_beta = [Filepath,Filestem,'_beta_function.mat'];
    %                    eval(['save ',filename_beta,' beta'])
                    end

                    % Change 7 Feb. 2014
        %              Fsub =             -(2-acrossface_in)*ny/4/2/pi*beta.*dz_expjkmloverml*thinplaneboost; disp('Prop OLD')
                    Fsub = -ny/4/2/pi*beta.*dz_expjkmloverml*thinplaneboost;  % disp('Prop NEW')
        %             acrossface_in


                else % ny must be 1/2 to end up here
        %            dbetadthetaS = -4*ch*cos(thetaout/2).*(cos(thetaout) + 2*ch.^2 - 3)./(cos(thetaout)^2 - 2*ch.^2 + 1)./(1+j*k*mdist);
        %            Fsub = -1/8/2/pi*dbetadthetaS.*dz_expjkmloverml*thinplaneboost;   % Version 0
                    dbetadthetaS = -4*ch*cos(thetaout/2).*(cos(thetaout) + 2*ch.^2 - 3)./(cos(thetaout)^2 - 2*ch.^2 + 1);
                    Fsub = -1/8/2/pi*dbetadthetaS.*dz_expjkmloverml./mdist*thinplaneboost;   % Version 1
                end

            else % 'accurate'

                if isthinhole == 0 && acrossface_in == 1
                    % Validated case: far edge coplanar with edge2 (the
                    % geometry guaranteed for convex polyhedra, currently
                    % the only fully-supported tier). Against a fine
                    % independent reference, this consistently reduced
                    % error substantially (e.g. 24%->1.8%, 5.8%->0.03%
                    % across two grid resolutions) and matched legacy to
                    % ~1e-18 in a case with nothing to correct -- no
                    % regression found.
                    [Fsub,ivaddition,diagsub] = EDcalcpropagatematrix_accuratebeta( ...
                        ny,thetain,thetaout,acrossface_in,[], ...
                        ze1_re2,re1_re2,ze2_re2,zR,rR,len2,n1vertmat,n2vertmat, ...
                        mdist,dz_expjkmloverml,dzvec1,dzvec2,thinplaneboost,edge1,edge2, ...
                        k,Rstart);
                        propdiagnostics.npairs_eligible = propdiagnostics.npairs_eligible + 1;
                        propdiagnostics.nrows           = propdiagnostics.nrows + diagsub(2);
                        propdiagnostics.nrows_corrected = propdiagnostics.nrows_corrected + diagsub(1);
                        if diagsub(1) > 0
                            propdiagnostics.npairs_corrected = propdiagnostics.npairs_corrected + 1;
                            propdiagnostics.correctedpairs(end+1,:) = [edge1 edge2 diagsub];
                        end                    
%                        receivernumber_fordebug,debugtarget_fordebug,k,Rstart);
                   elseif isthinhole == 0 % acrossface_in == 0
                    % NOT validated -- falls back to legacy unconditionally.
                    % On a synthetic non-coplanar (skewed far edge) test
                    % case, the same correction sometimes reduced error
                    % substantially (43%->5%) but in another, only
                    % moderately near-singular case, it made the result
                    % WORSE than legacy (8%->12-13%, robust across several
                    % threshold settings) -- i.e. the "degrades gracefully
                    % away from the true singularity" assumption behind
                    % the I_1+I_2 split did not hold there. Root cause not
                    % yet understood; needs more work (and a real
                    % EDcoordtrans1-based test geometry from the toolbox,
                    % not just the synthetic one used so far) before this
                    % branch can be trusted. See chat history / delivery
                    % notes for the validation details.
                    ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
                    ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
                    ch = ( ch + 1./ch)/2;
                    beta =  sin(ny*(pi + thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) + thetaout   ))) + ...
                            sin(ny*(pi + thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) - thetaout   ))) + ...
                            sin(ny*(pi - thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) + thetaout   ))) + ...
                            sin(ny*(pi - thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) - thetaout   )));
                    if cos(ny*(pi + thetain + thetaout   )) > 0.99 || cos(ny*(pi + thetain - thetaout   )) > 0.99
                        ivaddition = find(ch<1.01);
                    else
                        ivaddition = [];
                    end
                    Fsub = -ny/4/2/pi*beta.*dz_expjkmloverml*thinplaneboost;
                else
                    % isthinhole~=0 (ny==1/2): the za(z1)/B1(z1) derivation
                    % used for the isthinhole==0 case above hasn't been
                    % worked out for this branch's dbetadthetaS formula.
                    % Flagged, not fixed -- falls back to the legacy
                    % formula unconditionally.
                    ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
                    ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
                    ch = ( ch + 1./ch)/2;
                    dbetadthetaS = -4*ch*cos(thetaout/2).*(cos(thetaout) + 2*ch.^2 - 3)./(cos(thetaout)^2 - 2*ch.^2 + 1);
                    Fsub = -1/8/2/pi*dbetadthetaS.*dz_expjkmloverml./mdist*thinplaneboost;
                    ivaddition = [];
                end

            end

            if showtext >= 3
                % 11 Feb 2025 by PS: fivec was moved to always be
                % calculated, not only of showtext >= 3.
%                 fivec = pi + thetain*[1 1 -1 -1] + thetaout*[1 -1 1 -1];
%                 absnyfivec = abs(ny*fivec);
                useserialexp1 = absnyfivec < 0.1;
                useserialexp2 = abs(absnyfivec - 2*pi) < 0.1;
                useserialexp = useserialexp1 | useserialexp2;
                if any(useserialexp)
                    disp('WARNING! Analytical approximation should be used')
                    format long
                    absnyfivec
                    figure(77)
                    plot(abs(Fsub),'-o')
                    pause
                end
            end

            propdiagnostics.npairs = propdiagnostics.npairs + 1;
            Fmatrix(Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) =     Fmatrix(Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) + Fsub.';

            if ~isempty(ivaddition)
                sectionofQvec = [Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)].';
                ivproblematic = [ivproblematic;sectionofQvec(ivaddition)];
            end
    end
end

end


% ========================================================================
% Local functions for the 'accurate' implementation.
% ========================================================================

function [Fsub,ivaddition,diagsub] = EDcalcpropagatematrix_accuratebeta( ...
        ny,thetain,thetaout,acrossface_in,thetae1_re2, ...
        ze1_re2,re1_re2,ze2_re2,zR,rR,len2,n1vertmat,n2vertmat, ...
        mdist,dz_expjkmloverml,dzvec1,dzvec2,thinplaneboost,edge1,edge2, ...
        k,Rstart)
%        receivernumber_fordebug,debugtarget_fordebug,k,Rstart)
    % Reproduces the legacy beta/Fsub computation term-by-term (still
    % using the exact closed-form beta terms everywhere), then, for
    % far-edge-node/beta-term combinations (i,k) flagged as near-singular,
    % replaces that term's contribution with a modified-GL-weight
    % (product-integration) treatment on the EXISTING z2 GL nodes, using
    % the per-i "virtual source" za(z1_i)/B1(z1_i) Lorentzian
    % approximation of that beta term near its own near-singular point.
    %
    % See the za(z1)/B1(z1) derivation: for fixed far-edge node i, edge2's
    % z2-direction near-singularity has exactly the same m(z2),l(z2),ch(z2)
    % structure as the single-edge (EDwedge1st_fd-style) problem, with a
    % virtual source at (r=re1_re2(i), z=ze1_re2(i)) and the fixed
    % receiver (rR,zR). This gives, for each of the (2 or 4) beta terms k:
    %   term_k(z2) ~= num_ik(i,k) / ( (z2-za(i))^2 + B1_ik(i,k) )   near za(i)
    % which is then integrated against the still-only-known-at-fixed-nodes
    % smooth factor h(z2)=kernelfactor(i,z2) (folded straight into the
    % modified weight, no special treatment needed for it -- it is fully
    % known, just evaluated at each node).
    %
    % Uses the I = I_1 + I_2 split: the exact term_k(i,j) values are kept
    % at every node (I_1, handled by the ordinary GL sum, unchanged from
    % legacy), and only the LORENTZIAN part of term_k near za(i) is
    % corrected (I_2): the plain-GL estimate of the Lorentzian-only
    % integral is subtracted off and replaced by the modified-GL estimate.
    % An earlier version replaced term_k(i,:) outright by the Lorentzian
    % for flagged (i,k); that was simpler but, when B1 wasn't extremely
    % small, could bias the tails (away from za(i)) enough to make the
    % result WORSE than legacy -- confirmed on a validation case. This
    % split fixes that: when the Lorentzian is a good local model, I_2's
    % correction is exactly what's needed near za(i) and I_1 is untouched
    % elsewhere; when B1 is only moderately small, the correction is small
    % by construction (plain and modified estimates of I_2 already agree),
    % so this degrades gracefully instead of introducing new bias.

    n1 = length(re1_re2);
    n2 = length(ze2_re2);

    % --- exact (legacy-identical) ch and beta terms ---
    ldist = sqrt( (zR - ze2_re2(n2vertmat)).^2 + rR.^2 );  % same as caller's ldist
    ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
    ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
    ch = ( ch + 1./ch)/2;

    if acrossface_in == 1
        thetaS_i = thetain*ones(n1,1);
        s1 = [1 1];  s2 = [1 -1];  termboost = 2;
        term1 = sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   )));
        term2 = sin(ny*(pi + thetain - thetaout  ))./(ch - cos(ny*(pi + thetain - thetaout   )));
        betaterms = {term1,term2};
        beta = termboost*(term1+term2);
    else
        thetaS_i = thetae1_re2(:);
        s1 = [1 1 -1 -1];  s2 = [1 -1 1 -1];  termboost = 1;
        term1 = sin(ny*(pi + thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) + thetaout   )));
        term2 = sin(ny*(pi + thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi + thetae1_re2(n1vertmat) - thetaout   )));
        term3 = sin(ny*(pi - thetae1_re2(n1vertmat) + thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) + thetaout   )));
        term4 = sin(ny*(pi - thetae1_re2(n1vertmat) - thetaout  ))./(ch - cos(ny*(pi - thetae1_re2(n1vertmat) - thetaout   )));
        betaterms = {term1,term2,term3,term4};
        beta = term1+term2+term3+term4; %#ok<NASGU> (kept for parity/debug; patched sum used below)
    end

    if cos(ny*(pi + thetain + thetaout   )) > 0.99 || cos(ny*(pi + thetain - thetaout   )) > 0.99
        ivaddition = find(ch<1.01);
    else
        ivaddition = [];
    end

    % --- geometry-only apex point / Lorentzian half-width per far-edge
    % node i, one column per beta term k (za(z1)/B1(z1) derivation) ---
    R0vec = sqrt( (ze1_re2 - zR).^2 + (re1_re2 + rR).^2 );            % [n1x1]
    zavec = ( zR*re1_re2 + ze1_re2*rR ) ./ ( re1_re2 + rR );          % [n1x1]
    kappavec = ny^2*(re1_re2+rR).^4 ./ (2*re1_re2.^2*rR^2.*R0vec.^2); % [n1x1]

    % NOTE: num_ik/B1_ik approximate each RAW (undoubled) term_k(z2), i.e.
    % exactly what betaterms{k} holds -- termboost must NOT be baked in
    % here. It is applied uniformly to both the exact and the modified-
    % weight contributions below, exactly mirroring how legacy applies it
    % once, uniformly, when forming beta = termboost*(term1+term2).
    phi_ik  = pi + thetaS_i*s1 + thetaout*s2;               % [n1 x nk]
    B1_ik   = 2*sin(ny*phi_ik/2).^2 ./ kappavec;             % [n1 x nk]
    num_ik  = sin(ny*phi_ik) ./ kappavec;                    % [n1 x nk]

    % --- flag which (i,k) entries are under-resolved by the existing,
    % fixed z2 grid ---
    if n2 > 1
        % Change on 7 Oct. 2026 as suggested by Claude
        %localspacing = min(diff(sort(ze2_re2)));
        localspacing = max(diff(sort(ze2_re2)));
    else
        localspacing = len2;
    end
    %zrelthreshold = 1.0;     % heuristic -- tune like the other resolution
                              % thresholds used elsewhere in this codebase.
                              % Deliberately generous: sqrt(B1) is the
                              % Lorentzian's half-width, and a fixed GL
                              % rule needs several nodes actually inside
                              % that half-width to resolve it well, so
                              % flagging only when sqrt(B1) is much
                              % *smaller* than one node spacing is already
                              % too late -- see the worked n1=n2=16 case in
                              % the accompanying validation notes.
     % Change on 7 Oct. 2026 as suggested by Claude
     zrelthreshold = 5;        
                              
    insiderange = (zavec >= -0.05*len2) & (zavec <= 1.05*len2);   % [n1x1]
    narrowenough = sqrt(max(B1_ik,0)) < zrelthreshold*localspacing;             % [n1 x nk]
    flag_ik = insiderange & narrowenough;  % [n1 x nk]
    diagsub = [nnz(any(flag_ik,2)) size(flag_ik,1)];   % [corrected rows, all rows]

    % ---------------------------------------------------------------
    % TEMPORARY DIAGNOSTIC -- prints ONE line per (i,k) the moment the
    % correction actually engages (flag_ik true), i.e. exactly the
    % "singularity handling kicks in here" signal, mirroring the disp()
    % you added in the legacy branch. Silent otherwise, so a 629-receiver
    % sweep only produces output for the receivers/edge-pairs where
    % something is actually being corrected. thetaout is the receiver's
    % angle in edge2's local frame -- tag your outer receiver loop with
    % its index/global angle if you need to map these lines back to a
    % specific one of your 629 positions. NOTE: flag_ik depends only on
    % geometry, not on frequency, so if you loop over frequencies inside
    % the same receiver/edge-pair, expect the identical line to repeat
    % once per frequency -- that's expected, not a bug.
    % Set singularitydebugprint=false to silence this once you're done.
    % ---------------------------------------------------------------
    singularitydebugprint = false;


    if singularitydebugprint && any(flag_ik(:))
        [flagi,flagk] = find(flag_ik);
        for hh = 1:numel(flagi)
            irow = flagi(hh); kcol = flagk(hh);
            fprintf(['[SINGULARITY HIT] edge1=%d edge2=%d  i=%d k=%d  ' ...
                     'thetain=%.4f thetaout=%.4f  za=%.4f (edge2 span [0,%.4g])  ' ...
                     'sqrt(B1)=%.4e (localspacing=%.4g)  rR=%.4g zR=%.4g\n'], ...
                    edge1,edge2,irow,kcol,thetain,thetaout,zavec(irow),len2, ...
                    sqrt(max(B1_ik(irow,kcol),0)),localspacing,rR,zR);
        end
    end

    % --- build the patched per-node beta*z2-weight sum, term by term ---
    dzvec2_row = dzvec2(:).';                     % [1 x n2]
    patched_betaweight = zeros(n1,n2);            % holds sum_k [dzvec2(j)*term_k(i,j)], patched

    if any(flag_ik(:))
        % 8 Oct. 2026: nodes, weights and the scaled Legendre matrix depend
        % only on n2, so they are computed once and cached.
        [xnodes_std,wnodes_std,Sstd] = glStandardCached(n2); %#ok<ASGLU>
    end

    ze2_row = ze2_re2(:).';                       % [1 x n2]

    for kk = 1:numel(betaterms)
        termk_mat = reshape(betaterms{kk},n1,n2);      % exact term_k(i,j), RAW/undoubled
        contrib = termboost * termk_mat .* dzvec2_row;  % legacy per-(i,j) contribution, [n1 x n2] -- LEFT UNCHANGED
        flaggedi = find(flag_ik(:,kk));

        if ~isempty(flaggedi)
            % 8 Oct. 2026: vectorised over all flagged rows of this term.
            % Same I_1 + I_2 split as before: the exact term_k(i,j) values
            % stay untouched; only the Lorentzian part is corrected, by
            % subtracting its plain-GL estimate and adding the modified-
            % weight (product-integration) estimate.
            zaf  = reshape(zavec(flaggedi),[],1);              % [nf x 1]
            B1f  = reshape(max(B1_ik(flaggedi,kk),0),[],1);    % [nf x 1]
            numf = reshape(num_ik(flaggedi,kk),[],1);          % [nf x 1]
            xa   = 2*zaf/len2 - 1;
            Bhat = 4*B1f/len2^2;
            WT = modifiedWeightsLorentzianVec(Sstd,wnodes_std,xa,Bhat);   % [n2 x nf]
            wt_physical = termboost*(2/len2) * (numf .* WT.');            % [nf x n2], modified-GL estimate
            lorentzian_at_nodes = numf ./ ((ze2_row - zaf).^2 + B1f);     % [nf x n2]
            plain_lorentzian_est = termboost * lorentzian_at_nodes .* dzvec2_row;
            contrib(flaggedi,:) = contrib(flaggedi,:) + (wt_physical - plain_lorentzian_est);
        end

        patched_betaweight = patched_betaweight + contrib;
    end

    kernelfactor = dz_expjkmloverml ./ (dzvec2(n2vertmat).*dzvec1(n1vertmat));  % [n1*n2 x 1]
    kernelfactor_mat = reshape(kernelfactor,n1,n2);

    Fsub_mat = -ny/4/2/pi*thinplaneboost * (dzvec1(:) .* patched_betaweight) .* kernelfactor_mat;
    Fsub = reshape(Fsub_mat,n1*n2,1);
end


function [x,w] = gaussLegendreNodes(n)
    % Golub-Welsch algorithm: standard n-point Gauss-Legendre nodes and
    % weights on [-1,1]. Verified to match EDdistelements(n,2)'s hardcoded
    % tables (mapped to [0,1]) to ~1e-16.
    if n == 1
        x = 0; w = 2; return
    end
    ivec = (1:n-1);
    betavec = ivec ./ sqrt(4*ivec.^2 - 1);
    J = diag(betavec,1) + diag(betavec,-1);
    [V,D] = eig(J);
    x = diag(D);
    [x,idx] = sort(x);
    V = V(:,idx);
    w = 2*(V(1,:).^2)';
    x = x(:); w = w(:);
end

function P = legendreVandermonde(x,Nmax)
    % P(j,n+1) = P_n(x(j)), n=0..Nmax-1
    x = x(:);
    M = length(x);
    P = zeros(M,Nmax);
    P(:,1) = 1;
    if Nmax > 1
        P(:,2) = x;
    end
    for n = 2:Nmax-1
        P(:,n+1) = ((2*n-1)*x.*P(:,n) - (n-1)*P(:,n-1))/n;
    end
end
function [x,w,S] = glStandardCached(n)
    % Standard n-point Gauss-Legendre nodes x and weights w on [-1,1]
    % (ascending), and S(j,m+1) = (2m+1)*P_m(x_j), m = 0..n-1.
    % These depend only on n, so they are cached between calls.
    persistent cachedn cachedx cachedw cachedS
    if isempty(cachedn)
        cachedn = [];  cachedx = {};  cachedw = {};  cachedS = {};
    end
    idx = find(cachedn == n,1);
    if isempty(idx)
        [x,w] = gaussLegendreNodes(n);
        P = legendreVandermonde(x,n);
        S = P .* (2*(0:n-1)+1);
        cachedn(end+1) = n;
        cachedx{end+1} = x;
        cachedw{end+1} = w;
        cachedS{end+1} = S;
    else
        x = cachedx{idx};  w = cachedw{idx};  S = cachedS{idx};
    end
end

function Q = neumannQ(z,Nmax)
    % Legendre functions of the second kind Q_0..Q_{Nmax-1} at the complex
    % values z (off the interval [-1,1]). Returns [nz x Nmax].
    % Vectorised over z, 8 Oct. 2026.
    %  - Poles close to [-1,1] (rho^(2*Nmax) < 1e4): upward recurrence,
    %    which is stable there.
    %  - Poles further away: downward (Miller) recurrence with a common
    %    padding for the whole group, rescaled against the closed form Q_0.
    zcol = z(:);
    nz = length(zcol);
    Q = zeros(nz,Nmax);
    if nz == 0
        return
    end
    Q0true = 0.5*log((zcol+1)./(zcol-1));
    rho = abs(zcol + sqrt(zcol.^2-1));
    rho = max(rho,1./rho);
    useupward = (2*Nmax)*log(rho) < log(1e4);

    % --- upward recurrence ---
    iv = find(useupward);
    if ~isempty(iv)
        zz = zcol(iv);
        QA = zeros(length(iv),Nmax);
        QA(:,1) = Q0true(iv);
        if Nmax > 1
            QA(:,2) = zz.*Q0true(iv) - 1;
        end
        for n = 1:Nmax-2
            QA(:,n+2) = ( (2*n+1)*zz.*QA(:,n+1) - n*QA(:,n) )/(n+1);
        end
        Q(iv,:) = QA;
    end

    % --- downward (Miller) recurrence ---
    iv = find(~useupward);
    if ~isempty(iv)
        zz = zcol(iv);
        nb = length(iv);
        pad = ceil(37/(2*log(min(rho(iv))))) + 20;
        Ntop = Nmax + pad;
        QB = zeros(nb,Nmax);
        qnp2 = zeros(nb,1);        % plays the role of Qd(n+2)
        qnp1 = ones(nb,1);         % plays the role of Qd(n+1)
        for n = Ntop-1:-1:1
            qn = ( (2*n+1)*zz.*qnp1 - (n+1)*qnp2 )/n;
            qnp2 = qnp1;
            qnp1 = qn;
            if n <= Nmax
                QB(:,n) = qn;
            end
            % Guard against overflow for poles far from the interval
            ivbig = find(abs(qn) > 1e150);
            if ~isempty(ivbig)
                qnp1(ivbig) = qnp1(ivbig)*1e-150;
                qnp2(ivbig) = qnp2(ivbig)*1e-150;
                if n <= Nmax
                    QB(ivbig,n:Nmax) = QB(ivbig,n:Nmax)*1e-150;
                end
            end
        end
        Q(iv,:) = QB .* (Q0true(iv)./QB(:,1));
    end
end

function WT = modifiedWeightsLorentzianVec(S,wnodes,za,B)
    % Modified weights for g(t) = 1/((t-za)^2 + B) on the standard GL
    % nodes, for a whole vector of (za,B) pairs at once.
    %   S       [N x N], S(j,m+1) = (2m+1)*P_m(x_j)   (from glStandardCached)
    %   wnodes  [N x 1] standard GL weights
    %   za,B    [nz x 1], real, B >= 0
    % Returns WT [N x nz], one column per (za,B) pair.
    %
    % With the pole p = za + i*sqrt(B):
    %   1/((t-za)^2+B) = Im( 1/(t-p) ) / sqrt(B)
    % and the single-pole weights follow from Neumann's formula
    %   int P_m(t)/(t-z) dt = -2 Q_m(z):
    %   w_j(z) = -w_j * sum_m (2m+1) P_m(x_j) Q_m(z)
    % The weights for the conjugate pole are the complex conjugate, so one
    % evaluation is enough.
    za = za(:);
    sB = sqrt(max(B(:),eps));
    N = size(S,1);
    Q = neumannQ(za + 1i*sB,N);                 % [nz x N]
    WT = -(S*imag(Q).') .* wnodes(:) ./ (sB.');  % [N x nz]
end
