function [Qfirstterm,soudiagnostics] = EDinteg_souterm(envdata,edgedata,edgetoedgedata,...
    Hsubmatrixdata,controlparameters,...
    vispartedgesfroms,vispartedgesfroms_start,vispartedgesfroms_end,frequency,gaussvectors,rSvec,thetaSvec,zSvec,showtext)
% EDinteg_souterm - Does the first step of the ED int eq: compute the 
% first term of the edge source signals.
%
% Input parameters:
%   envdata,edgedata,edgetoedgedata,
%   Hsubmatrixdata,controlparameters   structs
%   vispartedgesfroms
%   vispartedgesfroms_start
%   vispartedgesfroms_end
%   frequency
%   gaussvectors            A cell variable, with nodes and weights up to
%                           the max. ngauss that could occur.
%   rSvec,thetaSvec,zSvec
%   showtext (optional)     0 -> no text displayed on screen. Default: 0
%
% Optional field in the struct controlparameters (8 Oct. 2026):
%   .ESIEsourcetermimplementation   'legacy' (default) or 'accurate'.
%       'legacy'   : original code path, including the single-sample
%                    LCN replacement via EDintegratebetaoverm.
%       'accurate' : when the source is close to the plane of a face, the
%                    source term q0(z2->z3) has, for each fixed z3, a
%                    narrow Lorentzian peak along z2 which the GL nodes
%                    cannot resolve. The sampled values are then replaced
%                    by "effective" values, such that every later GL sum
%                    over z2 (in the H-matrix product and in the
%                    propagation to the receiver) integrates the
%                    Lorentzian part exactly against the polynomial
%                    through the nodal values of the remaining factors
%                    (product integration; this is a locally corrected
%                    Nystrom scheme, applied to all nodes of the edge).
%                    Same method as the 'accurate' mode of
%                    EDcalcpropagatematrix. Applied to edge pairs with
%                    non-skewed edges, a fully visible via-edge and
%                    Gauss-Legendre nodes; other edge pairs use 'legacy'.
%
% Output parameters:
%   Qfirstterm      Column vector with all edge source signals caused by
%                   the external source. All "via-edge" to "to-edge" combos
%                   are stacked after each other. For each such combo, the
%                   sequence is:
%                           (edge3)             (edge2)
%                       to-edge-element     via-edge-element
%                            1                      1
%                            1                      2
%                            1                      3
%                   etc
%   soudiagnostics  (optional) struct with fields
%                   .implementation     'legacy' or 'accurate'
%                   .npairs             edge pairs processed
%                   .npairs_eligible    pairs sent to the accurate branch
%                   .npairs_corrected   pairs with at least one corrected slice
%                   .nslices            z3-slices in the eligible pairs
%                   .nslices_corrected  z3-slices that were corrected
%                   .correctedpairs     [edge2 edge3 corrected_slices all_slices]
%                   .npairs_skew, .npairs_partlyvisible, .npairs_nodesnotGL
%                                       pairs left to the legacy code, by reason
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [Qfirstterm,soudiagnostics] = EDinteg_souterm(envdata,edgedata,...
%    edgetoedgedata,Hsubmatrixdata,controlparameters,vispartedgesfroms,...
%    vispartedgesfroms_start,vispartedgesfroms_end,frequency,gaussvectors,...
%    rSvec,thetaSvec,zSvec,showtext);

% 20120217 Functioning version
% 20121025 Added the isthinhole output parameter
% 20130129 Cleaned up thoroughly
% 20130130 Removed the isthinhole output parameter
% 20140207 Introduced the factor 1/2 scaling (and removed it from the
% propagation). Also introduced the input parameter vispartedgesfromIS
% 20140923 Fixed a bug for concave geometries and skewed edge pairs
% 20141201 Added the input parameter edgerelatedcoordsysmatrices
% 20141202 Replaced the input parameter edgerelatedcoordsysmatrices with
%          rSvec,thetaSvec,zSvec
% 24 Oct. 2017 Had not implemented the propagation from edge 2 to edge 3
%              for the case that the two edges are non-inplane
% 4 Nov. 2017 Implemented the partly-visible edge (take from the
%             ESIE2integ_souterm_ir version), which implied introducing the
%             new input parameters vispartedgesfroms_start and
%             vispartedgesfroms_end.
% 17 Nov. 2017 Fixed a bug: used uint8, so could not handle more than 255
%               points per edge.
% 27 Nov. 2017 Copied from ESIE2toolbox. Trimmed down, used more structs.
% 28 Nov. 2017 Introduced the non-global input parameter showtext
% 15 Dec. 2017 Experiments with detecting near singularities
% 7 Mar 2018 Hid the singularity information printouts behind "if showtext
% >= 2"
% 22 May 2018 Fixed a mistake "in2str" instead of "int2str".
% 29 May 2018 Implemented LCN for the highest sample (if it stands out).
% 2 Nov. 2023 Changed so that showtext = 3 is needed instead of 2 to show
% the "Startrow...." messages
% 20 Nov. 2023 Changed so that showtext = 3 is needed instead of 2 to show
% the "Building the G matrix" messages and showtext = 4 needed for the
% "Startrow" meassages.
% 24 Nov. 2023 The input parameter gaussvectors was changed to include all
% gauss order values up to the maximum one.
% 8 Oct. 2026 Code by Claude, prompted by PS: added the 'accurate'
% implementation (controlparameters.ESIEsourcetermimplementation) and the
% optional output soudiagnostics.
% 9 Oct. 2026 Code cleanup, changed bigmatrix... to Qvec...
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

if nargin < 14
    showtext = 0;
end

discretizationtype = controlparameters.discretizationtype; %#ok<NASGU>

useaccurate = 0;
if isfield(controlparameters,'ESIEsourcetermimplementation') && ~isempty(controlparameters.ESIEsourcetermimplementation)
    useaccurate = strcmpi(controlparameters.ESIEsourcetermimplementation,'accurate');
end
soudiagnostics = struct('implementation','legacy','npairs',0,'npairs_eligible',0,...
    'npairs_corrected',0,'nslices',0,'nslices_corrected',0,'correctedpairs',zeros(0,4),...
    'npairs_skew',0,'npairs_partlyvisible',0,'npairs_nodesnotGL',0);
if useaccurate
    soudiagnostics.implementation = 'accurate';
end

%  symmetrycompression = 1;

% ----------------------------------------------------------------------------------------------
% Create the G-matrix

G = zeros(Hsubmatrixdata.Qvecendnums(end),1);
refto = edgetoedgedata.reftoshortlistE;

k = 2*pi*frequency/envdata.cair;

if showtext >= 3
    disp(' ')
    disp('      Calculating the ESIE source term')%: symmetrycompression = ',int2str(symmetrycompression)])    
    disp(['      Building the G matrix of size (',int2str(Hsubmatrixdata.Qvecendnums(end)),',1)'])
end

% ----------------------------------------------------------------------------------------------
% Go through all the edge-pair combos

nyvec = pi./(2*pi-edgedata.closwedangvec);

n2previous = 0;
n3previous = 0;

for ii = 1:size(Hsubmatrixdata.edgepairlist,1)    
    edge2 = Hsubmatrixdata.edgepairlist(ii,2);
    
    if vispartedgesfroms(edge2) ~= 0
        edge3 = Hsubmatrixdata.edgepairlist(ii,1);
        if showtext >= 4
%         if showtext >= 1
            disp(' ')
            disp(['         From S via edge ',int2str(edge2),' to edge ',int2str(edge3)])
            disp(['         Startrow in G-matrix   = (',int2str(Hsubmatrixdata.Qvecstartnums(ii)),'). Endrow in G-matrix   = (',int2str(Hsubmatrixdata.Qvecendnums(ii)),')'])
%             nfrom = Hsubmatrixdata.nedgeelems(edge2);
%             nto   = Hsubmatrixdata.nedgeelems(edge3);
        end
        
        ny = nyvec(edge2);
        len2 = edgedata.edgelengthvec(edge2);
        n2 = Hsubmatrixdata.nedgeelems(edge2);
        n3 = Hsubmatrixdata.nedgeelems(edge3);

        % Changed 24 Nov. 2023: All gauss nodes and weights up to a max
        % value are stored in a cell variable
%         if size(gaussvectors,1) == n2
%             n2vec = gaussvectors(:,1);
%             weight2vec = gaussvectors(:,2);
%         else
%             disp('   Calling EDdistelements')
%            [n2vec,weight2vec] = EDdistelements(n2,discretizationtype);
%         end
        n2vec = gaussvectors{n2}(:,1);
        weight2vec = gaussvectors{n2}(:,2);
        % New 4 Nov. 2017: here we zero out the edge points that are not
        % visible.
        if vispartedgesfroms(edge2) ~= 1
            n2vec = n2vec.*(n2vec >= vispartedgesfroms_start(edge2)).*(n2vec <= vispartedgesfroms_end(edge2));
        end
        % Change 24 Nov. 2023 gaussvectors is now a cell variable
        % containing gauss nodes and weights for all orders up to the
        % maximum one.
%         if size(gaussvectors,1) == n3
%             n3vec = gaussvectors(:,1);
%         else
%            [n3vec,~] = EDdistelements(n3,discretizationtype);
%         end
        n3vec = gaussvectors{n3}(:,1);
        
        % FOund bug 23 Sep. 2014: did not take into account that the other
        % end of the edge might have a different theta value.
        thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
        thetaout_otherend = edgetoedgedata.thetae2sho(refto(edge3,edge2));
        
%         if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-closwedangvec(edge2))) < 1e-10)
        if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetaout-thetaout_otherend)<1e-10
            acrossface_out = 1;
        else
            acrossface_out = 0;
        end     
        if abs(thetaout - thetaout_otherend) < 1e-6
            skewedges = 0;
        else            
            skewedges = 1;
        end
        if showtext >= 4
            disp(['         Acrossface_out = ',int2str(acrossface_out)])
        end

%         edgecoords = [edgedata.edgestartcoords(edge2,:);edgedata.edgeendcoords(edge2,:)];
        rS = rSvec(edge2);
        thetaS = thetaSvec(edge2);
        zS = zSvec(edge2);
        
        singwarning = 0;
        singterm = [0 0 0 0];
        if thetaout == 0 && max([  abs(cos(ny*(pi + thetaS   )))  abs(cos(ny*(pi - thetaS   )))]) > 0.999
            singwarning = 1;
            if showtext >= 2
                disp(['      Singularity warning, location 1'])
            end
        end
        if thetaout ~= 0 && max([ abs(cos(ny*(pi + thetaS + thetaout   )))  abs(cos(ny*(pi + thetaS - thetaout   )))  ...
                          abs(cos(ny*(pi - thetaS + thetaout   )))  abs(cos(ny*(pi - thetaS - thetaout   )))  ]) > 0.999
            singwarning = 2;                       
            if showtext >= 2
                disp(['      Singularity warning, location 2'])
            end
        end
        
        % Build vertical [n3,n2] matrices and expand them horizontally
        if n3 ~= n3previous || n2 ~= n2previous 

            if n3 <= 255
                n3vertmat = uint8(1:n3);
            else
                n3vertmat = uint16(1:n3);
            end
%             n3vertmat = uint8(1:n3);
            n3vertmat = reshape(repmat(n3vertmat,n2,1),n2*n3,1);
            n3previous = n3;
            
             if n2 <= 255
                   n2vertmat = uint8(1:n2).';
             else
                   n2vertmat = uint16(1:n2).';                   
             end
%             n2vertmat = uint8(1:n2).';
             n2vertmat = repmat(n2vertmat,n3,1);        
             n2previous = n2;            
        end
                
        % z will always be linearly distributed whether edges are skewed or
        % not
        ze2_re2 = n2vec*len2;
%         ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2))   );
        ze3_start_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2));
        ze3_end_re2   = edgetoedgedata.ze2sho(refto(edge3,edge2));
        ze3_re2 = ze3_start_re2 + n3vec*( ze3_end_re2 - ze3_start_re2   );        

        % .... so mdist will be simple, also for skewed edges

        mdist = (zS - ze2_re2(n2vertmat)).^2   + rS.^2;
        mdist = sqrt(mdist);
        
%         dz_expjkmoverm = exp(-j*k*mdist)./mdist.*dzvec(n2vertmat);
        dz_expjkmoverm = exp(-1i*k*mdist)./mdist;

        % The r- and theta-values will not be linearly distributed for
        % skewed edges

        re3_start_re2 = edgetoedgedata.re1sho(refto(edge3,edge2));
        re3_end_re2   = edgetoedgedata.re2sho(refto(edge3,edge2));

        if skewedges == 0
        
    %         re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );
            re3_re2 = re3_start_re2 + n3vec*( re3_end_re2 - re3_start_re2   );

        else
            
             % First we need to reconvert the cylindrical coordinates of
            % edge2-re1 into cartesian! Then we can use the ESIE2coordtrans.
            % We define our own cartesian coord syst such that the reference
            % edge has its starting point in [0 0 0], and the x-axis is along
            % the reference plane.

            xe3_start_re2 = re3_start_re2*cos(thetaout);
            ye3_start_re2 = re3_start_re2*sin(thetaout);
            xe3_end_re2   = re3_end_re2*cos(thetaout_otherend);
            ye3_end_re2   = re3_end_re2*sin(thetaout_otherend);

            xe3_vec_re2 = xe3_start_re2 + (xe3_end_re2 - xe3_start_re2)*n3vec;
            ye3_vec_re2 = ye3_start_re2 + (ye3_end_re2 - ye3_start_re2)*n3vec;
            ze3_vec_re2 = ze3_start_re2 + (ze3_end_re2 - ze3_start_re2)*n3vec;

            % NB!!! Below we need to check if really the full lengths of both
            % edges should be used. If S or R see only part of the respective
            % edge then the edge-to-edge contribution should be "windowed" too.

            [re3_re2,thetae3_re2,ze3_re2] = EDcoordtrans1([xe3_vec_re2 ye3_vec_re2 ze3_vec_re2],[0 0 0;0 0 len2],[0 1 0]);
            
        end

        ldist = (ze3_re2(n3vertmat) - ze2_re2(n2vertmat)).^2   + re3_re2(n3vertmat).^2;
        ldist = sqrt(ldist);

        ch = ((ze2_re2(n2vertmat) - zS).*(ze2_re2(n2vertmat) - ze3_re2(n3vertmat)) + mdist.*ldist)./rS./re3_re2(n3vertmat);
        ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
        ch = ( ch + 1./ch)/2;               
            
        if thetaout == 0
           beta = 2*( sin(ny*(pi + thetaS  ))./(ch - cos(ny*(pi + thetaS   ))) + ...
                       sin(ny*(pi - thetaS  ))./(ch - cos(ny*(pi - thetaS   ))));        

         else

            if skewedges == 0

                beta = ( sin(ny*(pi + thetaS + thetaout  ))./(ch - cos(ny*(pi + thetaS + thetaout   ))) + ...
                         sin(ny*(pi + thetaS - thetaout  ))./(ch - cos(ny*(pi + thetaS - thetaout   ))) + ...
                         sin(ny*(pi - thetaS + thetaout  ))./(ch - cos(ny*(pi - thetaS + thetaout   ))) + ...
                         sin(ny*(pi - thetaS - thetaout  ))./(ch - cos(ny*(pi - thetaS - thetaout   ))));  
                    
            else

                 beta = ( sin(ny*(pi + thetaS + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetaS + thetae3_re2(n3vertmat)   ))) + ...
                 sin(ny*(pi + thetaS - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetaS - thetae3_re2(n3vertmat)   ))) + ...
                 sin(ny*(pi - thetaS + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetaS + thetae3_re2(n3vertmat)   ))) + ...
                 sin(ny*(pi - thetaS - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetaS - thetae3_re2(n3vertmat)   )))); 
             
            end
        end
                        
        Gsub = -(2-acrossface_out)*ny/4/pi*beta.*dz_expjkmoverm;
        
        soudiagnostics.npairs = soudiagnostics.npairs + 1;

        % 8 Oct. 2026: 'accurate' treatment of the near-singular source term
        accuratedone = 0;
        % 8 Oct. 2026 (later): "fully visible" is tested on the nodes
        % themselves (unchanged by the windowing above), not on the value
        % of vispartedgesfroms, which is not 1 for a fully visible edge.
        fullyvisible = isequal(n2vec,gaussvectors{n2}(:,1));
        if useaccurate
            soudiagnostics.npairs_skew = soudiagnostics.npairs_skew + (skewedges ~= 0);
            soudiagnostics.npairs_partlyvisible = soudiagnostics.npairs_partlyvisible + (fullyvisible == 0);
        end
        if useaccurate && skewedges == 0 && fullyvisible && n2 > 1
            Gsubfactor = -(2-acrossface_out)*ny/4/pi;
            [Gsub,nslicescorrected,accuratedone] = EDinteg_souterm_accurate(Gsub,ny,thetaS,thetaout,...
                rS,zS,re3_re2,ze3_re2,ze2_re2,weight2vec,len2,dz_expjkmoverm(1:n2),Gsubfactor);
            soudiagnostics.npairs_nodesnotGL = soudiagnostics.npairs_nodesnotGL + (accuratedone == 0);
            if accuratedone == 1
                soudiagnostics.npairs_eligible = soudiagnostics.npairs_eligible + 1;
                soudiagnostics.nslices = soudiagnostics.nslices + n3;
                soudiagnostics.nslices_corrected = soudiagnostics.nslices_corrected + nslicescorrected;
                if nslicescorrected > 0
                    soudiagnostics.npairs_corrected = soudiagnostics.npairs_corrected + 1;
                    soudiagnostics.correctedpairs(end+1,:) = [edge2 edge3 nslicescorrected n3];
                end
            end
        end

        if accuratedone == 0    % legacy treatment (unchanged)

          % For each value of z1 = the receiver edge point (which there
          % are n3 values of), select the z2-range (which there are n2
          % values of). We should integrate numerically over z2.
          n2countervec = [1:n2].';
          for jj = 1:n3
            ivselect = (jj-1)*n2 + n2countervec;
                        
            dGrel = (diff(abs(Gsub(ivselect)))./abs(Gsub(ivselect(1:end-1))));

            npeaks = 0;
            ivpeaks = [];
            dofixsingularity = 0;
            iv1 = find( abs(dGrel(1)+1) < 2e-3 );
            iv2 = find( abs(dGrel) > 0.8e3 );
            if any(iv1) || any(iv2)
                dofixsingularity = 1;
            end
                        
            % Now we need to compute Iaccurate = integrate beta from S,
            % via edge2, to edge3. Integration should be over the
            % entire edge2, to one point on edge3.

            if dofixsingularity == 1
               [maxval,ivpeaks] = max(abs(Gsub(ivselect)));
               npeaks = 1;
               if showtext >= 2
                    disp(['   Integrating source term, from S, via edge ',int2str(edge2),' to edge ',int2str(edge3)])
               end
                              
               % Source parameters: rS, thetaS, zS
               % Receiver parameters: re3_re2(jj),
               % thetaout, ze3_re2(jj)
               % Edgelength = [0,len2], ny
               % The function EDwedge1st_fd uses an integrand which is
               % -nu/(4*pi)*beta*exp(-jkm)/m*exp(-jkl)/l
               % We want only part of the integrand: beta*exp(-jkm)/m
               % and call the function EDintegratebetaoverm
                  
               [Iaccurate,singularterm] = EDintegratebetaoverm(envdata.cair,...
                   frequency,edgedata.closwedangvec(edge2),...
                   rS,thetaS,zS,re3_re2(jj),thetaout,ze3_re2(jj),...
                   [0,len2]);
                Iaccurate = -(2-acrossface_out)*ny/4/pi*Iaccurate;
                
                      
                % Here comes the LCN part, but with a single sample
                % replacement.
                % Iaccurate should be = len2*(w1*Gsub1 + w2*Gsub2 + ... +
                % wN*GsubN)
                % but due to the singularity, one term is imprecise: wL*GsubL             
                % Therefore we replace GsubL by L such that
                % w1*Gsub1 + ... + wL*L + ... + wN*GsubN = Iaccurate/len2 ->
                % L = (Iaccurate/len2 - dot(wvector,Gsubvector_mod))/wL
                % where Gsubvector_mod is [Gsub1,..., 0, ..., GsubN]
                % that is, with GsubL replaced by zero.

                Gsubvector_mod = Gsub(ivselect);
                Gsubvector_mod(ivpeaks) = 0;
                Lvalue = (Iaccurate/len2 - sum(Gsubvector_mod.*weight2vec))/weight2vec(ivpeaks);

                Gsuborig = Gsub(ivselect);
%                 reldiff = abs( Lvalue - Gsuborig(ivpeaks) )/abs( Gsuborig(ivpeaks) );
                    Gsub(ivselect(ivpeaks)) = Lvalue;
                Gsubmod = Gsub(ivselect);                

            end
            
          end

        end    % if accuratedone == 0
         
        G(Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) =     G(Hsubmatrixdata.Qvecstartnums(ii):Hsubmatrixdata.Qvecendnums(ii)) + Gsub;
                
    end
    
end

Qs = 1;
Qfirstterm = G*Qs;

if useaccurate && showtext >= 2
    disp(['      Source term, accurate: ',int2str(soudiagnostics.nslices_corrected),' of ',...
        int2str(soudiagnostics.nslices),' slices corrected, in ',int2str(soudiagnostics.npairs_corrected),...
        ' of ',int2str(soudiagnostics.npairs),' edge pairs'])
    disp(['         Pairs left to legacy: skewed ',int2str(soudiagnostics.npairs_skew),...
        ', partly visible ',int2str(soudiagnostics.npairs_partlyvisible),...
        ', nodes not GL ',int2str(soudiagnostics.npairs_nodesnotGL)])
end

% ========================================================================
% Local functions for the 'accurate' implementation. (No function in this
% file is closed with "end", so these are not either.)
% ========================================================================

function [Gsub,nslicescorrected,nodesok] = EDinteg_souterm_accurate(Gsub,ny,thetaS,thetaout,...
    rS,zS,re3_re2,ze3_re2,ze2_re2,weight2vec,len2,fvec,Gsubfactor)
% For each fixed point z3 on the to-edge (one "slice"), the source term
%   Gsub(z2) = Gsubfactor*beta(z2)*exp(-jkm)/m
% is a function of z2 along the via-edge. Near the apex point za of the
% path S -> edge2 -> z3, each beta term k behaves as
%   term_k(z2) ~ num_k/( (z2-za)^2 + B1_k )
% with za, B1, num as in EDcalcpropagatematrix (source S in the role of the
% receiver, z3 in the role of the far-edge node).
% For flagged (slice,term) combinations the nodal values are changed to
%   Gsub_eff(j) = Gsub(j) + Gsubfactor*boost*( W_j/dz_j - L_j )*f_j
% where L_j is the Lorentzian at node j, dz_j the GL weight (in metres),
% f_j = exp(-jkm_j)/m_j and W_j = integral of L(z)*l_j(z) dz (l_j = Lagrange
% polynomial of node j). A GL sum over z2 with Gsub_eff then equals the
% I_1 + I_2 split: the exact beta values integrated by plain GL, with the
% plain-GL estimate of the Lorentzian part replaced by its product-
% integration estimate.

n2 = length(ze2_re2);
n3 = length(ze3_re2);
nslicescorrected = 0;

[xstd,wstd,Sstd,tfine,wfine,Pint] = glStandardCached(n2);

% The modified weights assume ascending Gauss-Legendre nodes on the edge.
nodesok = double( max(abs(2*ze2_re2(:)/len2 - 1 - xstd)) < 1e-10 && max(abs(2*weight2vec(:) - wstd)) < 1e-10 );
if nodesok == 0
    return
end

if thetaout == 0
    phivec = pi + thetaS*[1 -1];
    termboost = 2;
else
    phivec = pi + thetaS*[1 1 -1 -1] + thetaout*[1 -1 1 -1];
    termboost = 1;
end

re3 = re3_re2(:);
ze3 = ze3_re2(:);
R0vec    = sqrt( (zS - ze3).^2 + (rS + re3).^2 );                  % [n3 x 1]
zavec    = ( ze3*rS + zS*re3 )./( rS + re3 );                      % [n3 x 1]
kappavec = ny^2*(rS + re3).^4 ./ (2*rS^2*re3.^2.*R0vec.^2);        % [n3 x 1]
B1_ik    = 2*sin(ny*phivec/2).^2 ./ kappavec;                      % [n3 x nk]
num_ik   = sin(ny*phivec) ./ kappavec;                             % [n3 x nk]

% Same flag criterion as in EDcalcpropagatematrix
localspacing = max(diff(sort(ze2_re2(:))));
zrelthreshold = 5;
insiderange  = (zavec >= -0.05*len2) & (zavec <= 1.05*len2);
narrowenough = sqrt(max(B1_ik,0)) < zrelthreshold*localspacing;
flag_ik = insiderange & narrowenough;                              % [n3 x nk]
nslicescorrected = nnz(any(flag_ik,2));
if nslicescorrected == 0
    return
end

Gmat  = reshape(Gsub,n2,n3);           % column jj = slice for to-edge point jj
z2col = ze2_re2(:);                    % [n2 x 1]
dzcol = weight2vec(:)*len2;            % [n2 x 1]
fcol  = fvec(:);                       % [n2 x 1]

for kk = 1:length(phivec)
    flaggedj = find(flag_ik(:,kk));
    if ~isempty(flaggedj)
        zaf  = reshape(zavec(flaggedj),1,[]);                % [1 x nf]
        B1f  = reshape(max(B1_ik(flaggedj,kk),0),1,[]);      % [1 x nf]
        numf = reshape(num_ik(flaggedj,kk),1,[]);            % [1 x nf]
        xa   = 2*zaf/len2 - 1;
        Bhat = 4*B1f/len2^2;
        WT = modifiedWeightsLorentzianVec(Sstd,wstd,tfine,wfine,Pint,xa,Bhat);   % [n2 x nf]
        modifiedest = (2/len2)*(numf.*WT)./dzcol;            % W_j/dz_j
        plainest    = numf./( (z2col - zaf).^2 + B1f );      % L_j
        Gmat(:,flaggedj) = Gmat(:,flaggedj) + Gsubfactor*termboost*(modifiedest - plainest).*fcol;
    end
end

Gsub = reshape(Gmat,n2*n3,1);

function [x,w] = gaussLegendreNodes(n)
% Golub-Welsch algorithm: standard n-point Gauss-Legendre nodes and
% weights on [-1,1], ascending.
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

function P = legendreVandermonde(x,Nmax)
% P(j,n+1) = P_n(x(j)), n=0..Nmax-1
x = x(:);
P = zeros(length(x),Nmax);
P(:,1) = 1;
if Nmax > 1
    P(:,2) = x;
end
for n = 2:Nmax-1
    P(:,n+1) = ((2*n-1)*x.*P(:,n) - (n-1)*P(:,n-1))/n;
end

function [x,w,S,tfine,wfine,Pint] = glStandardCached(n)
% Standard n-point Gauss-Legendre nodes x and weights w on [-1,1]
% (ascending), S(j,m+1) = (2m+1)*P_m(x_j), a finer 4n-point rule
% (tfine,wfine) and the matrix Pint [4n x n] which interpolates nodal
% values at x to the fine nodes. Cached, since they depend only on n.
persistent cachedn cachedx cachedw cachedS cachedt cachedwf cachedP
if isempty(cachedn)
    cachedn = [];  cachedx = {};  cachedw = {};  cachedS = {};
    cachedt = {};  cachedwf = {};  cachedP = {};
end
idx = find(cachedn == n,1);
if isempty(idx)
    [x,w] = gaussLegendreNodes(n);
    P = legendreVandermonde(x,n);
    S = P .* (2*(0:n-1)+1);
    [tfine,wfine] = gaussLegendreNodes(4*n);
    Pfine = legendreVandermonde(tfine,n);
    Pint = (Pfine*S.') .* (w(:).'/2);
    cachedn(end+1) = n;
    cachedx{end+1} = x;   cachedw{end+1} = w;   cachedS{end+1} = S;
    cachedt{end+1} = tfine;  cachedwf{end+1} = wfine;  cachedP{end+1} = Pint;
else
    x = cachedx{idx};  w = cachedw{idx};  S = cachedS{idx};
    tfine = cachedt{idx};  wfine = cachedwf{idx};  Pint = cachedP{idx};
end

function Q = neumannQupward(z,Nmax)
% Legendre functions of the second kind Q_0..Q_{Nmax-1} at complex z close
% to [-1,1], by upward recurrence. Returns [nz x Nmax].
zcol = z(:);
Q = zeros(length(zcol),Nmax);
Q(:,1) = 0.5*log((zcol+1)./(zcol-1));
if Nmax > 1
    Q(:,2) = zcol.*Q(:,1) - 1;
end
for n = 1:Nmax-2
    Q(:,n+2) = ( (2*n+1)*zcol.*Q(:,n+1) - n*Q(:,n) )/(n+1);
end

function WT = modifiedWeightsLorentzianVec(S,wnodes,tfine,wfine,Pint,za,B)
% Modified weights WT(j,q) = integral over [-1,1] of L_j(t)/((t-za_q)^2+B_q) dt
% for the standard GL nodes (L_j = Lagrange polynomial of node j).
% Narrow peaks (rho^(2N) < 1e4): Neumann's formula with Legendre functions
% of the second kind. Wide peaks: a 4N-point GL rule, which is then exact
% to rounding accuracy. Same routine as in EDcalcpropagatematrix.
za = za(:);
sB = max(sqrt(max(B(:),0)),1e-100);   % (a floor at eps would distort peaks narrower than 1e-8)
N = size(S,1);
WT = zeros(N,length(za));
p = za + 1i*sB;
rho = abs(p + sqrt(p.^2-1));
rho = max(rho,1./rho);
isnarrow = (2*N)*log(rho) < log(1e4);
iv = find(isnarrow);
if ~isempty(iv)
    Q = neumannQupward(p(iv),N);
    WT(:,iv) = -(S*imag(Q).') .* wnodes(:) ./ (sB(iv).');
end
iv = find(~isnarrow);
if ~isempty(iv)
    gfine = 1./( (tfine(:) - za(iv).').^2 + (sB(iv).').^2 );
    WT(:,iv) = Pint.'*(wfine(:).*gfine);
end
