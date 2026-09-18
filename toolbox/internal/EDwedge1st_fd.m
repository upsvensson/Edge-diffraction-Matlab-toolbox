function [tf,singularterm,zfirst,mlnearedge] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc,implementation)
% EDwedge1st_fd - Gives the 1st order diffraction transfer function.
% Gives the 1st order diffraction transfer function
% for a point source irradiating a finite wedge. A free-field tf amplitude
% of 1/r is used. An accurate integration method is used
% so receivers can be placed right at the zone boundaries.
%
% Output parameters:
%   tf              the transfer function
%   singularterm    a vector, [1 4], of 0 and 1 that indicates if one or more of the four
%                   beta terms is singular, i.e., right on a zone boundary.
%                   If there is one or more elements with the value 1, the
%                   direct sound or a specular reflection needs to be given
%                   half its specular value.
%   zfirst          A value between 0 and L, which tells which point along
%                   the edge gives the first arrival. A value of 0 or L
%                   implies that the apex point is outside the finite edge,
%                   and one of the edge end points causes the first arrival.
%                   A value between 0 and L is the value of the apex point
%                   position along the finite edge. zfirst could also be
%                   empty, in case the tf doesn't exist (for an edge of
%                   zero length). L is the length of the edge =
%                   abs(zw(2)-zw(1))
%   mlnearedge      1 if the source and/or receiver are close enough to
%                   the edge that the 1/(m*l) factor itself is sharply
%                   varying right around the apex point (in addition to,
%                   or independently of, any beta zone-boundary
%                   singularity). The analytical approximation below
%                   (tf2) only treats the beta-singularity and assumes
%                   1/(m*l) is ~constant over the tiny apex neighborhood
%                   zrangespecial; that assumption's accuracy is fine
%                   even in this case (see comment at its computation),
%                   but the flanking numerical integration (tf1, tf3)
%                   can become slow in the 'legacy' implementation, since
%                   it must resolve a sharp, frequency-independent bump
%                   in 1/(m*l) at every frequency. The 'fast'
%                   implementation (see below) grades its fixed mesh
%                   toward the source/receiver location specifically to
%                   address this. 0 otherwise.
%
% Input parameters:
%   cair            the speed of sound
%   frequencies	  	list of frequencies to compute the result for
%   closwedang	    closed wedge angle
%   rs, thetas, zs	cyl. coordinates of the source pos. (zs = 0)
%   rr, thetar, zr	cyl. coordinates of the receiver pos
%   zw		        the edge ends given as [z1 z2]
%   Method	        'New' (the new method by Svensson-Andersson-Vanderkooy) or 'Van'
%                   (Vanderkooys old method) or 'KDA' (Kirchhoff approximation).
%   Rstart (optional)	if this is included, the impulse response time zero will
%				    correspond to this distance. Otherwise, the time zero
%				    corresponds to the distance R0 via the apex point.
%   bc (optional) 	boundary condition, [1 1] (default) => rigid-rigid, [-1 -1] => soft-soft
%                   [1 -1] => rigid-soft. First plane should be the ref. plane.
%   implementation (optional) 'legacy' (default) uses the original,
%                   per-frequency adaptive-quadrature implementation,
%                   unchanged. 'fast' uses a fixed (frequency-independent)
%                   quadrature mesh, built once and reused across all
%                   frequencies:
%                     - the near-apex analytical term (tf2) is computed
%                       once outside the frequency loop instead of once
%                       per frequency, since it never depended on k in
%                       the first place;
%                     - the flanking numerical integrals (tf1, tf3) are
%                       replaced by a fixed composite Gauss-Legendre mesh
%                       (graded toward the source/receiver location, and
%                       sized to resolve oscillation at max(frequencies)),
%                       evaluated once per frequency via a single
%                       vectorized call + dot product instead of a fresh
%                       adaptive-quadrature solve.
%                   For the rigid boundary condition, tf1/tf3 in 'fast'
%                   evaluate the beta/(m*l)*exp(-i*k*(m+l-Rstart))
%                   integrand via the local function evalBetaOverMl,
%                   copied term-for-term from the formula the real
%                   EDquadgk.m inlines internally (it does not call a
%                   separate EDbetaoverml_fd for this calling
%                   convention) -- so 'fast' reproduces exactly what
%                   'legacy' integrates, not an assumed external
%                   function. The rigid apexincluded==0 and
%                   apexincluded==1 branches both use it; the Dirichlet
%                   (bc(1)==-1) branches still call
%                   EDbetaoverml_fd_dirichlet directly, unverified (see
%                   the comment in the 'legacy' Dirichlet branch --
%                   that path looks broken independently of this
%                   implementation switch). Compare 'legacy' and 'fast'
%                   on your own cases before relying on 'fast' in
%                   production -- the mesh-grading parameters inside
%                   buildGradedFlankMesh are heuristic starting points,
%                   not yet validated against a broad test set the way
%                   zrelmax was.
%
% Uses the built-in function QUADGK and the function EDbetaoverml_fd for numerical integration.
%
% Peter Svensson (peter.svensson@ntnu.no) 27 Aug. 2026
%
% [tf,singularterm,zfirst,mlnearedge] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc,implementation);

% 18 Jul. 2010   Stable version
% 19 Feb. 2013   Fixed the Dirichlet case
% 1 Nov. 2016    Implemented the analytic approximation for the case at the
%                zone boundaries, with Sara Martin
% 28 Nov. 2017 Copied to the EDtoolbox. Introduced the non-global input
% parameter cair.
% 15 Jan. 2018 Turned off the automatic printing out of "Singularity for
% term ..."
% 16 Jan. 2018 Implemented (de-commented) the serial expansion around the
% apex point.
% 17 Jan. 2018 Made some experiments with the zrelrange.
% 24 Jan 2018 Inserted a local switch for debugging: select analytical
% integration or not. Set the analytical integration to 0 for now.
% 25 Jan 2018 New implementation of the analytical integration, where only
% the cosh-factor is seen as a varying part of the integrand.
% 26 Jan 2018 Fixed bug: the case useserialexp2 had not been implemented.
% 12 Feb 2018 Introduced new output parameter: zfirst. This is a value
% between 0 and L, which tells which point along the edge gives the first
% arrival.
% 15 March 2021 Fixed bug: the code was tf2 = tf2/2; but should be
% tf2(jj) = tf2(jj)/2;
% 22 Nov. 2023 Changed the criterion for singularity
% 19 April 2024 Changed the numerical integration function from Matlab's
% built-in quadgk to a stripped-down version, EDquadgk.
% 27 Aug. 2026 Code by Claude, prompted by PS
% - Added output mlnearedge: flags (but does not yet specially
% handle in 'legacy') geometries where the source/receiver are close
% enough to the edge that 1/(m*l), not just beta, varies sharply near
% the apex point.
% - Added input implementation ('legacy'/'fast') and the 'fast'
% code path: tf2 hoisted out of the frequency loop (frequency-independent),
% and tf1/tf3 evaluated on a fixed, geometry- and frequency-graded
% Gauss-Legendre mesh built once, reused across all frequencies.
% - After seeing the real EDquadgk.m: it inlines the rigid-BC
% integrand itself rather than calling a separate EDbetaoverml_fd, so
% 'fast' now inlines the same exact formula (evalBetaOverMl) instead of
% calling an assumed external function. Also noted, but not changed: the
% legacy apexincluded==0/bc(1)==-1 branch passes EDquadgk a function
% handle, but the real EDquadgk.m always expects a struct -- that branch
% looks broken independently of this implementation switch.

localshowtext = 0;

localselectanalytical = 1;

% doublecheck = 1;

% if localshowtext >= 2
%     disp(' ')
%     disp('Reference solution for wedge 1st')
% end

if nargin < 12
	Rstart = 0;
end
if nargin < 13
	bc = [1 1];
else
	if bc(1)*bc(2)~=1
		error('ERROR: Only all-rigid, or all-soft, wedges are implemented');
	end
end
% % implementation = 'fast';
% % disp('WARNING!! implementation forced to fast inside EDwedge1st_fd')
if nargin < 14 || isempty(implementation)
    implementation = 'legacy';
end
if ~(strcmpi(implementation,'legacy') || strcmpi(implementation,'fast'))
    error(['ERROR: implementation must be ''legacy'' or ''fast'', got ''',implementation,'''']);
end


%------------------------------------------------------------------
% If the wedge has the length zero, return a zero TF immediately

if zw(2)-zw(1) == 0
	tf = zeros(size(frequencies));
    singularterm = [0 0 0 0];
    zfirst = [];
    mlnearedge = 0;
	return
end

%----------------------------------------------
% Method + Some edge constants

ny = pi/(2*pi-closwedang);

if Method(1) == 'N' || Method(1) == 'n'
	Method = 'New';
else
	error(['ERROR: The method ',Method,' has not been implemented yet!']);
end

nfrequencies = length(frequencies);

tol = 1e-11;                        % The relative accuracy for the numerical integration
                                    % It was found that 1e-6 generated a
                                    % substantial error for some cases (PC
                                    % 041129)

zrelmax = 0.001;                    % The part of the edge that should use the analyt. approx.
                                    % Some simple tests indicated that
                                    % 0.001 => rel.errors smaller than 4e-7
                                    % 0.005 => rel. errors smaller than 3e-6
                                    %

za = (rs*zr+rr*zs)/(rs+rr);         % The apex point: the point on the edge with the shortest path

% Move the z-values so that za ends up in z = 0.

% Dz = abs(zs-zr);
zs = zs - za;
zr = zr - za;
zw = zw - za;                       % zw is a two-element vector
% za = 0;
% disp(['Positive edge-end z-value is ',num2str(zw(2))])
if localshowtext == 1
    disp(['            Shifted edge extension is ',num2str(zw(1)),' to ',num2str(zw(2))])
end

% rs2 = rs^2;
% rr2 = rr^2;

R0 = sqrt( (rs+rr)^2 + (zr-zs)^2 );
m0 = sqrt(rs^2 + zs^2);
l0 = sqrt(rr^2 + zr^2);

% ---------------------------------------------------------------------
% Flag geometries where the source and/or receiver sit close enough to
% the edge that 1/(m*l) itself is sharply varying right around the apex
% point, as distinct from the beta-factor zone-boundary singularity
% handled below.
%
% Below (tf2), 1/(m*l) is approximated by its constant value 1/(m0*l0)
% over the small apex-neighbourhood zrangespecial = zrelmax*min(rs,rr).
% That approximation's relative error is governed by zrangespecial/m0
% and zrangespecial/l0, and since m0 >= rs, l0 >= rr, and zrangespecial
% = zrelmax*min(rs,rr), both ratios are bounded by zrelmax regardless of
% how small rs or rr get -- so tf2 itself stays accurate (this is the
% same conclusion as the beta-only Lorentzian approximation derived
% separately: beta_i(z) ~ B1_i*cot(ny*phi_i/2)/((z-za)^2+B1_i), with
% 1/(ml) folded in as a constant because it varies negligibly over
% zrangespecial).
%
% The actual cost falls on the flanking numerical integration (tf1,
% tf3): if rs or rr is small, m(z) or l(z) has a sharp, but
% frequency-independent, near-singular bump of width ~min(rs,rr) that
% sits just outside zrangespecial (since zrangespecial is 1/zrelmax =
% 1000x narrower than that bump). In the 'legacy' implementation,
% EDquadgk must adaptively resolve this bump -- and because it is called
% inside the frequency loop, that adaptive work is repeated at every
% single frequency even though the bump's shape and location never
% change with frequency. The 'fast' implementation's mesh grades toward
% zs/zr specifically to handle this without repeating the adaptive work.

mlsharpthreshold = 0.05;             % fraction of R0; tune empirically,
                                      % analogous to zrelmax above
mlnearedge = min(m0,l0) < mlsharpthreshold*R0;
% % % if mlnearedge
% % %     warning('EDwedge1st_fd:mlNearEdge', ...
% % %         ['Source and/or receiver are close to the edge (min(m0,l0) = ', ...
% % %         num2str(min(m0,l0)),', R0 = ',num2str(R0), ...
% % %         '). The 1/(m*l) factor is sharply varying near the apex; ', ...
% % %         'results are still accurate but the ''legacy'' numerical ', ...
% % %         'integration may be slow across many frequencies. See ', ...
% % %         'mlnearedge output and consider implementation = ''fast''.']);
% % % end

% Calculate the sinnyfi and cosnyfi values for use in the four beta-terms
% later

fivec = pi + thetas*[1 1 -1 -1] + thetar*[1 -1 1 -1];

absnyfivec = abs(ny*fivec);

sinnyfivec = sin(ny*fivec);
cosnyfivec = cos(ny*fivec);

%------------------------------------------------------------------
% Check which category we have.
%
% Is the apex point inside the physical wedge or not? We want to use an
% analytic approximation for a small part of the apex section (and
% an ordinary numerical integration for the rest of the wedge).
%
%       apexincluded        0 or 1

apexincluded = 1;
if sign( zw(1)*zw(2) ) == 1
    apexincluded = 0;
end
if localshowtext == 1
    disp(['            apexincluded = ',int2str(apexincluded)])
end

%------------------------------------------------------------------
%

if apexincluded == 1

    zfirst = za;

    % Finally the endpoint of the integral of the analytical approximation
	% (=zrangespecial) should be only up to z = 0.01
	% or whatever is specified as zrelmax for the analytic approximation

	zrangespecial = zrelmax*min([rs,rr]);

    if zw(1) > 0
       error('ERROR! Problem: zw(1) > 0')
    end

    zanalyticalstart = zrangespecial;
    zanalyticalend = zrangespecial;
    if zw(1) < -zrangespecial
        includepart1 = 1;
    else
        includepart1 = 0;
        zanalyticalstart = -zw(1); % NB! zanalyticalstart gets a pos.value
    end
    if zw(2) > zrangespecial
        includepart2 = 1;
    else
        includepart2 = 0;
        zanalyticalend = zw(2);
    end
    if zw(1) == 0 || zw(2) == 0
       exacthalf = 1;
    else
        exacthalf = 0;
    end
    if localshowtext == 1
        disp(['            zrangespecial = ',num2str(zrangespecial)])
        disp(['            exacthalf = ',int2str(exacthalf)])
    end
    if includepart1*includepart2 == 0 && exacthalf == 0
       analyticalsymmetry = 0;
    else
       analyticalsymmetry = 1;
    end
    if localshowtext == 1
        disp(['            analyticalsymmetry = ',num2str(analyticalsymmetry)])
        disp(['             '])
    end

end

%----------------------------------------------------------------
% Compute the transfer function by carrying out a numerical integration for
% the entire wedge extension.

% singularterm = [0 0 0 0];
% 22 Nov. 2023 Changed the criterion for singularity by a factor of 1e3
singularterm = absnyfivec < 1e4*eps | abs(absnyfivec - 2*pi) < 1e4*eps;
useterm = 1 - singularterm;
if any(singularterm) && localshowtext
     disp(['            Singularity for term ',int2str(find(singularterm))])
end

tf = zeros(nfrequencies,1);
integralparameters = struct;
integralparameters.rs = rs;
integralparameters.rr = rr;
integralparameters.zs = zs;
integralparameters.zr = zr;
integralparameters.ny = ny;
integralparameters.sinnyfivec = sinnyfivec;
integralparameters.cosnyfivec = cosnyfivec;
integralparameters.useterm = useterm;


Rextra = R0 - Rstart;
Rstarttemp = R0;
integralparameters.Rstart = Rstarttemp;

if strcmpi(implementation,'legacy')

%==================================================================
% 'legacy': original, unmodified per-frequency implementation
%==================================================================

if apexincluded == 0
    if bc(1) == 1
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            integralparameters.k = k;
%              tf(ii) = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),'RelTol',tol)*(-ny/4/pi);
           tf(ii) = EDquadgk(integralparameters,zw(1),zw(2),tol)*(-ny/4/pi);
            tf(ii) = tf(ii)*exp(-1i*k*Rextra);
        end
    else
        % [observation, not changed] The real EDquadgk.m always reads its
        % first argument as integralparameters.<field> (rs, rr, k, ...);
        % it never checks isstruct or accepts a function handle. The
        % call below passes a function handle, not a struct, so this
        % branch looks like it would error if actually exercised. Left
        % untouched here since 'legacy' must stay byte-identical to the
        % original for a fair baseline -- flagging it for you to check
        % independently of anything in this implementation switch.
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            integralparameters.k = k;
%              tf(ii) = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),'RelTol',tol)*(-ny/4/pi);
           tf(ii) = EDquadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),tol)*(-ny/4/pi);
            tf(ii) = tf(ii)*exp(-1i*k*Rextra);
        end
    end

    Rend1  = sqrt( rs^2 + (zw(1)-zs)^2 ) + sqrt( rr^2 + (zw(1)-zr)^2 );
    Rend2  = sqrt( rs^2 + (zw(2)-zs)^2 ) + sqrt( rr^2 + (zw(2)-zr)^2 );

    if Rend1 < Rend2
        zfirst = 0;
    else
       zfirst = zw(2)-zw(1);
    end

else  % if apexincluded == 0

    if bc(1) == 1
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            integralparameters.k = k;
            if localselectanalytical == 0
%                 tf(ii) = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),'RelTol',tol)*(-ny/4/pi)*exp(-1i*k*Rextra);
               tf(ii) = EDquadgk(integralparameters,zw(1),zw(2),tol)*(-ny/4/pi)*exp(-1i*k*Rextra);

            else
                if includepart1 == 1
%                     tf1 = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),-zrangespecial,'RelTol',tol);
                    tf1 = EDquadgk(integralparameters,zw(1),-zrangespecial,tol);
                    tf1 = tf1*(-ny/4/pi);
                    if localshowtext
                       disp(['            Part1 included, from ',num2str(zw(1)),' to ',num2str(-zrangespecial)]);
                       disp(['            tf1 = ',num2str(tf1)])
                    end
                else
                   tf1 = 0;
                   if localshowtext
                       disp('            Part1 not included');
                   end
                end
                if includepart2 == 1
%                     tf3 = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zrangespecial,zw(2),'RelTol',tol);
                   tf3 = EDquadgk(integralparameters,zrangespecial,zw(2),tol);
                    tf3 = tf3*(-ny/4/pi);
                    if localshowtext
                       disp(['            Part2 included, from ',num2str(zrangespecial),' to ',num2str(zw(2))]);
                       disp(['            tf3 = ',num2str(tf3)])
                   end
                else
                   tf3 = 0;
                   if localshowtext
                       disp('            Part2 not included');
                   end
                end

                %-----------------------------------------------------------
                % The analytical approximation

                useserialexp1 = absnyfivec < 0.01;
                useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
                useserialexp = (useserialexp1 | useserialexp2) & (~singularterm);
                if localshowtext==1 && any(useserialexp)==1
                       disp('            Using series expansion for the cosnyfivec');
                end

                if ~any(singularterm)
                    cosnyfifactor = ny./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...
                                1./abs(fivec).*(useserialexp1==1) + ...
                                1./sqrt( (fivec-2*pi/ny).^2 ).*(useserialexp2==1);
                    sinovercosnyfifactor =  sinnyfivec./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...
                                sign(fivec).*(1- (ny*fivec).^2/6).*(useserialexp1==1) + ...
                                sign(ny*fivec-2*pi) .*(1- (ny*fivec-2*pi).^2/6).*(useserialexp2==1);
                    if analyticalsymmetry == 1
                        tf2 = -sinovercosnyfifactor/pi/R0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                        if exacthalf == 1
                            tf2 = tf2/2;
                        end
                        if localshowtext == 1
                            disp(['            Analytical approx from -zrangespecial to +zrangespecial'])
                            if exacthalf == 1
                                disp(['               with exacthalf = 1 -> half the amplitude'])
                            end
                        end
                    else
                        tf2 = -sinovercosnyfifactor/pi/R0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;
                        if localshowtext == 1
                            disp(['            Analytical approx from ',num2str(-zanalyticalstart),' to ',num2str(zanalyticalend)])
                        end
                    end
                    if localshowtext == 1
                        disp(['            tf2 = ',num2str(tf2)])
                    end
                else   % ~any(singularterm)
                  if localshowtext
                       disp('            Some term is zero due to singularity -> for-loop over the 4 terms');
                   end
                   tf2 = zeros(1,4);
                   for jj = 1:4
                       if ~singularterm(jj)
                            cosnyfifactor = ny./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...
                                        1./abs(fivec(jj)).*(useserialexp1(jj)==1);
                            sinovercosnyfifactor =  sinnyfivec(jj)./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...
                                        sign(fivec(jj)).*(1- (ny*fivec(jj)).^2/6).*(useserialexp1(jj)==1);
                            if analyticalsymmetry == 1
                                tf2(jj) = -sinovercosnyfifactor/pi/R0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                                if exacthalf == 1
                                    % Fixed a bug 15.3.2021; the old wrong
                                    % version was the line below, which changed
                                    % all 4 beta-values
                                    % tf2 = tf2/2;
                                    tf2(jj) = tf2(jj)/2;
                                end
                                if localshowtext == 1
                                    disp(['            Term ',int2str(jj),': Analytical approx from -zrangespecial to +zrangespecial'])
                                    if exacthalf == 1
                                        disp(['               with exacthalf = 1 -> half the amplitude'])
                                    end
                                end

                            else  % if analyticalsymmetry == 1
                                tf2(jj) = -sinovercosnyfifactor/pi/R0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;
                                if localshowtext == 1
                                    disp(['            Term ',int2str(jj),': Analytical approx from ',num2str(-zanalyticalstart),' to ',num2str(zanalyticalend)])
                                end

                            end   % if analyticalsymmetry == 1
                       else
                            if localshowtext == 1
                                disp(['            Term ',int2str(jj),': zero due to singularity'])
                            end
                       end   % if ~singularterm(jj)
                       if localshowtext == 1
                            disp(['            tf2(',int2str(jj),') = ',num2str(tf2(jj))])
                       end

                   end  % for jj = 1:4
                end   % if ~any(singularterm)

                tf2 = sum(tf2);

                tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*Rextra);

            end
        end
    else   % if bc(1) == 1
        error('ERROR: Not implemented the integration properly for the Dirichlet wedge yet');
%         for ii = 1:nfrequencies
%             k = 2*pi*frequencies(ii)/cair;
%             tf1 = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),-zrangespecial,'RelTol',tol)*(-ny/4/pi);
%             tf2 = 0;
%             tf3 = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zrangespecial,zw(2),'RelTol',tol)*(-ny/4/pi);
%             tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*Rextra);
%         end
    end
end

else

%==================================================================
% 'fast': tf2 hoisted out of the frequency loop (tier 1); tf1/tf3
% evaluated on a fixed, geometry- and frequency-graded mesh reused
% across all frequencies (tier 2).
%==================================================================

    kmax = 2*pi*max(frequencies)/cair;   % worst case for oscillation resolution;
                                          % a mesh that resolves this resolves
                                          % every lower frequency too.

    meshopts = struct;
    meshopts.Ngl = 12;             % Gauss-Legendre points per panel
    meshopts.phaselimit = pi/6;    % max phase change of exp(-i*k*R) allowed per panel at kmax
    meshopts.gradefactor = 0.05;   % width of the smallest panel next to a
                                    % feature point, as a fraction of that
                                    % feature's length scale (rs or rr)
    meshopts.growfactor = 1.7;     % panel-to-panel growth factor when
                                    % grading away from a feature point
    meshopts.maxpanels = 60;       % safety cap on panels per graded run
    % These are heuristic starting points, in the same spirit as
    % zrelmax above -- tune against your own accuracy/speed tests.

    if apexincluded == 0

        featurepts = [zs, zr];
        featurescales = [rs, rr];
        [zn,wn] = buildGradedFlankMesh(zw(1),zw(2),featurepts,featurescales,kmax,meshopts);

        if bc(1) == 1
            for ii = 1:nfrequencies
                k = 2*pi*frequencies(ii)/cair;
                fvals = evalBetaOverMl(zn,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm);
                tf(ii) = sum(wn.*fvals(:))*(-ny/4/pi);
                tf(ii) = tf(ii)*exp(-1i*k*Rextra);
            end
        else
            for ii = 1:nfrequencies
                k = 2*pi*frequencies(ii)/cair;
                fvals = EDbetaoverml_fd_dirichlet(zn,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm);
                tf(ii) = sum(wn.*fvals(:))*(-ny/4/pi);
                tf(ii) = tf(ii)*exp(-1i*k*Rextra);
            end
        end

        Rend1  = sqrt( rs^2 + (zw(1)-zs)^2 ) + sqrt( rr^2 + (zw(1)-zr)^2 );
        Rend2  = sqrt( rs^2 + (zw(2)-zs)^2 ) + sqrt( rr^2 + (zw(2)-zr)^2 );

        if Rend1 < Rend2
            zfirst = 0;
        else
           zfirst = zw(2)-zw(1);
        end

    else  % if apexincluded == 0

        if bc(1) == 1

            % ---- tier 1: tf2 computed once, identical formula to legacy ----

            useserialexp1 = absnyfivec < 0.01;
            useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
            useserialexp = (useserialexp1 | useserialexp2) & (~singularterm);

            if ~any(singularterm)
                cosnyfifactor = ny./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...
                            1./abs(fivec).*(useserialexp1==1) + ...
                            1./sqrt( (fivec-2*pi/ny).^2 ).*(useserialexp2==1);
                sinovercosnyfifactor =  sinnyfivec./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...
                            sign(fivec).*(1- (ny*fivec).^2/6).*(useserialexp1==1) + ...
                            sign(ny*fivec-2*pi) .*(1- (ny*fivec-2*pi).^2/6).*(useserialexp2==1);
                if analyticalsymmetry == 1
                    tf2 = -sinovercosnyfifactor/pi/R0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                    if exacthalf == 1
                        tf2 = tf2/2;
                    end
                else
                    tf2 = -sinovercosnyfifactor/pi/R0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;
                end
            else
               tf2 = zeros(1,4);
               for jj = 1:4
                   if ~singularterm(jj)
                        cosnyfifactor = ny./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...
                                    1./abs(fivec(jj)).*(useserialexp1(jj)==1);
                        sinovercosnyfifactor =  sinnyfivec(jj)./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...
                                    sign(fivec(jj)).*(1- (ny*fivec(jj)).^2/6).*(useserialexp1(jj)==1);
                        if analyticalsymmetry == 1
                            tf2(jj) = -sinovercosnyfifactor/pi/R0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                            if exacthalf == 1
                                tf2(jj) = tf2(jj)/2;
                            end
                        else
                            tf2(jj) = -sinovercosnyfifactor/pi/R0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;
                        end
                   end   % if ~singularterm(jj)
               end  % for jj = 1:4
            end   % if ~any(singularterm)
            tf2 = sum(tf2);

            % ---- tier 2: fixed meshes for tf1, tf3, built once ----

            featurepts = [zs, zr];
            featurescales = [rs, rr];

            if includepart1 == 1
                [zn1,wn1] = buildGradedFlankMesh(zw(1),-zrangespecial,featurepts,featurescales,kmax,meshopts);
            else
                zn1 = []; wn1 = [];
            end
            if includepart2 == 1
                [zn3,wn3] = buildGradedFlankMesh(zrangespecial,zw(2),featurepts,featurescales,kmax,meshopts);
            else
                zn3 = []; wn3 = [];
            end

            for ii = 1:nfrequencies
                k = 2*pi*frequencies(ii)/cair;

                if ~isempty(zn1)
                    fvals1 = evalBetaOverMl(zn1,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm);
                    tf1 = sum(wn1.*fvals1(:))*(-ny/4/pi);
                else
                    tf1 = 0;
                end
                if ~isempty(zn3)
                    fvals3 = evalBetaOverMl(zn3,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm);
                    tf3 = sum(wn3.*fvals3(:))*(-ny/4/pi);
                else
                    tf3 = 0;
                end

                tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*Rextra);
            end

        else   % if bc(1) == 1
            error('ERROR: Not implemented the integration properly for the Dirichlet wedge yet');
        end
    end

end   % if strcmpi(implementation,'legacy')

end   % main function


% =========================================================================
function [zn,wn] = buildGradedFlankMesh(zleft,zright,featurepts,featurescales,kmax,opts)
% Build a fixed composite Gauss-Legendre mesh on [zleft,zright], graded
% geometrically toward any feature point (e.g. the source or receiver
% z-projection, zs or zr) that falls strictly inside the interval, and
% otherwise sized so that kmax*panelwidth <= opts.phaselimit -- this
% keeps oscillation at the highest frequency in play well resolved by an
% opts.Ngl-point Gauss rule per panel. The mesh does not depend on k
% except through kmax, so it can be built once and reused for every
% frequency <= max(frequencies).
%
% featurepts / featurescales: equal-length vectors of candidate grading
% points and their local length scales; a feature point outside
% (zleft,zright) is ignored.

    if zright <= zleft
        zn = []; wn = [];
        return
    end

    tol = 1e-9*max(1,abs(zright-zleft));

    keep = featurepts > zleft+tol & featurepts < zright-tol;
    featurepts = featurepts(keep);
    featurescales = featurescales(keep);

    breakpts = unique([zleft, featurepts(:).', zright]);

    [xg,wg] = gaussLegendreNodes(opts.Ngl);   % fixed nodes/weights on [-1,1]

    zn = [];
    wn = [];

    for bb = 1:numel(breakpts)-1
        a = breakpts(bb);
        b = breakpts(bb+1);
        if b <= a
            continue
        end

        sa = featurescales(abs(featurepts-a) < tol);
        sb = featurescales(abs(featurepts-b) < tol);
        wmax = min(b-a, opts.phaselimit/max(kmax,eps));

        if ~isempty(sa) || ~isempty(sb)
            if ~isempty(sa) && ~isempty(sb)
                mid = (a+b)/2;
                panels = [gradedPanels(a,mid,sa(1),wmax,opts,'up'); ...
                          gradedPanels(mid,b,sb(1),wmax,opts,'down')];
            elseif ~isempty(sa)
                panels = gradedPanels(a,b,sa(1),wmax,opts,'up');
            else
                panels = gradedPanels(a,b,sb(1),wmax,opts,'down');
            end
        else
            npan = max(1,ceil((b-a)/wmax));
            edges = linspace(a,b,npan+1);
            panels = [edges(1:end-1)', edges(2:end)'];
        end

        for pp = 1:size(panels,1)
            pa = panels(pp,1); pb = panels(pp,2);
            if pb <= pa
                continue
            end
            zn = [zn; (pb-pa)/2*xg + (pa+pb)/2]; %#ok<AGROW>
            wn = [wn; (pb-pa)/2*wg]; %#ok<AGROW>
        end
    end
end


% =========================================================================
function panels = gradedPanels(a,b,scale,wmax,opts,direction)
% Geometrically graded panel edges on [a,b]. direction='up' places the
% smallest panel next to a, growing toward b; direction='down' places
% the smallest panel next to b, growing toward a.
%
% Two phases, so panel count is never capped by opts.maxpanels while
% oscillation resolution is still outstanding:
%   1. Ramp phase: geometric growth from opts.gradefactor*scale by
%      opts.growfactor per panel, only while that width is still finer
%      than the phase-limited cap wmax (i.e. only while grading toward
%      the feature is actually buying extra resolution beyond what
%      oscillation already demands everywhere). Bounded by
%      opts.maxpanels as a sanity cap on this phase only.
%   2. Uniform phase: the remainder of [a,b] is covered by panels of
%      width exactly wmax (rounded down to fit evenly) -- as many as
%      needed, unconditionally, so oscillation resolution is never
%      silently sacrificed to a panel-count cap.

    L = b-a;
    w0 = min(max(opts.gradefactor*scale, L*1e-6), wmax);
    edges = 0;
    w = w0;
    nramp = 0;
    while w < wmax - 1e-12*max(1,wmax) && edges(end) < L - 1e-9*max(1,L) && nramp < opts.maxpanels
        remaining = L - edges(end);
        wstep = min([w, remaining]);
        edges(end+1) = edges(end) + wstep; %#ok<AGROW>
        w = w*opts.growfactor;
        nramp = nramp+1;
    end

    remaining = L - edges(end);
    if remaining > 1e-9*max(1,L)
        npan = max(1,ceil(remaining/wmax));
        step = remaining/npan;
        edges = [edges, edges(end)+step*(1:npan)]; %#ok<AGROW>
    end
    edges(end) = L;   % clean up any round-off

    if strcmp(direction,'up')
        e = a + edges;
    else
        e = b - fliplr(edges);
    end
    e = sort(e);
    panels = [e(1:end-1)', e(2:end)'];
end


% =========================================================================
function fval = evalBetaOverMl(zvec,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm)
% Rigid-wedge beta/(m*l)*exp(-i*k*(m+l-Rstart)) integrand, vectorized
% over zvec. Written to match, term for term, the formula the real
% EDquadgk.m inlines internally for efficiency (it does not call a
% separate EDbetaoverml_fd function for the struct/rigid-BC calling
% convention used here) -- so 'fast' reproduces exactly what 'legacy'
% integrates, rather than relying on an assumed external function with
% an unverified signature.

    zvec = zvec(:);
    if ~any(useterm)
        fval = zeros(size(zvec));
        return
    end
    m = sqrt( (zvec-zs).^2 + rs^2 );
    l = sqrt( (zvec-zr).^2 + rr^2 );
    ml = m.*l;
    y = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
    expnyeta = real((real(sqrt(y.^2-1)) + y).^ny);
    coshnyeta = 0.5*(expnyeta + 1./expnyeta);
    fval = zeros(size(zvec));
    for ii = 1:4
        if useterm(ii)
            fval = fval + sinnyfivec(ii)./(coshnyeta - cosnyfivec(ii));
        end
    end
    fval = fval./ml .* exp(-1i*k*(m+l-Rstart));
end


% =========================================================================
function [x,w] = gaussLegendreNodes(n)
% n-point Gauss-Legendre nodes/weights on [-1,1] via the Golub-Welsch
% algorithm (eigenvalues/eigenvectors of the Jacobi matrix).
    beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(T);
    [x,idx] = sort(diag(D));
    w = 2*(V(1,idx).^2)';
    x = x(:);
    w = w(:);
end
