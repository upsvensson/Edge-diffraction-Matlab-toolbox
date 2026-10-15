function [tf,singularterm,zfirst] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc,diff1method,diff1GLorder)
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
%   diff1method     (optional) 'accurate' or 'GL'. Default = 'accurate'.
%                   With 'GL', the 'accurate' method is still used for
%                   receivers close to a zone boundary (when the apex
%                   point is on the edge, and the peak of the integrand
%                   there is narrower than 5 GL node spacings).
%   diff1GLorder    (optional) GL order; only relevant if diff1method =
%                   'GL'.
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [tf,singularterm,zfirst] = EDwedge1st_fd(cair,frequencies,closwedang,rs,...
%    thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc,diff1method,diff1GLorder);

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
% 24 Sep. 2026 Introduced GL quadrature as an option, controlled by the
% optional input parameters diff1method and diff1GLorder.
% 14 Oct. 2026 (by Claude, prompted by PS) The defaults of the optional
% input parameters work also when fewer than 13 parameters are given.
% For diff1method = 'GL', a case which is singular (a receiver exactly on
% a zone boundary) is computed with the 'accurate' method instead; before,
% the check for this was in the 'accurate' branch, so it was never used.
% An unknown diff1method gives an error.
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).
% 15 Oct. 2026 (by Claude, prompted by PS) For diff1method = 'GL', also
% near-singular cases are computed with the 'accurate' method: when the
% apex point is on the edge and the peak of the integrand at the apex
% point is narrower than 5 GL node spacings. Before, GL gave errors of up
% to 100% for receivers within a few degrees of a zone boundary.

localshowtext = 0;

localselectanalytical = 1;

% doublecheck = 1;

% if localshowtext >= 2
%     disp(' ')
%     disp('Reference solution for wedge 1st')
% end

if nargin < 12 || isempty(Rstart)
	Rstart = 0;
end
if nargin < 13 || isempty(bc)
	bc = [1 1];
end
if bc(1)*bc(2)~=1
	error('ERROR: Only all-rigid, or all-soft, wedges are implemented');
end
if nargin < 14 || isempty(diff1method)
    diff1method = 'accurate';
end
if nargin < 15 || isempty(diff1GLorder)
    if strcmp(diff1method,'GL')
        diff1GLorder = 16;
    else
        diff1GLorder = NaN;
    end
end
if ~strcmp(diff1method,'accurate') && ~strcmp(diff1method,'GL')
    error(['ERROR: diff1method must be ''accurate'' or ''GL'', not ''',diff1method,''''])
end

%------------------------------------------------------------------
% If the wedge has the length zero, return a zero TF immediately

if zw(2)-zw(1) == 0
	tf = zeros(size(frequencies));
    singularterm = [0 0 0 0];
    zfirst = [];
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

% Calculate the sinnyfi and cosnyfi values for use in the four beta-terms
% later

fivec = pi + thetas*[1 1 -1 -1] + thetar*[1 -1 1 -1];

absnyfivec = abs(ny*fivec);

sinnyfivec = sin(ny*fivec);
cosnyfivec = cos(ny*fivec);

if strcmp(diff1method,'accurate')
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
    %
    % When there is a truly singular case, very, very close to a zone
    % boundary, then tfdiff should be zero (for 1-4 of the terms).
    
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
elseif strcmp(diff1method,'GL')
%    [nvec,wvec] = EDdistelements(diff1GLorder,2);
    [nvec,wvec] = EDgetGLrule(diff1GLorder);
    zwlength = zw(2) - zw(1);
    zvec = nvec*zwlength + zw(1);

    % Close to a zone boundary, the integrand has a narrow peak at the
    % apex point (z = 0 after the shift above). As in the analytical
    % approximation of the 'accurate' branch, the integral over the peak
    % behaves as atan(z/wpeak), with the half-width
    %   wpeak = m0*l0*sqrt(2*(1-cos(ny*fi)))/(ny*R0)
    % for each of the four terms. GL quadrature can not resolve a peak
    % that is narrower than a few node spacings, so if the apex point is
    % on the edge (with a 5% margin) and wpeak is smaller than 5 node
    % spacings near the apex point, the 'accurate' method is used instead.
    % This includes the singular case (a receiver exactly on a zone
    % boundary), where wpeak = 0. Away from the zone boundaries, the
    % relative error of the GL result is then below about 1e-8.
    wpeak = m0*l0*sqrt(2*(1-cosnyfivec))/(ny*R0);
    apexonedge = zw(1) <= 0.05*zwlength && zw(2) >= -0.05*zwlength;
    if apexonedge
        zsorted = sort([zw(1);zvec(:);zw(2)]);
        zapexclamped = min(max(0,zw(1)),zw(2));
        kinterval = find(zsorted <= zapexclamped,1,'last');
        kinterval = min(kinterval,numel(zsorted)-1);
        nodespacings = diff(zsorted);
        hnearapex = max(nodespacings(max(1,kinterval-1):min(numel(nodespacings),kinterval+1)));
        if any(wpeak < 5*hnearapex) || any(absnyfivec < 1e4*eps | abs(absnyfivec - 2*pi) < 1e4*eps)
            [tf,singularterm,zfirst] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs+za,rr,thetar,zr+za,zw+za,Method,Rstart,bc,'accurate',NaN);
            return
        end
    end
    zfirst = NaN;
    singularterm = [0 0 0 0];

    tf = zeros(nfrequencies,1);
    kvec = 2*pi*frequencies/cair;
    useterm = [1 1 1 1];

    for ii = 1:nfrequencies
        outval = EDbetaoverml_fd(zvec,kvec(ii),rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstart,useterm);
        tf(ii) = sum(outval(:).*wvec(:));
    end
    tf = tf*(-ny/4/pi);
end

