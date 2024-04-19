function [tf,singularterm,zfirst] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc)
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
%
% Uses the built-in function QUADGK and the function EDbetaoverml_fd for numerical integration.
%
% Peter Svensson (peter.svensson@ntnu.no) 19 April 2024
%
% [tf,singularterm] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc);

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


