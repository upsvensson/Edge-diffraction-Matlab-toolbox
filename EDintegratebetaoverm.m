function [tf,singularterm] = EDintegratebetaoverm(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw)
% EDintegratebetaoverm - Integrates the function beta*exp(-jkm)/m.
% This is similar to first-order diffraction but lacks the exp(-jkl)/l factor.
% An accurate integration method is used
% so receivers can be placed right at the zone boundaries.
%
% Output parameters:
%   tf              the transfer function
%   singularterm    a vector, [1 4], of 0 and 1 that indicates if one or more of the four
%                   beta terms is singular, i.e., right on a zone boundary.
%                   If there is one or more elements with the value 1, the
%                   direct sound or a specular reflection needs to be given
%                   half its specular value.
% 
% Input parameters:
%   cair            the speed of sound
%   frequencies	  	list of frequencies to compute the result for
%   closwedang	    closed wedge angle
%   rs, thetas, zs	cyl. coordinates of the source pos. (zs = 0)
%   rr, thetar, zr	cyl. coordinates of the receiver pos
%   zw		        the edge ends given as [z1 z2]
%
% Uses the built-in function QUADGK and the function EDbetaoverm_fd for numerical integration.
%
% Peter Svensson (peter.svensson@ntnu.no) 28 May 2018
%
% [tf,singularterm] = EDintegratebetaoverm(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw);

% 26 May 2018 Modified the EDwedge1st_fd by removing some parts. Method,
% Rstart and bc are not used here.
% 28 May 2018 Did some corrections, switching R0 to m0

localshowtext = 0;

localselectanalytical = 1;

% doublecheck = 0;

% if localshowtext >= 2
%     disp(' ')
%     disp('Reference solution for wedge 1st')
% end

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

nfrequencies = length(frequencies);

tol = 1e-11;                         % The relative accuracy for the numerical integration
                                     % The relative accuracy for the numerical integration
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

% if localshowtext >= 2
%     disp(' ')
%     disp(['      absnyfivec = ',num2str(absnyfivec)])
%     disp(['      cosnyfivec = ',num2str(cosnyfivec)])
% end

%------------------------------------------------------------------
% Check which category we have.
%
% Is the apex point inside the physical wedge or not? We want to use an
% analytic approximation for a small part of the apex section (and
% an ordinary numerical integration for the rest of the wedge).
%
%       apexincluded        0 or 1
%
%

apexincluded = 1;
if sign( zw(1)*zw(2) ) == 1
    apexincluded = 0;
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

    if includepart1*includepart2 == 0 && exacthalf == 0
       analyticalsymmetry = 0;
    else
       analyticalsymmetry = 1;        
    end

%     if localshowtext
%          disp(['      Analyticalsymmetry =  ',int2str(analyticalsymmetry),' and exacthalf = ',int2str(exacthalf)])   
%     end
    
end

%----------------------------------------------------------------
% Compute the transfer function by carrying out a numerical integration for
% the entire wedge extension. 

% singularterm = [0 0 0 0];

singularterm = absnyfivec < 10*eps | abs(absnyfivec - 2*pi) < 10*eps;
useterm = 1 - singularterm;
if any(singularterm) && localshowtext
     disp(['      Singularity for term ',int2str(find(singularterm))])   
end

tf = zeros(nfrequencies,1);

% Rextra = R0 - Rstart;
% Rstarttemp = R0;

if apexincluded == 0
    for ii = 1:nfrequencies
        k = 2*pi*frequencies(ii)/cair;
        tf(ii) = quadgk(@(x)EDbetaoverm_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm),zw(1),zw(2),'RelTol',tol)*exp(-1i*k*m0);
%         *(-ny/4/pi);
%         tf(ii) = tf(ii)*exp(-1i*k*R0);
    end
    
    Rend1  = sqrt( rs^2 + (zw(1)-zs)^2 ) + sqrt( rr^2 + (zw(1)-zr)^2 );
    Rend2  = sqrt( rs^2 + (zw(2)-zs)^2 ) + sqrt( rr^2 + (zw(2)-zr)^2 );

    if Rend1 < Rend2
        zfirst = 0;
    else
       zfirst = zw(2)-zw(1); 
    end
    
else
    for ii = 1:nfrequencies
        k = 2*pi*frequencies(ii)/cair;
        if localselectanalytical == 0
            tf(ii) = quadgk(@(x)EDbetaoverm_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm),zw(1),zw(2),'RelTol',tol)*exp(-1i*k*m0);
%                 *(-ny/4/pi)*exp(-1i*k*Rextra);
        else
             
            if includepart1 == 1
                tf1 = quadgk(@(x)EDbetaoverm_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm),zw(1),-zrangespecial,'RelTol',tol);
%                     tf1 = tf1*(-ny/4/pi);
            else
               tf1 = 0; 
               if localshowtext
                   disp('   Part1 not included');
               end
            end
            if includepart2 == 1
                tf3 = quadgk(@(x)EDbetaoverm_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm),zrangespecial,zw(2),'RelTol',tol);
%                     tf3 = tf3*(-ny/4/pi);
            else
               tf3 = 0; 
               if localshowtext
                   disp('   Part2 not included');
               end
            end
%             if doublecheck == 1
%                 tf2orig = quadgk(@(x)EDbetaoverm_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,m0,useterm),-zrangespecial,zrangespecial,'RelTol',tol);                                        
%             end

            %-----------------------------------------------------------
            % The analytical approximation

            useserialexp1 = absnyfivec < 0.01;
            useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
            useserialexp = (useserialexp1 | useserialexp2) & (~singularterm);
            if localshowtext && useserialexp
                   disp('   Using serial expansion for the cosnyfivec');
            end

            if ~any(singularterm)
                cosnyfifactor = ny./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...
                            1./abs(fivec).*(useserialexp1==1) + ...
                            1./sqrt( (fivec-2*pi/ny).^2 ).*(useserialexp2==1);
                sinovercosnyfifactor =  sinnyfivec./sqrt(2*(1-cosnyfivec)).*(useserialexp==0) + ...          
                            sign(fivec).*(1- (ny*fivec).^2/6).*(useserialexp1==1) + ...
                            sign(ny*fivec-2*pi) .*(1- (ny*fivec-2*pi).^2/6).*(useserialexp2==1);
                if analyticalsymmetry == 1
                    tf2 = -sinovercosnyfifactor/pi/m0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                    if exacthalf == 1
                        tf2 = tf2/2;
                    end
                else
                    tf2 = -sinovercosnyfifactor/pi/m0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;        
                end
            else
              if localshowtext 
                   disp('   Some term is zero due to singularity');
               end
              tf2 = zeros(1,4); 
               for jj = 1:4
                   if ~singularterm(jj)
                    cosnyfifactor = ny./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...
                                1./abs(fivec(jj)).*(useserialexp1(jj)==1);
                    sinovercosnyfifactor =  sinnyfivec(jj)./sqrt(2*(1-cosnyfivec(jj))).*(useserialexp(jj)==0) + ...          
                                sign(fivec(jj)).*(1- (ny*fivec(jj)).^2/6).*(useserialexp1(jj)==1);
                    if analyticalsymmetry == 1  
                        tf2(jj) = -sinovercosnyfifactor/pi/m0.*atan(R0*zrangespecial/m0/l0*cosnyfifactor);
                        if exacthalf == 1
                            tf2 = tf2/2;
                        end
                    else
                        tf2(jj) = -sinovercosnyfifactor/pi/m0.*(atan(R0*zanalyticalstart/m0/l0*cosnyfifactor) + atan(R0*zanalyticalend/m0/l0*cosnyfifactor))/2;                            
                    end
                   end   
               end
            end
            
            tf2 = sum(tf2);
            tf2 = tf2*(-4*pi/ny);

%                 if doublecheck 
%                    relerr = abs(tf2-tf2orig)/abs(tf2orig);
%                    if relerr > 1e-4
%                        disp(['   Rel err = ',num2str(relerr)]) 
%                        end
%                 end

            tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*m0);

        end
    end
end


