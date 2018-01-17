function [tf,singularterm] = EDwedge1st_fd(cair,frequencies,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,Rstart,bc)
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
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (svensson@iet.ntnu.no) 17 Jan. 2018
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

localshowtext = 0;

if localshowtext
    disp(' ')
    disp('Reference solution for wedge 1st')
end

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

tol = 1e-11;                         % The relative accuracy for the numerical integration
                                     % The relative accuracy for the numerical integration
                                    % It was found that 1e-6 generated a
                                    % substantial error for some cases (PC
                                    % 041129)

zrelmax = 0.1;                    % The part of the edge that should use the analyt. approx.
% zrelmax = 0.00006;                    % The part of the edge that should use the analyt. approx.

za = (rs*zr+rr*zs)/(rs+rr);         % The apex point: the point on the edge with the shortest path

% Move the z-values so that za ends up in z = 0.

Dz = abs(zs-zr);
zs = zs - za;
zr = zr - za;
zw = zw - za;                       % zw is a two-element vector
za = 0;
% disp(['Positive edge-end z-value is ',num2str(zw(2))])

rs2 = rs^2;
rr2 = rr^2;

R0 = sqrt( (rs+rr)^2 + (zr-zs)^2 );

% Calculate the sinnyfi and cosnyfi values for use in the four beta-terms
% later

fivec = pi + thetas*[1 1 -1 -1] + thetar*[1 -1 1 -1];

absnyfivec = abs(ny*fivec);

sinnyfivec = sin(ny*fivec);
cosnyfivec = cos(ny*fivec);

if localshowtext
    disp(' ')
    disp(['      absnyfivec = ',num2str(absnyfivec)])
    disp(['      cosnyfivec = ',num2str(cosnyfivec)])
end

%------------------------------------------------------------------
% Check which category we have.
%
% Is the apex point inside the physical wedge or not? We want to use an
% analytic approximation for a small part of the apex section (and
% possibly an ordinary numerical integration for the rest of the apex
% section. The ordinary numerical integration will be needed if the apex
% section extends further than what allows the analytic approximation to
% be used.)
%
%       apexincluded        0 or 1
%
%
% In addition to this choice, we also want to check if we have a
% perpendicular case or a symmetrical case because these two cases can use
% a slightly simpler anaytical approximation.
%
% We call it a perpendicular case when zs = zr (which can happen only if
% apexincluded = 1!)
%
%       perpcase           0 or 1
%
% We call it a symmetrical case when rs = rr (which could also be a
% perpendicular case. A symmetrical case might not necessarily
% include the apex point. If the apex point is not included then we are not
% really interested in whether we have a symmetrical case or not since the
% analytic approximation is used only for the apex section!)
%
%       symmcase            0 or 1
%
% Since both the perpendicular case and the symmetric case can use the same
% simpler analytic approximation, we also denote those cases that are
% neither "skew cases". Again, this is relevant only when the apex is
% included.
%
%       skewcase            0 or 1

perpcase = 0;
if zs == zr
    perpcase = 1;    
end

symmcase = 0;
if rs == rr
    symmcase = 1;    
end

skewcase = 0;
if perpcase == 0 && symmcase == 0
     skewcase = 1;
end

apexincluded = 1;
if sign( zw(1)*zw(2) ) == 1
    apexincluded = 0;
end

%------------------------------------------------------------------
%

if apexincluded == 1

    % Finally the endpoint of the integral of the analytical approximation
	% (=zrangespecial) should be only up to z = 0.01
	% or whatever is specified as zrelmax for the analytic approximation

	zrangespecial = zrelmax*min([rs,rr]);

end

%----------------------------------------------------------------
% Compute the transfer function by carrying out a numerical integration for
% the entire wedge extension. 

singularterm = [0 0 0 0];

singularterm = absnyfivec < 10*eps | abs(absnyfivec - 2*pi) < 10*eps;
useterm = 1 - singularterm;
if any(singularterm) && localshowtext
     disp(['      Singularity for term ',int2str(find(singularterm))])   
end

tf = zeros(nfrequencies,1);

Rextra = R0 - Rstart;
Rstarttemp = R0;

if apexincluded == 0
    if bc(1) == 1
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            tf(ii) = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),'RelTol',tol)*(-ny/4/pi);
            tf(ii) = tf(ii)*exp(-1i*k*Rextra);
        end
    else
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            tf(ii) = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),zw(2),'RelTol',tol)*(-ny/4/pi);
            tf(ii) = tf(ii)*exp(-1i*k*Rextra);
        end        
    end
else
    if bc(1) == 1
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            tf1 = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),-zrangespecial,'RelTol',tol);
            tf1 = tf1*(-ny/4/pi);
            tf3 = quadgk(@(x)EDbetaoverml_fd(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zrangespecial,zw(2),'RelTol',tol);
            tf3 = tf3*(-ny/4/pi);
                        
            
            %-----------------------------------------------------------
            % First the analytical approximation

            rho = rr/rs;
            sinpsi = (rs+rr)/R0;
            cospsi = (zr-zs)/R0;
            
            tempfact = (1+rho)^2*sinpsi^2 - 2*rho;        

            useserialexp1 = absnyfivec < 0.01;
            useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
            useserialexp = (useserialexp1 | useserialexp2) & (~singularterm);
              
            sqrtB1vec    = ( sqrt( 2*(1-cosnyfivec) )*R0*rho/(1+rho)^2/ny    ).*(useserialexp==0) + ...
                               ( fivec.*sqrt(1-(ny*fivec).^2/12)*R0*rho/(1+rho)^2 ).*(useserialexp1==1) + ...
                               ( (fivec-2*pi/ny).*sqrt(1-(ny*fivec-2*pi).^2/12)*R0*rho/(1+rho)^2 ).*(useserialexp2==1);

            if skewcase == 0               % Use the slightly simpler formulation of the analytical approx.

                sqrtB3 = sqrt(2)*R0*rho/(1+rho)/sqrt(tempfact);

                usespecialcase = abs(sqrtB3 - sqrtB1vec) < 1e-14;

                fifactvec    = ( (1-cosnyfivec)./ny^2             ).*(useserialexp==0) + ...
                               ( fivec.^2/2.*(1-(ny*fivec).^2/12) ).*(useserialexp1==1) + ...
                               ( (fivec-2*pi/ny).^2/2.*(1-(ny*fivec-2*pi).^2/12) ).*(useserialexp2==1);

                temp1vec = sinnyfivec./( (1+rho)^2 - tempfact.*fifactvec + eps*10);


                temp1_2vec = (sinnyfivec+ 10*eps)./( (1+rho)^2 - tempfact.*fifactvec)./(sqrtB1vec+10*eps).*atan( zrangespecial./(sqrtB1vec+eps) );

                temp3vec = -1./sqrtB3.*atan( zrangespecial./sqrtB3 );

                approxintvalvec = 2/ny^2*rho*(temp1_2vec + temp1vec*temp3vec);	
                approxintvalvec = approxintvalvec.*(sinnyfivec~=0).*(usespecialcase==0);

                if any(usespecialcase)
%                     disp('SPECIAL CASE!')
                    specialcasevalue = 1/(2*sqrtB3*sqrtB3).*( zrangespecial./(zrangespecial^2 + sqrtB3*sqrtB3) + 1./sqrtB3.*atan(zrangespecial./(sqrtB3+eps)) );
                    approxintvalvec = approxintvalvec + 4*R0*R0*rho*rho*rho*sinnyfivec/ny/ny/(1+rho)^4./tempfact.*specialcasevalue.*(usespecialcase==1);
                end


            else                           % Use the more general formulation of the analytical approx.
                B3 = 2*R0*R0*rho*rho/(1+rho)/(1+rho)/tempfact;
                B1vec = sqrtB1vec.^2;
                B2 = -2*R0*(1-rho)*rho*cospsi/(1+rho)/tempfact;
                E1vec = 4*R0^2*rho^3*sinnyfivec./ny^2/(1+rho)^4./tempfact;
        % Error found 050117 by PS. The line below is wrong and should be
        % replaced by the next one.
        %         multfact = E1vec.*B2./(B1vec*B2^2 + B1vec.^2 - B3.^2);
                multfact = -E1vec.*B2./(B1vec*B2^2 + (B1vec - B3).^2);
        
        % Error found 050118 by PS: the abs in the P1 equation were missing!
                P1 = 0.5*log( abs(B3)*abs(zrangespecial^2 + B1vec)./abs(B1vec)./abs(zrangespecial^2 + B2*zrangespecial + B3 ) );
                P2 = (B1vec-B3)./(sqrtB1vec+10*eps)./B2.*atan(zrangespecial./(sqrtB1vec+10*eps));
                q = 4*B3-B2^2;

                if q > 0
                    sqrtq = sqrt(q);
                    F = 2/sqrtq*( atan( (2*zrangespecial+B2)/sqrtq ) - atan( B2/sqrtq ) );  
                elseif q < 0
                    sqrtminq = sqrt(-q);
        % Error found 050118 by PS: the abs in the F equation were missing!
                    F = 1./sqrtminq.*log( abs(2*zrangespecial+B2-sqrtminq).*abs(B2+sqrtminq)./abs(2*zrangespecial+B2+sqrtminq)./abs(B2-sqrtminq) );            
                else   % q = 0
                    F = 4*zrangespecial./B2./(2*zrangespecial+B2);
                end
                P3 = (2*(B3-B1vec)-B2^2)/2/B2.*F;

        % Error found 050118 by PC. The line below is wrong and should be replaced
        % by the one after.
        %        approxintvalvec = sum(multfact.*(P1+P2+P3))      
                approxintvalvec = multfact.*(P1+P2+P3);      
            end

            tf2 = sum(approxintvalvec.*(1-singularterm)*(-ny/2/pi));
                               
            tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*Rextra);
            
        end
    else
        error('ERROR: Not implemented the integration properly for the Dirichlet wedge yet');
        for ii = 1:nfrequencies
            k = 2*pi*frequencies(ii)/cair;
            tf1 = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zw(1),-zrangespecial,'RelTol',tol)*(-ny/4/pi);
            tf2 = 0;
            tf3 = quadgk(@(x)EDbetaoverml_fd_dirichlet(x,k,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec,Rstarttemp,useterm),zrangespecial,zw(2),'RelTol',tol)*(-ny/4/pi);
            tf(ii) = (tf1+tf2+tf3)*exp(-1i*k*Rextra);
        end        
    end
end


