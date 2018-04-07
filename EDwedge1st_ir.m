function [ir,initdelay,singularterm,zfirst] = EDwedge1st_ir(fs,cair,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,R_irstart,bc)
% EDwedge1st_ir - Gives the 1st order diffraction impulse response.
% Gives the 1st order diffraction impulse response
% for a point source irradiating a finite wedge. A free-field impulse response amplitude
% of 1/r is used. An accurate integration method is used
% so receivers can be placed right at the zone boundaries.
%
% Output parameters:
%   ir              the impulse response
%   initdelay 	    the excluded initial time delay, from the point source switch-on until the start of ir.
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
%   fs	  	        sampling frequency
%   cair            speed of sound, in m/s
%   closwedang	    closed wedge angle
%   rs, thetas, zs	cyl. coordinates of the source pos. (zs = 0)
%   rr, thetar, zr	cyl. coordinates of the receiver pos
%   zw		        the edge ends given as [z1 z2]
%   Method	        'New' (the new method by Svensson-Andersson-Vanderkooy) or 'Van'
%                   (Vanderkooys old method) or 'KDA' (Kirchhoff approximation).
%   R_irstart (optional)	if this is included, the impulse response time zero will
%				    correspond to this distance. Otherwise, the time zero
%				    corresponds to the distance R0 via the apex point.
%   bc (optional) 	boundary condition, [1 1] (default) => rigid-rigid, [-1 -1] => soft-soft
%                   [1 -1] => rigid-soft. First plane should be the ref. plane.
%   CAIR (global)   The speed of sound
%
% Uses the functions EDquadstep and EDbetaoverml for numerical integration.
%
% Peter Svensson (peter.svensson@ntnu.no) 7 Apr 2018
%
% [ir,initdelay,singularterm,zfirst] =
% EDwedge1st_ir(fs,closwedang,rs,thetas,zs,rr,thetar,zr,zw,Method,R_irstart,bc);

% 7 June 2006 FUnctioning version in EDB toolbox
% 29 Dec. 2016 Renamed to ESIE2, without making any changes.
% 28 Jan 2018 Copied from ESIE2wedge1st_int
% 7 Apr 2018 Introduced the zfirst output parameter

localshowtext = 0;

if localshowtext
    disp(' ')
    disp('Reference solution for wedge 1st')
end

if nargin < 12
	R_irstart = 0;
end
if nargin < 13
	bc = [1 1];
else
	if bc(1)*bc(2)~=1 && bc(1) ~= 1
		error('ERROR: Only all-rigid wedges are implemented');
	end
end

%------------------------------------------------------------------
% If the wedge has the length zero, return a zero IR immediately

if zw(2)-zw(1) == 0
	ir = [0 0].';
	initdelay = 0;
    singularterm = [0 0 0 0];
    zfirst = [];
	return
end

%----------------------------------------------
% Method + Some edge constants

ny = pi/(2*pi-closwedang);

initdelay = R_irstart/cair;
sampconst = cair/fs;

if Method(1) == 'N' || Method(1) == 'n'
	Method = 'New';
else
	error(['ERROR: The method ',Method,' has not been implemented yet!']);
end

tol = 1e-11;                         % The relative accuracy for the numerical integration
                                    % It was found that 1e-6 generated a
                                    % substantial error for some cases (PC
                                    % 041129)

zrelmax = 0.1;                    % The part of the edge that should use the analyt. approx.

za = (rs*zr+rr*zs)/(rs+rr);         % The apex point: the point on the edge with the shortest path

zfirst = za;

% Move the z-values so that za ends up in z = 0.

Dz = abs(zs-zr);
zs = zs - za;
zr = zr - za;
zw = zw - za;                       % zw is a two-element vector
% za = 0;

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
% 1. Is the apex point inside the physical wedge or not? We want to use an
%    analytic approximation for a small part of the apex section (and
%    possibly an ordinary numerical integration for the rest of the apex
%    section. The ordinary numerical integration will be needed if the apex
%    section extends further than what allows the analytic approximation to
%    be used.)
%
%       apexincluded        0 or 1
%
% 2. Will the IR have a step or not? A step is caused if the path lengths
%    to the two edge end points are different. If these two path lengths
%    are different, then we say that the IR has a "tail". This tail has
%    the amplitude 0.5 times the infinite wedge-IR, because there is no
%    no corresponding wedge part on the opposite side of the apex point.
%    We will use an ordinary numerical integration for the tail section.
%
%       tailincluded       0 or 1
%
% Note that if at least one of apexincluded and tailincluded must be 1.
%
% In addition to these two choices, we also want to check if we have a
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

Rendnear = sqrt( rs2 + (zw(1)-zs)^2 ) + sqrt( rr2 + (zw(1)-zr)^2 );
Rendfar  = sqrt( rs2 + (zw(2)-zs)^2 ) + sqrt( rr2 + (zw(2)-zr)^2 );

if Rendnear > Rendfar
    temp     = Rendnear;
    Rendnear = Rendfar;
    Rendfar  = temp;
end

tailincluded = 1;
if abs(Rendfar-Rendnear) < eps*1e6
    tailincluded = 0;
end

%------------------------------------------------------------------
% Find out where the IR should start and end

arrivalsampnumb = round((R0-R_irstart)/sampconst);
endsampnumb_apex = round((Rendnear-R_irstart)/sampconst);
if tailincluded == 1
    endsampnumb_tail = round((Rendfar-R_irstart)/sampconst);    
    ir = zeros(endsampnumb_tail+1,1);    
else
    ir = zeros(endsampnumb_apex+1,1);        
end

%------------------------------------------------------------------
% Construct vectors with the integration intervals along the edge
%
% For the apex section, which will include one part along the negative
% branch and one part along the positive branch of the wedge
% we carry out the integration only along the positive branch and multiply
% by two.
% For the tail section of the wedge (if it is included), we will also carry
% out the numerical intregration along the positive branch - even if the
% physical wedge actually has its tail secction along the negative branch.
%
% This choice implies means that we must always find the positive solutions
% to the equation where we find z-values along the wedge that correspond to
% the integration time values.

% Calculate the integration intervals for the apex section

if apexincluded == 1
    
    % x is the path length (m+l) which corresponds to the end of each
    % integration interval = time sample.
    % This is the part of the edge which is symmetrical around the apex point.

    x = sampconst*(arrivalsampnumb+[0.5:1:(endsampnumb_apex-arrivalsampnumb)+0.5]) + R_irstart;
    nx = length(x);
    
    if x(nx) > Rendnear
        x(nx) = Rendnear;
        if nx > 1
            if x(nx) < x(nx-1)
                error('ERROR: Strange error number 1')    
            end
        end    
    end
    x2 = x.^2;
    % ptc note: x is the same as L which is used in other derivations
    
    % Convert the x-values to z-values along the edge. These z-values will
	% define the integration intervals.

    % ptc new code
    M02   = rs*rs + zs*zs;
    L02   = rr*rr + zr*zr;
    K     = M02 - L02 - x2;
    
    diffz = (zs - zr);
    denom = diffz^2 - x2;
    a     = (2*x2*zr - K*diffz) ./ denom; 
    b     = ((K.^2 / 4) - x2*L02) ./ denom;

    zrange_apex = -a/2 + sqrt((a.^2)/4 - b);    

    % Finally the endpoint of the integral of the analytical approximation
	% (=zrangespecial) should be only up to z = 0.01
	% or whatever is specified as zrelmax for the analytic approximation

	zrangespecial = zrelmax*min([rs,rr]);
            % The explicit values below were just tested to get identical results with
            % another implementation
            %      zrangespecial = 0.01629082265142;
            %     zrangespecial = 0.01530922080939;    %0.00521457034952;
            %     zrangespecial = 0.00510940311841;
            %   disp('End of apex range = :')
            %      zrange_apex(1)
    
    if zrangespecial > zrange_apex(1)
        zrangespecial = abs(zrange_apex(1));
        splitinteg = 0;
    else
        splitinteg = 1;
	end
end

% Calculate the integration intervals for the tail section

if tailincluded == 1

    % x is the path length (m+l) which corresponds to the end of each
	% integration interval = time sample.
    % This is the part of the edge which is outside the symmetrical part around the apex point.
    
    x = sampconst*(endsampnumb_apex+[-0.5:1:(endsampnumb_tail-endsampnumb_apex)+0.5]) + R_irstart;
    nx = length(x);
    if x(nx) > Rendfar
        x(nx) = Rendfar;
        if nx > 1
            if x(nx) < x(nx-1)
                error('ERROR: Strange error number 2')    
            end
        end
    end
    if x(1) < Rendnear
        x(1) = Rendnear;
        if nx > 1
            if x(2) < x(1)
               error('ERROR: Strange error number 2')    
            end
        end
    end
    x2 = x.^2;
    
    % ptc new code
    M02   = rs*rs + zs*zs;
    L02   = rr*rr + zr*zr;
    K     = M02 - L02 - x2;
    diffz = (zs - zr);
    denom = diffz^2 - x2;
    a     = (2*x2*zr - K*diffz) ./ denom; 
    b     = ((K.^2 / 4) - x2*L02) ./ denom;

    zrange_tail = -a/2 + sqrt((a.^2)/4 - b);
end   

if localshowtext
    if apexincluded
        if tailincluded
            disp([int2str(length(zrange_apex)+length(zrange_tail)),' elements'])    
        else
            disp([int2str(length(zrange_apex)),' elements'])    
        end
    else
        disp([int2str(length(zrange_tail)),' elements'])            
    end
end

%----------------------------------------------------------------
% Compute the impulse responses by carrying out a numerical integration for
% the little segments that represent each sample slot. 
%
% For the apex section, there is one little part of the first sample slot
% which can employ an explicit integration based on the analytical
% approximation. Often, part of the first sample slot is done with the
% explicit integration value and another part is done with usual
% integration.

singularterm = [0 0 0 0];

if apexincluded == 1

    %-----------------------------------------------------------
    % First the analytical approximation

    rho = rr/rs;
	sinpsi = (rs+rr)/R0;
	cospsi = (zr-zs)/R0;
    if localshowtext
        disp(['      cospsi = ',num2str(cospsi)])    
    end
    
    tempfact = (1+rho)^2*sinpsi^2 - 2*rho;

	useserialexp1 = absnyfivec < 0.01;
    useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
    useserialexp = useserialexp1 | useserialexp2;
    
    singularterm = absnyfivec < 10*eps | abs(absnyfivec - 2*pi) < 10*eps;
    if any(singularterm) && localshowtext
        disp(['      Singularity for term ',int2str(find(singularterm))])    
    end
    
    sqrtB1vec    = ( sqrt( 2*(1-cosnyfivec) )*R0*rho/(1+rho)^2/ny    ).*(useserialexp==0) + ...
                       ( fivec.*sqrt(1-(ny*fivec).^2/12)*R0*rho/(1+rho)^2 ).*(useserialexp1==1) + ...
                       ( (fivec-2*pi/ny).*sqrt(1-(ny*fivec-2*pi).^2/12)*R0*rho/(1+rho)^2 ).*(useserialexp2==1);

    if skewcase == 0              % Use the slightly simpler formulation of the analytical approx.
        
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
            disp('SPECIAL CASE!')
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
% Wrong version:        multfact = E1vec.*B2./(B1vec*B2^2 + B1vec.^2 - B3.^2);
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
% Wrong version:        approxintvalvec = sum(multfact.*(P1+P2+P3))      
        approxintvalvec = multfact.*(P1+P2+P3);      
    end

    approxintvalvec = sum(approxintvalvec.*(1-singularterm));
    if localshowtext
        disp(['      Analyt. int = (final IR value) ',num2str(sum(approxintvalvec*(-ny/2/pi)),15)])
    end
    
    %-----------------------------------------------------------
    % Numerical integration for the rest of the apex section

    % The first IR sample will either have part analytic, part numerical
    % integration, or just analytic
    
    ir(1+arrivalsampnumb) = ir(1+arrivalsampnumb) + approxintvalvec;     

    if splitinteg == 1
        h = 0.13579*(zrange_apex(1)-zrangespecial);
        x = [zrangespecial zrangespecial+h zrangespecial+2*h (zrangespecial+zrange_apex(1))/2 zrange_apex(1)-2*h zrange_apex(1)-h zrange_apex(1)];
        y = EDbetaoverml(x,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
        
        Q = [0 0 0];
        Q(:,1) = EDquadstep(x(:,1),x(:,3),y(:,1),y(:,2),y(:,3),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
        Q(:,2) = EDquadstep(x(:,3),x(:,5),y(:,3),y(:,4),y(:,5),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
        Q(:,3) = EDquadstep(x(:,5),x(:,7),y(:,5),y(:,6),y(:,7),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
        ir(1+arrivalsampnumb) = ir(1+arrivalsampnumb) + sum(Q);     

        if localshowtext
            disp(['      Split: num value = (final IR value) ',num2str(sum(Q)*(-ny/2/pi)),15])
            disp(['      Arrival sample = ',int2str(1+arrivalsampnumb)])
        end
	end
    if localshowtext
        disp(['      First IR value  = ',num2str(ir(1+arrivalsampnumb)*(-ny/2/pi),15)])
        disp(' ')
    end
    
    % The numerical integration for the rest of the apex section
    
    zstarts = zrange_apex(1:end-1).';
	zends = zrange_apex(2:end).';
    
	h = 0.13579*(zends-zstarts);
	nslots = length(h);
	x = [zstarts zstarts+h zstarts+2*h (zstarts+zends)/2 zends-2*h zends-h zends];
	x = reshape(x,7*nslots,1);
	y = EDbetaoverml(x,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	y = reshape(y,nslots,7);
	x = reshape(x,nslots,7);
	Q = zeros(nslots,3);

    Q(:,1) = EDquadstep(x(:,1),x(:,3),y(:,1),y(:,2),y(:,3),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	Q(:,2) = EDquadstep(x(:,3),x(:,5),y(:,3),y(:,4),y(:,5),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	Q(:,3) = EDquadstep(x(:,5),x(:,7),y(:,5),y(:,6),y(:,7),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	ir(arrivalsampnumb+[2:length(zrange_apex)]) = ir(arrivalsampnumb+[2:length(zrange_apex)]) + sum(Q.').';
	ir = ir*(-ny/2/pi);   % Mult by 2 because analyt int. is half the wedge

else

    Rend1  = sqrt( rs^2 + (zw(1)-zs)^2 ) + sqrt( rr^2 + (zw(1)-zr)^2 );
    Rend2  = sqrt( rs^2 + (zw(2)-zs)^2 ) + sqrt( rr^2 + (zw(2)-zr)^2 );

    if Rend1 < Rend2
        zfirst = 0;
    else
       zfirst = zw(2)-zw(1); 
    end
    
end

%--------------------------------------------------------------------------
% If we have a non-symmetrical case, then we must integrate for the long
% end of the edge as well.

if tailincluded == 1
	zstarts = zrange_tail(1:end-1).';
	zends = zrange_tail(2:end).';    
    
	h = 0.13579*(zends-zstarts);
	nslots = length(h);
	x = [zstarts zstarts+h zstarts+2*h (zstarts+zends)/2 zends-2*h zends-h zends];
	x = reshape(x,7*nslots,1);
	y = EDbetaoverml(x,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	y = reshape(y,nslots,7);
	x = reshape(x,nslots,7);
	Q = zeros(nslots,3);
	Q(:,1) = EDquadstep(x(:,1),x(:,3),y(:,1),y(:,2),y(:,3),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	Q(:,2) = EDquadstep(x(:,3),x(:,5),y(:,3),y(:,4),y(:,5),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
	Q(:,3) = EDquadstep(x(:,5),x(:,7),y(:,5),y(:,6),y(:,7),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
    
    ir(endsampnumb_apex+[1:length(zrange_tail)-1]) = ir(endsampnumb_apex+[1:length(zrange_tail)-1]) + (sum(Q.').')*(-ny/4/pi);

end

ir = real(ir);

if localshowtext
    disp(['      and finally,'])
    disp(['      the first IR values  = ',num2str(ir(0+arrivalsampnumb),15)])
    disp(['                             ',num2str(ir(1+arrivalsampnumb),15)])
    disp(['                             ',num2str(ir(2+arrivalsampnumb),15)])
    disp(' ')
end
