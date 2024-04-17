function [ir,ninit] = EDwedge2nd(cylS,cylR,cylE2_r1,cylE1_r2,...
nyveclist,edgelengthlist,dzvec,method,pathalongplane,R_irstart,bc,cair,fs)
% EDwedge2nd - Gives the 2nd-order diffraction IR.
%
% NB!!! Only edges that are in-plane with each other can be handled with
% this version. A small extension is needed to handle non-in-plane edge
% pairs.
%
% NB!!! Known limitation/error: geometries like reflector clouds with
% several surfaces in the same plane will give some errors: obstruction
% checks are not done properly for edge-to-edge pairs if they are in the
% same plane. Some removal of non-allowed edge-to-edge pairs is done: if a
% path from one edge to another edges starts by crossing the first edge's
% own plane, then that edge-to-edge combination is shut off. This could be
% expressed as if a second edge is "behind" the first edge. It corresponds
% somewhat to the check for the image source method if an image source is
% behind or in front of its last potential reflection surface -  a
% "validity" check. But, the subsequent "obstruction" check is not done
% yet.
%
%   ir     The 2nd order diffraction impulse response for two
%          edges that share one plane.
%   ninit  The initial delay that has been cut away before the
%          beginning of ir. Note that the time zero of the ir (including ninit)
%		   corresponds to R_irstart.
%
% with input parameters
%	cylS	vector containing [rS,thetaS,zS] = cyl. coordinates of the source rel. to the first edge
%	cylR	vector containing [rR,thetaR,zR] = cyl. coordinates of the receiver rel. to the second edge
%	cylE2_r1	matrix containing [rstart_r1,thetastart_r1,zstart_r1;rend_r1,thetaend_r1,zend_r1]
%				= cyl. coordinates of the start and end points of the second edge rel. to the first edge
%	cylE1_r2	matrix containing [rstart_r2,thetastart_r2,zstart_r2;rend_r2,thetaend_r2,zend_r2]
%				= cyl. coordinates of the start and end points of the first edge rel. to the second edge
%   nyveclist
%          	vector with the two wedge indices of the two edges
%   edgelengthlist 
%           Matrix, [2 2], with the start and end values of the two edges:
%           Row 1 has the start and end values of edge 1
%           Row 2 has the start and end values of edge 2
%           For a fully visible/active edge, the start value = 0 and the
%           end value equals the length (in meters) of the edge.
%	elemsize_2nd
%		   	the number of elements per wavelength at the sampling frequency.
%		   	Recommended value 0.5, and this value is chosen if no other value
%		   	has been set. One can experiment with lower values (0.25, 0.1, ...)
%		   	for comparisons - the error increases with frequency but if first
%		   	order diffraction and specular response is masking the second order
%		   	diffraction, low values of elemsize_2nd might be perfectly acceptable.
%	method	'n' for the New method or 'v' for Vanderkooys method
%   pathalongplane      1 or 0 depending on if the path runs along a plane
%           or not
%   R_irstart (optional)
%          	a distance to which the ir time zero will correspond
%	bc (optional)
%			a matrix containing [refl11 refl12;refl21 refl22] where the values are +1 or -1
%			indicating the reflection coefficients of the four wedge planes involved.
%			+1 means rigid and -1 means soft (pressure-release)
%			refl11 is edge 1, plane 1 (refplane)
%			refl12 is edge 1, plane 2
%			refl21 is edge 2, plane 1 (refplane)
%			refl22 is edge 2, plane 2
%			default is all-rigid
%   CAIR (global)
%          	the speed of sound
%	FSAMP (global)
%			the sampling frequency
%   BIGEDGESTEPMATRIX (global)
%           a matrix, size [N,2], where N = nedge1*nedge2. nedge1 is the number of edge
%           elements that edge1 is sub-divided into and nedge2 is the
%           number of edge elements that edge2 is sub-divided into. The
%           matrix contains values between 0 and 1, representing all
%           combinations of relative positions along the two edges for the
%           mid-points of these edge elements. They are organized like
%           below, if nedge1 = 10 and nedge2 = 25:
%           BIGEDGESTEPMATRIX = 
%                  [ 0.05   0.02;
%                    0.05   0.06;
%                    0.05   0.10;
%                    ...
%                    0.05   0.98;
%                    0.15   0.02;
%                    ...
%                    0.95   0.98];
%           These values are multiplied with the lengths of the respective
%           edges to give the actual position, in meters along each edge.
%
% Uses the subroutine EDcoordtrans1
%
% Peter Svensson (peter.svensson@ntnu.no) 17 Apr 2024
%
% [ir,ninit] = EDwedge2nd(cylS,cylR,cylE2_r1,cylE1_r2,...
% nyveclist,edgelengthlist,dzvec,method,pathalongplane,BigB,R_irstart,bc,cair,fs);

% 8 May 2008: Functioning version
% 1 Sep. 2014: Fixed probable bug for edgesareinplane
% 24 March 2017 Major speedup by using accumarray for filling the
% impulse response at the end.
% 7 April 2017: tweaked the accumarray a bit; significant speedup.
% 16 Mar 2018 Copied to EDtoolbox
% 21 May 2019 Clarified a bit that theta is computed correctly for
% non-plane edge combinations - but still not if swapbigmatrix = 1?
% (which never seems to occur??)
% 17 Apr 2024 Removed (commented out) 3 lines that cleared variables, which
% cut the calculation time drastically for one text case (time 44% after
% the change for on-axis receivers and a polygonal disc).

global BIGEDGESTEPMATRIX

multfac = 1/(double(pathalongplane) + 1);

%-----------------------------------------------------------------------
% Extract the individual parameters

rS     = cylS(1);
thetaS = cylS(2);
zS     = cylS(3);

rR     = cylR(1);
thetaR = cylR(2);
zR     = cylR(3);

rE1_r1     = cylE2_r1(1,1);
thetaE1_r1 = cylE2_r1(1,2);
zE1_r1     = cylE2_r1(1,3);
rE2_r1     = cylE2_r1(2,1);
thetaE2_r1 = cylE2_r1(2,2);
zE2_r1     = cylE2_r1(2,3);

rE1_r2     = cylE1_r2(1,1);
thetaE1_r2 = cylE1_r2(1,2);
zE1_r2     = cylE1_r2(1,3);
rE2_r2     = cylE1_r2(2,1);
thetaE2_r2 = cylE1_r2(2,2);
zE2_r2     = cylE1_r2(2,3);

% Probable bug fixed Sep. 2014
% edgesareinplane  = ( abs(thetaE1_r1 - thetaE2_r1) < 1e-6 );
closwedang = 2*pi-pi/nyveclist(2);
if (abs(thetaE1_r1) < 1e-10 || abs(thetaE1_r1-(2*pi-closwedang)) < 1e-10) && abs(thetaE1_r1 - thetaE2_r1) < 1e-10
    edgesareinplane = 1;
else
    edgesareinplane = 0;
end
    
swapbigmatrix = 0;
% zmid1 = (edgelengthlist(1,1)+edgelengthlist(1,2))/2;
% zmid2_re1 = ( zE1_r1 + zE2_r1  )/2;
% midpointdist = sqrt( rE1_r1.^2  + (zmid2_re1-zmid1).^2); 
% reledgetoedgedist = midpointdist/max(dzvec);

% E2E, first part

if swapbigmatrix == 0
    B4 = edgelengthlist(1,1) + (edgelengthlist(1,2)-edgelengthlist(1,1))*BIGEDGESTEPMATRIX(:,1);      % B4 is zedge1_re1 which goes, per def., from 0 to edgelength1
    zedge2_re1 = zE1_r1 + BIGEDGESTEPMATRIX(:,2)*( zE2_r1 - zE1_r1 );     % zedge2_re1 are the z-values of edge2, expressed in edge 1
    B5 = edgelengthlist(2,1) + (edgelengthlist(2,2)-edgelengthlist(2,1))*BIGEDGESTEPMATRIX(:,2);      % B5 is zedge2_re2 which goes, per def., from 0 to edgelength2
    zedge1_re2 = zE1_r2 + BIGEDGESTEPMATRIX(:,1)*( zE2_r2 - zE1_r2 ); 
else
    B4 = edgelengthlist(1,1) + (edgelengthlist(1,2)-edgelengthlist(1,1))*BIGE(:,1);      % B4 is zedge1_re1 which goes, per def., from 0 to edgelength1
    zedge2_re1 = zE1_r1 + BIGE(:,2)*( zE2_r1 - zE1_r1 );     % zedge2_re1 are the z-values of edge2, expressed in edge 1
    B5 = edgelengthlist(2,1) + (edgelengthlist(2,2)-edgelengthlist(2,1))*BIGE(:,2);      % B5 is zedge2_re2 which goes, per def., from 0 to edgelength2
    zedge1_re2 = zE1_r2 + BIGE(:,1)*( zE2_r2 - zE1_r2 ); 
end

% S2E
S2Edist = sqrt( rS^2 + ( B4 - zS ).^2);

% E2R
E2Rdist = sqrt( rR^2 + ( B5 - zR ).^2);

% E2E, second part
if swapbigmatrix == 0
    if edgesareinplane == 0      
        
        % First we need to reconvert the cylindrical coordinates of
        % edge2-re1 into cartesian! Then we can use the EDB1coordtrans.
        % We define our own cartesian coord syst such that the reference
        % edge has its starting point in [0 0 0], and the x-axis is along
        % the reference plane.
        
        xE1_re1 = rE1_r1*cos(thetaE1_r1);
        yE1_re1 = rE1_r1*sin(thetaE1_r1);
        xE2_re1 = rE2_r1*cos(thetaE2_r1);
        yE2_re1 = rE2_r1*sin(thetaE2_r1);
        
        xvec_re1 = xE1_re1 + (xE2_re1 - xE1_re1)*BIGEDGESTEPMATRIX(:,2);
        yvec_re1 = yE1_re1 + (yE2_re1 - yE1_re1)*BIGEDGESTEPMATRIX(:,2);
        zvec_re1 = zE1_r1 + (zE2_r1 - zE1_r1)*BIGEDGESTEPMATRIX(:,2);       
        
        % NB!!! Below we need to check if really the full lengths of both
        % edges should be used. If S or R see only part of the respective
        % edge then the edge-to-edge contribution should be "windowed" too.
        [redge2_re1,thetaedge2_re1,znewedge2_re1] = EDcoordtrans1([xvec_re1 yvec_re1 zvec_re1],[0 0 0;0 0 edgelengthlist(1,2)],[0 1 0]);
         
        xE1_re2 = rE1_r2*cos(thetaE1_r2);
        yE1_re2 = rE1_r2*sin(thetaE1_r2);
        xE2_re2 = rE2_r2*cos(thetaE2_r2);
        yE2_re2 = rE2_r2*sin(thetaE2_r2);
        
        xvec_re2 = xE1_re2 + (xE2_re2 - xE1_re2)*BIGEDGESTEPMATRIX(:,1);
        yvec_re2 = yE1_re2 + (yE2_re2 - yE1_re2)*BIGEDGESTEPMATRIX(:,1);
        zvec_re2 = zE1_r2 + (zE2_r2 - zE1_r2)*BIGEDGESTEPMATRIX(:,1);
                
        % NB!!! Below we need to check if really the full lengths of both
        % edges should be used. If S or R see only part of the respective
        % edge then the edge-to-edge contribution should be "windowed" too.
        [redge1_re2,thetaedge1_re2,znewedge1_re2] = EDcoordtrans1([xvec_re2 yvec_re2 zvec_re2],[0 0 0;0 0 edgelengthlist(2,2)],[0 1 0]);   
         
    else
        redge2_re1 = rE1_r1 + BIGEDGESTEPMATRIX(:,2)*( rE2_r1 - rE1_r1 );
        redge1_re2 = rE1_r2 + BIGEDGESTEPMATRIX(:,1)*( rE2_r2 - rE1_r2 );                
        thetaedge1_re2 = thetaE1_r2;
        thetaedge2_re1 = thetaE1_r1;
    end
else
    
    if edgesareinplane == 0
         error('ERROR: r and theta calculation for non-in-plane edges not implemented yet, for swapbigmatrix=1!')

    else
        redge2_re1 = rE1_r1 + BIGE(:,2)*( rE2_r1 - rE1_r1 );
        redge1_re2 = rE1_r2 + BIGE(:,1)*( rE2_r2 - rE1_r2 );    
        thetaedge1_re2 = thetaE1_r2;
        thetaedge2_re1 = thetaE1_r1;
    end
end

if method == 'n'

	%-----------------------------------------------------------------------
	% In directivity function 2, DF2, we need the quantity
	%
	%   (  ( ze2_re2 - ze1_re2 ) * ( ze2_re2 - zr_re2 ) + n*l )/re1_re2/rr

    % First, E2E is relative to edge 2
    E2Edist = sqrt( (B5 - zedge1_re2).^2 + redge1_re2.^2 );
        
	% B2 will get the coshnyeta values for DF2
	B2 =  ( ( B5 - zedge1_re2 ).*( B5 - zR ) + E2Edist.*E2Rdist )./redge1_re2/rR;
	B2 = ( sqrt( B2.^2-1) + B2 ).^nyveclist(2);
	B2 = real( B2 + 1./B2)/2;
    
%    clear redge1_re2 zedge1_re2
    
	%-----------------------------------------------------------------------
	% In directivity function 1, DF1, we need the quantity
	%
	%   (  ( ze1_re1  -zs_re1 )*( ze1_re1 - ze2_re1 ) + n*m )/rs/re2_re1

    % Now, we need E2E relative to edge 1
    E2Edist = sqrt( (B4 - zedge2_re1).^2 + redge2_re1.^2 );

	% B1 will get the coshnyeta values for DF1
	
 	B1 = (( B4 - zS ).*(  B4 - zedge2_re1 ) + S2Edist.*E2Edist)/rS./redge2_re1;
 	B1 = ( real(sqrt( B1.^2-1)) + B1 ).^nyveclist(1);
 	B1 = ( B1 + 1./B1)/2;

 %   clear zedge2_re1 redge2_re1

else
	B1 = 1;   B2 = 1;
end
	
% Calculate DF1.*DF2

if pathalongplane == 0

    % Note that here we use E2E re. edge 1!
    B3 = multfac*nyveclist(1)*nyveclist(2)/16/pi^2*dzvec(1)*dzvec(2)*...
	              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi + thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS - thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))) + ...
					sin(nyveclist(1)*(pi - thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS - thetaedge2_re1   ))))./S2Edist./E2Edist;

    B3 = 	  B3.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))) + ...
				    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))) + ...
				    sin(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))) + ...
					sin(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))))./E2Rdist;
                
else
    
    if thetaS == 0
        B3 = multfac*nyveclist(1)*nyveclist(2)/2/pi^2*dzvec(1)*dzvec(2)*...
		              ( sin(nyveclist(1)*(pi + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaedge2_re1   ))))./S2Edist./E2Edist;
    else
        
        B3 = multfac*nyveclist(1)*nyveclist(2)/4/pi^2*dzvec(1)*dzvec(2)*...
		              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
					    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))))./S2Edist./E2Edist;
    end    

    
    B3 = 	      B3.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))) + ...
					    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))))./E2Rdist;                    
end

%-------------------------------------------------------------------------
% Determine in which sample slots the amplitude contributions should
% be added, based on the total distance.

% B2 is number of sample slots as non-integers

B2 = ( S2Edist + E2Edist + E2Rdist - R_irstart)/(cair/fs)+1;
% clear S2Edist E2Edist E2Rdist

% B1 is the sampleslot numbers in integer number
B1 = floor(B2);

% B2 is a value between 0 and 1 stating how much of dh that should be added
% to the first sample slot, i.e. the one given on B1. The rest of dh should
% be added to the following sample slot.
B2 = B2 - B1;

% New approach: simply state that the amplitudevalue, mult. by 1-B2, should
% be placed in the integer slot given in B1. In addition, the
% amplitudevalues, mult by B2, should be placed in the integer slot given
% by B1+1.

% Previous version
% B1 = [B1;B1+1];
% B3 = [B3.*(1-B2);B3.*B2];
% B1 = [B1;B1+1];

lastslot = max(B1);
firstslot = min(B1);
ir = zeros(lastslot+1,1 );

ir(firstslot:lastslot) = accumarray(B1-firstslot+1,B3.*(1-B2));

lastslot = lastslot + 1;
firstslot = firstslot+1;

ir(firstslot:lastslot) = ir(firstslot:lastslot) + accumarray(B1-firstslot+2,B3.*B2);

% Here we could have the possibility to save space by cutting out initial zeros and
% set the value ninit correspondingly to the number of removed zeros. However, the rest of
% the EDB1 functions have not implented support for this.

ninit = 0;
