function [ir,ninit] = EDwedgeN(cylS,cylR,cylE,ncylrows,nyveclist,edgelengthlist,dzvec,method,...
pathalongplane,nedgeelcombs,R_irstart,bc,cair,fs,BIGEDGE1stvalue)
% EDwedgeN - Gives the Nth-order diffraction impulse response.
%   ir     The Nth order diffraction impulse response for a combination of edges
%   ninit  The initial delay that has been cut away before the
%          beginning of ir. Note that the time zero of the ir (including ninit)
%		   corresponds to R_irstart.
%
% with input parameters
%	cylS	        List containing [rS,thetaS,zS] = cyl. coord. of the source rel. to the first edge
%	cylR	        List containing [rR,thetaR,zR] = cyl. coord. of the receiver rel. to the last edge
%   cylE            Matrix containing consecutive blocks of 
%                       [rstartE2_r1,thetastartE2_r1,zstartE2_r1;rendE2_r1,thetaendE2_r1,zendE2_r1]
%                       [rstartE1_r2,thetastartE1_r2,zstartE1_r2;rendE1_r2,thetaendE1_r2,zendE1_r2]
%                       etc (E3_r2,E2_r3)
%   ncylrows        The number of rows in cylE
%   nyveclist       List of the wedge indices of the N edges
%   edgelengthlist  List of the lengths of the N edges
%   dzvec           List of the N edge element sizes
%	method	        'n' for the New method or 'v' for Vanderkooys method
%   pathalongplane  List of 1 or 0, depending on if the (N-1) paths run along a plane or not
%   nedgeelcombs    The value prod(ndiv)
%   R_irstart       A distance to which the ir time zero will correspond
%	bc              Matrix containing [refl11 refl12 refl21 refl22] where the values are +1 or -1
%			            indicating the reflection coefficients of the four wedge planes involved.
%		            	+1 means rigid and -1 means soft (pressure-release)
%		            	refl11 is edge 1, plane 1 (refplane)
%			            refl12 is edge 1, plane 2
%		            	refl21 is edge 2, plane 1 (refplane)
%			            refl22 is edge 2, plane 2
%			            default is all-rigid
%   cair            The speed of sound
%   fs              The sampling frequency
%   BIGEDGE1stvalue The constant value that should be in the first column
%                   of BIGEDGESTEPMATRIX
%
%   BIGEDGESTEPMATRIX (global)
%           a matrix, size [Nelem,N], where Nelem = nedge1*nedge2*nedge3*....
%           nedge1 is the number of edge
%           elements that edge1 is sub-divided into and nedge2 is the
%           number of edge elements that edge2 is sub-divided into, etc. The
%           matrix contains values between 0 and 1, representing all
%           combinations of relative positions along the two edges for the
%           mid-points of these edge elements. They are organized like
%           below, if nedge1 = 10 and nedge2 = 10 and nedge3 = 25:
%           BIGEDGESTEPMATRIX = 
%                  [ 0.05   0.05  0.02;
%                    0.05   0.05  0.06;
%                    0.05   0.05  0.10;
%                    ...
%                    0.05   0.05  0.98;
%                    0.05   0.15  0.02;
%                    ...
%                    0.95   0.95  0.98];
%           These values are multiplied with the lengths of the respective
%           edges to give the actual position, in meters along each edge.
%           New version: The first column isn't included in the matrix, but
%           transferred separately as BIGEDGE1stvalue (see above).
%
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@ntnu.no) 21 May 2019
%
% [ir,ninit] = EDwedgeN(cylS,cylR,cylE,ncylrows,nyveclist,edgelengthlist,dzvec,method,...
% pathalongplane,nedgeelcombs,R_irstart,bc,cair,fs,BIGEDGE1stvalue);

% 1 Nov. 2006 Functioning version
% 24 Mar. 2017 Implemented accumarray for filling the impulse response in
%              the end. Gives major speedup.
% 7 April 2017 Small tweaking of the accumarray code.
% 16 Mar 2018 Copied over to EDtoolbox
% 21 May 2019 Implemented the handling of non-planar edge sequences (the
%             recalculation of theta angles).

global BIGEDGESTEPMATRIX

Ndifforder = length(edgelengthlist);

multfac = prod(1./(double(pathalongplane) + 1));

%-----------------------------------------------------------------------
% Create matrices with the z-, r- and theta-values
%
% BIGEDGESTEPMATRIX is a matrix that contains relative steps, ]0,1[, along each edge
%
% All the big matrices will have one column for each "edge n rel. to edge
% m" combination:
%   Column 1    E2_RE1
%   Column 2    E1_RE2
%   Column 3    E3_RE2
%   Column 4    E2_RE3
%   Column 5    E4_RE3
%   Column 6    E3_RE4
%   Column 7    E5_RE4
%   Column 8    E4_RE5
%       etc
% For Ndifforder = 2, there are 2 columns
% For Ndifforder = 3, there are 4 columns
% For Ndifforder = 4, there are 6 columns
%       etc

onesvec = uint8(ones(nedgeelcombs,1));

colvec = [ [2:Ndifforder].' [1:Ndifforder-1].'];
colvec = reshape(colvec.',1,(Ndifforder-1)*2);
ivecwithhole = [1 [3:length(colvec)]];
colvec = colvec(ivecwithhole)-1;

%-----------------------------------------------------------------------
% The z-values for each edge, relative to each own edge can be calculated
% from BIGEDGESTEPMATRIX.

zEn_REn = BIGEDGESTEPMATRIX.*edgelengthlist(onesvec,2:end);
zEn_REn1stvalue = BIGEDGE1stvalue.*edgelengthlist(1);

% The z-values for each edge, relative to the two neighbour edges

zreSTA = cylE(1:2:ncylrows,3).';
zreEND = cylE(2:2:ncylrows,3).';
zreDIF = zreEND - zreSTA;

zEn_REm = zreSTA(onesvec,ivecwithhole) + BIGEDGESTEPMATRIX(:,colvec).*zreDIF(onesvec,ivecwithhole);
zEn_REm2ndvalue = zreSTA(2) + BIGEDGE1stvalue.*zreDIF(2);

%-----------------------------------------------------------------------
% The S2E and E2R distances

% S2E
S2Edist = sqrt( cylS(1)^2 + ( zEn_REn1stvalue - cylS(3) ).^2);

% E2R
E2Rdist = sqrt( cylR(1)^2 + ( zEn_REn(:,Ndifforder-1) - cylR(3) ).^2);

%-----------------------------------------------------------------------
% The r-values for each edge, relative to the two neighbour edges

rreSTA = cylE(1:2:ncylrows,1).';
rreEND = cylE(2:2:ncylrows,1).';
rreDIF =  rreEND - rreSTA;

rEn_REm = rreSTA(onesvec,ivecwithhole) + BIGEDGESTEPMATRIX(:,colvec).*rreDIF(onesvec,ivecwithhole);
rEn_REm2ndvalue = rreSTA(2) + BIGEDGE1stvalue.*rreDIF(2);

%-----------------------------------------------------------------------
% The edge-to-edge distances

% E2Edist will be "edge2 re 1","edge3 re 2" etc
modcols = [1 [3:2:Ndifforder*2-3]-1];
E2Edist = sqrt( ([zEn_REn1stvalue*double(onesvec) zEn_REn(:,1:Ndifforder-2)] - zEn_REm(:,modcols)).^2 + rEn_REm(:,modcols).^2 );

%-----------------------------------------------------------------------
% Calculate the coshnyeta values first, to use them in the directivity
% functions

if method(1) == 'n'

    % CH will get the coshnyeta values for the directivity functions
    
    CH = zeros(nedgeelcombs,Ndifforder-1);

    %-----------------------------------------------------------------------
	% For the coshnyeta values we need the quantity
	%
	%   (  ( ze1_RE1  -zs_RE1 )*( ze1_RE1 - ze2_RE1 ) + n*m )/rs/re2_RE1

    CH(:,1) = (( zEn_REn1stvalue - cylS(3) ).*(  zEn_REn1stvalue - zEn_REm(:,1) ) + S2Edist.*E2Edist(:,1))/cylS(1)./rEn_REm(:,1);
 	CH(:,1) = ( real(sqrt( CH(:,1).^2-1)) + CH(:,1) ).^nyveclist(1);
 	CH(:,1) = ( CH(:,1) + 1./CH(:,1))/2;

    for ii = 2:Ndifforder-1
        mcol1 = ii*2-3;
        mcol2 = ii*2-2;
        if ii == 2
            CH(:,ii) =  ( ( zEn_REn(:,ii-1) - zEn_REm2ndvalue ).*(  zEn_REn(:,ii-1) - zEn_REm(:,mcol2)  ) + E2Edist(:,ii-1).*E2Edist(:,ii) )./rEn_REm2ndvalue./rEn_REm(:,mcol2);
        else
             CH(:,ii) =  ( ( zEn_REn(:,ii-1) - zEn_REm(:,mcol1) ).*(  zEn_REn(:,ii-1) - zEn_REm(:,mcol2)  ) + E2Edist(:,ii-1).*E2Edist(:,ii) )./rEn_REm(:,mcol1)./rEn_REm(:,mcol2);
            
        end
		CH(:,ii) = ( sqrt( CH(:,ii).^2-1) + CH(:,ii) ).^nyveclist(ii);
		CH(:,ii) = real( CH(:,ii) + 1./CH(:,ii))/2;
    end    
        
    CHN =  ( ( zEn_REn(:,Ndifforder-1) - zEn_REm(:,2*Ndifforder-3) ).*( zEn_REn(:,Ndifforder-1) - cylR(3) ) + E2Edist(:,Ndifforder-1).*E2Rdist )./rEn_REm(:,2*Ndifforder-3)./cylR(1);
	CHN = ( sqrt( CHN.^2-1) + CHN ).^nyveclist(Ndifforder);
	CHN = real( CHN + 1./CHN)/2;

elseif method(1) == 'v'
    CH = ones(1,Ndifforder-1);
    CHN = 1;
else
    error('Method not allowed')
end

%-----------------------------------------------------------------------
% The theta-values for each edge, relative to the two neighbour edges

thetareSTA = cylE(1:2:ncylrows,2).';
thetareDIF = cylE(2:2:ncylrows,2).' - thetareSTA;

if all(thetareDIF==0)
    thetaEn_REm = thetareSTA(ivecwithhole);
else    
    % For each edge-to-edge leg (that is, each column), we need to 
    % reconvert the cylindrical coordinates of
    % edge2-re1 into cartesian. Then we can use the EDB1coordtrans.
    % We define our own cartesian coord syst such that the reference
    % edge has its starting point in [0 0 0], and the x-axis is along
    % the reference plane.       
    % We should in principle check if really the full lengths of both
    % edges should be used.
    xSTA_REm = rreSTA(ivecwithhole).*cos(thetareSTA(ivecwithhole));
    ySTA_REm = rreSTA(ivecwithhole).*sin(thetareSTA(ivecwithhole));
    xEND_REm = rreEND(ivecwithhole).*cos(thetareSTA(ivecwithhole));
    yEND_REm = rreEND(ivecwithhole).*sin(thetareSTA(ivecwithhole));

    xvec_REm = xSTA_REm(onesvec,:) + (xEND_REm(onesvec,:) - xSTA_REm(onesvec,:)).*BIGEDGESTEPMATRIX(:,colvec);
    yvec_REm = ySTA_REm(onesvec,:) + (yEND_REm(onesvec,:) - ySTA_REm(onesvec,:)).*BIGEDGESTEPMATRIX(:,colvec);
    zvec_REm = zreSTA(onesvec,ivecwithhole)   + (zreEND(onesvec,ivecwithhole) - zreSTA(onesvec,ivecwithhole)).*BIGEDGESTEPMATRIX(:,colvec);        

    thetaEn_REm = zeros( size(BIGEDGESTEPMATRIX(:,colvec)) );

    edgerefnumbers = [2:Ndifforder;1:Ndifforder-1];
    edgerefnumbers = reshape(edgerefnumbers,prod(size(edgerefnumbers)),1);
    edgerefnumbers = edgerefnumbers(ivecwithhole);
    
    for colnumber = 1:size(thetaEn_REm,2)
        wedgecoords = [0 0 0;0 0 edgelengthlist(edgerefnumbers(colnumber))];
        [~,theta_onecol,~] = EDcoordtrans1([xvec_REm(:,colnumber) yvec_REm(:,colnumber) zvec_REm(:,colnumber)],wedgecoords,[0 1 0]);
        thetaEn_REm(:,colnumber) = theta_onecol;
    end
    
end

thetaEn_REm2ndvalue = thetareSTA(2) + BIGEDGE1stvalue.*thetareDIF(2);

%------------------------------------------------------
% Calculate the product of all directivity factors and inverse distances

% First "leg" from source, via the first edge, to the second edge, but
% including the DF only for the first edge

if pathalongplane(1) == 0
    AMP = multfac*prod(nyveclist)*prod(dzvec)*(-1/4/pi)^Ndifforder*(sin(nyveclist(1)*(pi + cylS(2) + thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi + cylS(2) + thetaEn_REm(:,1)   ))) + ...
	        sin(nyveclist(1)*(pi + cylS(2) - thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi + cylS(2) - thetaEn_REm(:,1)   ))) + ...
		    sin(nyveclist(1)*(pi - cylS(2) + thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi - cylS(2) + thetaEn_REm(:,1)   ))) + ...
		    sin(nyveclist(1)*(pi - cylS(2) - thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi - cylS(2) - thetaEn_REm(:,1)   ))))./S2Edist./E2Edist(:,1);
else
    if cylS(2) == 0
        AMP = 4*multfac*prod(nyveclist)*prod(dzvec)*(-1/4/pi)^Ndifforder*(sin(nyveclist(1)*(pi + thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi + thetaEn_REm(:,1)   ))))./S2Edist./E2Edist(:,1);
    else
         AMP = 2*multfac*prod(nyveclist)*prod(dzvec)*(-1/4/pi)^Ndifforder*(sin(nyveclist(1)*(pi + cylS(2) + thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi + cylS(2) + thetaEn_REm(:,1)   ))) + ...
			    sin(nyveclist(1)*(pi - cylS(2) + thetaEn_REm(:,1)  ))./(CH(:,1)- cos(nyveclist(1)*(pi - cylS(2) + thetaEn_REm(:,1)   ))))./S2Edist./E2Edist(:,1);
    end
end

% Second "leg" from edge 1, via the second edge, to the third edge
% including the DF only for the second edge. We keep this out of the
% for-loop because the matrix thetaEn_REm2ndvalue gives a special
% formulation.

for ii = 2:2
    iv2 = 2; 
    if pathalongplane(ii-1)*pathalongplane(ii) == 0
        AMP = AMP.*(sin(nyveclist(ii)*(pi + thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))) + ...
					sin(nyveclist(ii)*(pi + thetaEn_REm2ndvalue - thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm2ndvalue - thetaEn_REm(:,iv2)))) + ...
				    sin(nyveclist(ii)*(pi - thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi - thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))) + ...
				    sin(nyveclist(ii)*(pi - thetaEn_REm2ndvalue - thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi - thetaEn_REm2ndvalue - thetaEn_REm(:,iv2)))))./E2Edist(:,2);
    else
        AMP = 4*AMP.*(sin(nyveclist(ii)*(pi + thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm2ndvalue + thetaEn_REm(:,iv2)))))./E2Edist(:,2);    
    end            
end

% All the following "legs", except the last one that reaches the receiver.

for ii = 3:Ndifforder-1
    iv1 = ii*2-3;
    iv2 = iv1+1;
    if pathalongplane(ii-1)*pathalongplane(ii) == 0
        AMP = AMP.*(sin(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))) + ...
					sin(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) - thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) - thetaEn_REm(:,iv2)))) + ...
				    sin(nyveclist(ii)*(pi - thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi - thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))) + ...
				    sin(nyveclist(ii)*(pi - thetaEn_REm(:,iv1) - thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi - thetaEn_REm(:,iv1) - thetaEn_REm(:,iv2)))))./E2Edist(:,ii);
    else
        AMP = 4*AMP.*(sin(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))./(CH(:,ii)-cos(nyveclist(ii)*(pi + thetaEn_REm(:,iv1) + thetaEn_REm(:,iv2)))))./E2Edist(:,ii);    
    end            
end

% The last "leg" that reaches the receiver, and including the last
% directivity factor.

if pathalongplane(Ndifforder-1) == 0
    AMP = AMP.*(sin(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))) + ...
		    sin(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))) + ...
		    sin(nyveclist(Ndifforder)*(pi - thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi - thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))) + ...
		    sin(nyveclist(Ndifforder)*(pi - thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi - thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))))./E2Rdist;    
else
    AMP = 2*AMP.*(sin(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) + cylR(2)))) + ...
		          sin(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))./(CHN-cos(nyveclist(Ndifforder)*(pi + thetaEn_REm(:,2*Ndifforder-3) - cylR(2)))))./E2Rdist;        
end

%-------------------------------------------------------------------------
% Determine in which sample slots the amplitude contributions should
% be added, based on the total distance.

% B2 is number of sample slots as non-integers

if Ndifforder == 2
    B2 = ( S2Edist + E2Edist + E2Rdist - R_irstart)/(cair/fs)+1;
else
    B2 = ( S2Edist + sum(E2Edist.').' + E2Rdist - R_irstart)/(cair/fs)+1;
end

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

% B1 = [B1;B1+1];
% AMP = [AMP.*(1-B2);AMP.*B2];

% ir = zeros(max(B1),1 );

lastslot = max(B1);
firstslot = min(B1);
ir = zeros(lastslot+1,1 );

ir(firstslot:lastslot) = accumarray(B1-firstslot+1,AMP.*(1-B2));
% % % ir(firstslot:lastslot) = accumarray(B1-firstslot+1,AMP);

lastslot = lastslot + 1;
firstslot = firstslot+1;

ir(firstslot:lastslot) = ir(firstslot:lastslot) + accumarray(B1-firstslot+2,AMP.*B2);

% The four lines below were the original version, which were replaced by
% the single line using accumarray. This leads to a major speedup (one
% case, w 4 sources and 22 receivers, was tested, and gave a 35% speedup).
% [B1,sortvec] = sort(B1);
% AMP = cumsum(AMP(sortvec));
% ivstep = (find(diff([(B1);(B1(end))+1])));
% ir(B1(ivstep)) = diff([0;AMP(ivstep)]);

% ir = accumarray(B1,AMP);

% Here we could have the possibility to save space by cutting out initial zeros and
% set the value ninit correspondingly to the number of removed zeros. However, the rest of
% the EDB1 functions have not implented support for this.

ninit = 0;
