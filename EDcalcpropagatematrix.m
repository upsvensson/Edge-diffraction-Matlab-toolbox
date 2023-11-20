function [Fmatrix,ivproblematic] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
isthinhole,vispartedgesfromr,frequency,Rstart,rRvec,thetaRvec,zRvec,doesQsegmenthavevalues,showtext)
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
%
% Output parameters:
%   Fmatrix     
%   ivproblematic
%
% Uses EDcoordtrans1. 
%
% Peter Svensson (peter.svensson@ntnu.no) 2 Nov 2023
%
% [Fmatrix,ivproblematic] = EDcalcpropagatematrix(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,...
%     isthinhole,vispartedgesfromr,frequency,Rstart,...
%     rRvec,thetaRvec,zRvec,doesQsegmenthavevalues,showtext);

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

if nargin < 13
    showtext = 0;
end

% ----------------------------------------------------------------------------------------------
% Create the F-matrix by first creating a matrix showing which
% submatrices should be included.

nbig = Hsubmatrixdata.bigmatrixendnums(end);
Fmatrix = zeros(1,nbig);
refto = edgetoedgedata.reftoshortlistE;

% k = 2*pi*frequency/CAIR;
k = 2*pi*frequency/envdata.cair;

if showtext >= 3
    disp(' ')
    disp(['      Building the F matrix of size (1,',int2str(nbig),')'])
end

n1previous = 0;
n2previous = 0;
edge2previous = 0;
edge1previous = 0;

alledgepairsthin = prod(double(Hsubmatrixdata.isthinplaneedgepair));
noedgepairsthin  = prod(1-double(Hsubmatrixdata.isthinplaneedgepair));

ivrelevantedgepairs = find( vispartedgesfromr(Hsubmatrixdata.edgepairlist(:,1)) > 0 );

ivproblematic = [];

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
                disp(['         Startcol in F-matrix   = (',int2str(Hsubmatrixdata.bigmatrixstartnums(ii)),')'])
                disp(['         Endcol in F-matrix   = (',int2str(Hsubmatrixdata.bigmatrixendnums(ii)),')'])
            end

            ny = pi/(2*pi-edgedata.closwedangvec(edge2));

            if n2 ~= n2previous 
                n2vec = Hsubmatrixdata.quadraturematrix_pos(n2,1:n2);
                n2vec = n2vec(:);
                weightvec2 = Hsubmatrixdata.quadraturematrix_weights(n2,1:n2);
                weightvec2 = weightvec2(:);

                if edge2 ~= edge2previous
                    dzvec2 = weightvec2*edgedata.edgelengthvec(edge2);   
                    ze2_re2 = (0 + n2vec*len2);
                    edge2previous = edge2;
                end

            end

            if n1 ~= n1previous 
                n1vec = Hsubmatrixdata.quadraturematrix_pos(n1,1:n1);
                n1vec = n1vec(:);
                weightvec1 = Hsubmatrixdata.quadraturematrix_weights(n1,1:n1);
                weightvec1 = weightvec1(:);
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

            ldist = (zR - ze2_re2(n2vertmat)).^2   + rR.^2;
            ldist = sqrt(ldist);

            mdist = (ze1_re2(n1vertmat) - ze2_re2(n2vertmat)).^2   + re1_re2(n1vertmat).^2;
            mdist = sqrt(mdist);

            dz_expjkmloverml = exp(-1i*k*(ldist+mdist-Rstart))./ldist./mdist.*dzvec2(n2vertmat).*dzvec1(n1vertmat);       

            ch = ((ze2_re2(n2vertmat) - ze1_re2(n1vertmat)).*(ze2_re2(n2vertmat) - zR) + mdist.*ldist)./re1_re2(n1vertmat)./rR;
            ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
            ch = ( ch + 1./ch)/2;

            % 28 March 2011: thetain is always zero (or = theta_wedge) so terms
            % 1 and 3 will always be identical; and terms 2 and 4.

            if isthinhole == 0
                if acrossface_in == 1            
                    beta = 2* ( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))) + ...
                         sin(ny*(pi + thetain - thetaout  ))./(ch - cos(ny*(pi + thetain - thetaout   ))));  

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

            if showtext >= 3
                fivec = pi + thetain*[1 1 -1 -1] + thetaout*[1 -1 1 -1];
                absnyfivec = abs(ny*fivec); 
                useserialexp1 = absnyfivec < 0.01;
                useserialexp2 = abs(absnyfivec - 2*pi) < 0.01;
                useserialexp = useserialexp1 | useserialexp2;
                if any(useserialexp)
                    disp('WARNING! Analytical approximation should be used')
                end
                plot(abs(Fsub))

            end

            Fmatrix(Hsubmatrixdata.bigmatrixstartnums(ii):Hsubmatrixdata.bigmatrixendnums(ii)) =     Fmatrix(Hsubmatrixdata.bigmatrixstartnums(ii):Hsubmatrixdata.bigmatrixendnums(ii)) + Fsub.';

            if ~isempty(ivaddition)
                sectionofbigmatrix = [Hsubmatrixdata.bigmatrixstartnums(ii):Hsubmatrixdata.bigmatrixendnums(ii)].';            
                ivproblematic = [ivproblematic;sectionofbigmatrix(ivaddition)];
            end
    end
end

