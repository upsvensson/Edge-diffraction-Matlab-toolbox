function Qfirstterm = EDinteg_souterm(envdata,edgedata,edgetoedgedata,...
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
%   gaussvectors
%   rSvec,thetaSvec,zSvec
%   showtext (optional)     0 -> no text displayed on screen. Default: 0
%
% Output parameters:
%   Qfirstterm      Column vector with all edge source signals caused by
%                   the external source. All "via-edge" to "to-edge" combos
%                   are stacked after each other. For each such combo, the
%                   sequence is:
%                       to-edge-element     via-edge-element
%                            1                      1
%                            1                      2
%                            1                      3
%                   etc
%
% Uses functions EDdistelements  EDcoordtrans1
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by
%   the Free Software     
%   Foundation, either version 3 of the License, or (at your option) afny later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson 15 Dec. 2017 (peter.svensson@ntnu.no)
%
% Qfirstterm = EDinteg_souterm(envdata,edgedata,edgetoedgedata,Hsubmatrixdata,inteq_ngauss,inteq_discretizationtype,...
%     vispartedgesfroms,vispartedgesfroms_start,vispartedgesfroms_end,frequency,gaussvectors,rSvec,thetaSvec,zSvec,showtext);

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

if nargin < 14
    showtext = 0;
end

showtext = 3;

discretizationtype = controlparameters.discretizationtype;

%  symmetrycompression = 1;

% ----------------------------------------------------------------------------------------------
% Create the G-matrix

G = zeros(Hsubmatrixdata.bigmatrixendnums(end),1);
refto = edgetoedgedata.reftoshortlistE;

k = 2*pi*frequency/envdata.cair;

if showtext >= 2
    disp(' ')
    disp('      Calculating the ESIE source term')%: symmetrycompression = ',int2str(symmetrycompression)])    
    disp(['      Building the G matrix of size (',int2str(Hsubmatrixdata.bigmatrixendnums(end)),',1)'])
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
        if showtext >= 3
            disp(' ')
            disp(['         From S via edge ',int2str(edge2),' to edge ',int2str(edge3)])
            disp(['         Startrow in G-matrix   = (',int2str(Hsubmatrixdata.bigmatrixstartnums(ii)),'). Endrow in G-matrix   = (',int2str(Hsubmatrixdata.bigmatrixendnums(ii)),')'])
%             nfrom = Hsubmatrixdata.nedgeelems(edge2);
%             nto   = Hsubmatrixdata.nedgeelems(edge3);
        end
        
        ny = nyvec(edge2);
        len2 = edgedata.edgelengthvec(edge2);
        n2 = Hsubmatrixdata.nedgeelems(edge2);
        n3 = Hsubmatrixdata.nedgeelems(edge3);
                
        if size(gaussvectors,1) == n2
            n2vec = gaussvectors(:,1);
        else
           [n2vec,~] = EDdistelements(n2,discretizationtype);
        end
        % New 4 Nov. 2017: here we zero out the edge points that are not
        % visible.
        if vispartedgesfroms(edge2) ~= 1
            n2vec = n2vec.*(n2vec >= vispartedgesfroms_start(edge2)).*(n2vec <= vispartedgesfroms_end(edge2));
        end
        if size(gaussvectors,1) == n3
            n3vec = gaussvectors(:,1);
        else
           [n3vec,~] = EDdistelements(n3,discretizationtype);
        end
        
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
        if showtext >= 3
            disp(['         Acrossface_out = ',int2str(acrossface_out)])
        end

%         edgecoords = [edgedata.edgestartcoords(edge2,:);edgedata.edgeendcoords(edge2,:)];
        rS = rSvec(edge2);
        thetaS = thetaSvec(edge2);
        zS = zSvec(edge2);

        singwarning = 0;
        singterm = [0 0 0 0];
        if thetaout == 0 && max([  abs(cos(ny*(pi + thetaS   )))  abs(cos(ny*(pi - thetaS   )))]) > 0.999
            singwarning = 1 
        end
        if thetaout ~= 0 && max([ abs(cos(ny*(pi + thetaS + thetaout   )))  abs(cos(ny*(pi + thetaS - thetaout   )))  ...
                          abs(cos(ny*(pi - thetaS + thetaout   )))  abs(cos(ny*(pi - thetaS - thetaout   )))  ]) > 0.999
            singwarning = 2                       
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

                       beta1 = 2*sin(ny*(pi + thetaS  ))./(ch - cos(ny*(pi + thetaS   )));
                       beta2 = 2*sin(ny*(pi - thetaS  ))./(ch - cos(ny*(pi - thetaS   )));
                       iv1 = abs(beta1)>100;
                       iv2 = abs(beta2)>100;
                       
                       if sum(iv1) > 0
                          singterm(1) = 1; 
                       end
                       if sum(iv2) > 0
                          singterm(2) = 1; 
                       end
                       
                       if sum(iv1) > 0                       
                          plot([beta1 beta2],'-o')
                          disp('iv1')
                          sum(iv1)
                          pause
                       end
                       if sum(iv2) > 0                       
                          plot([beta1 beta2],'-o')
                          disp('iv2')
                          sum(iv2)
                          pause
                       end
            else
                
                if skewedges == 0
                
                    beta = ( sin(ny*(pi + thetaS + thetaout  ))./(ch - cos(ny*(pi + thetaS + thetaout   ))) + ...
                             sin(ny*(pi + thetaS - thetaout  ))./(ch - cos(ny*(pi + thetaS - thetaout   ))) + ...
                             sin(ny*(pi - thetaS + thetaout  ))./(ch - cos(ny*(pi - thetaS + thetaout   ))) + ...
                             sin(ny*(pi - thetaS - thetaout  ))./(ch - cos(ny*(pi - thetaS - thetaout   ))));           
                         beta1 = sin(ny*(pi + thetaS + thetaout  ))./(ch - cos(ny*(pi + thetaS + thetaout   )));
                         beta2 = sin(ny*(pi + thetaS - thetaout  ))./(ch - cos(ny*(pi + thetaS - thetaout   )));
                         beta3 = sin(ny*(pi - thetaS + thetaout  ))./(ch - cos(ny*(pi - thetaS + thetaout   )));
                         beta4 = sin(ny*(pi - thetaS - thetaout  ))./(ch - cos(ny*(pi - thetaS - thetaout   )));
                      iv1 = abs(beta1)>100;
                       iv2 = abs(beta2)>100;
                      iv3 = abs(beta3)>100;
                       iv4 = abs(beta3)>100;
                                                    
                      if sum(iv1) > 0                       
                          plot([beta1 beta2 beta3 beta4],'-o')
                          disp('iv1')
                          sum(iv1)
                          pause
                       end
                       if sum(iv2) > 0                       
                          plot([beta1 beta2 beta3 beta4],'-o')
                          disp('iv2')
                          sum(iv2)
                          pause
                       end
                      if sum(iv3) > 0                       
                          plot([beta1 beta2 beta3 beta4],'-o')
                          disp('iv3')
                          sum(iv3)
                          pause
                       end
                       if sum(iv2) > 0                       
                          plot([beta1 beta2 beta3 beta4],'-o')
                          disp('iv4')
                          sum(iv4)
                          pause
                       end
                         
                else
                    
                     beta = ( sin(ny*(pi + thetaS + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetaS + thetae3_re2(n3vertmat)   ))) + ...
                     sin(ny*(pi + thetaS - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetaS - thetae3_re2(n3vertmat)   ))) + ...
                     sin(ny*(pi - thetaS + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetaS + thetae3_re2(n3vertmat)   ))) + ...
                     sin(ny*(pi - thetaS - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetaS - thetae3_re2(n3vertmat)   )))); 
                    
                end
            end
                        
          Gsub = -(2-acrossface_out)*ny/4/pi*beta.*dz_expjkmoverm;
          
        G(Hsubmatrixdata.bigmatrixstartnums(ii):Hsubmatrixdata.bigmatrixendnums(ii)) =     G(Hsubmatrixdata.bigmatrixstartnums(ii):Hsubmatrixdata.bigmatrixendnums(ii)) + Gsub;
                
    end
    
end

Qs = 1;
Qfirstterm = G*Qs;
