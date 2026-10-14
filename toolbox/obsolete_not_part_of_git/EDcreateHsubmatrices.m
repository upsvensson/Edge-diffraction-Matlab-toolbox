function [Hsubdatalists,ivuselists,listofactiveHsubs,Hsubmatrixdata] = EDcreateHsubmatrices(Hsubmatrixdata,edgedata,edgetoedgedata,k,ngauss,newngauss,showtext)
% 
% Input parameters:
%
% Output parameters:
%
% Calls EDgetGLrule 
%
% Peter Svensson 9 Oct. 2026

% 26 Mar. 2026 Cut out from EDintegralequation_convex_tf
% 25 Sep. 2026 Changed to the new orientation of quadraturematrix, and
% calls EDgetGLrule rather than EDdistelements.
% Code improvement: changed nbig to nQ, bigmatrix... to Qvec...


    %
    % The number of sub-matrices is the same as the number of rows in
    % edgetripletlist (to-edge,via-edge,from-edge). Each submatrix 
    % contains transfer functions from the edge-source-pairs
    % (via-edge,from-edge) to the edge-source-pairs (to-edge,via-edge).
    % The submatrix thus multiplies with q-amplitudes in the rows given in
    % row jj
    % of bigmatrixstartnumbers and bigmatrixendnumbers and adds the result
    % to the q-amplitude vertical array in the rows given in row ii of 
    % bigmatrixstartnumbers and bigmatrixendnumbers. Row jj is the row where
    % edgepairlist has the values [via-edge,from-edge] and row ii is the row
    % where edgepairlist has the values [to-edge,via-edge].

    % In the for-loop, nn will be the Hsub number
    % The ii value will refer to the row-number in the edgepairlist for the 
    % receiving edge-pair.
    % The jj value will refer to the row-number in the edgepairlist for the 
    % sending edge-pair.


    
alltripletsthin = prod(double(Hsubmatrixdata.isthinplanetriplet));
notripletsthin = prod(1-double(Hsubmatrixdata.isthinplanetriplet));
 nQ = Hsubmatrixdata.Qvecendnums(end);
nsubmatrices = size(Hsubmatrixdata.listofsubmatrices,1);
 
listofactiveHsubs = zeros(size(Hsubmatrixdata.listofsubmatrices,1),1);
Hsubcounter = 1;

if alltripletsthin == 1
    thinplaneboostvec = 2*ones(nsubmatrices,1);
else
    if notripletsthin == 1
        thinplaneboostvec = ones(nsubmatrices,1);
    else
        thinplaneboostvec = double(Hsubmatrixdata.isthinplanetriplet(1:nsubmatrices)+1);
    end
end

Hsubdatalists = cell(nsubmatrices,1);
ivuselists = cell(nsubmatrices,1);

    n1prev = 0;
    n2prev = 0;
    n3prev = 0;
    shortlistprev = 0;
refto = edgetoedgedata.reftoshortlistE;
Hsubconstant = -1/4/pi/2;

for nn = 1:nsubmatrices

    edge1 = Hsubmatrixdata.edgetripletlist(nn,3);
    if newngauss(3) == ngauss
        n1 = Hsubmatrixdata.nedgeelems(edge1);
    else
        n1 = newngauss(3);
    end

    edge2 = Hsubmatrixdata.edgetripletlist(nn,2);                   
    if newngauss(2) == ngauss
        n2 = Hsubmatrixdata.nedgeelems(edge2);
    else
        n2 = newngauss(2);
    end
    len2 = edgedata.edgelengthvec(edge2);

    edge3 = Hsubmatrixdata.edgetripletlist(nn,1);
    if newngauss(1) == ngauss
        n3 = Hsubmatrixdata.nedgeelems(edge3);
    else
        n3 = newngauss(1);
%        [n1vec,weightvec] = EDdistelements(n3,2);
        [n1vec,weightvec] = EDgetGLrule(n3);
        % Hsubmatrixdata.quadraturematrix_pos(n3,1:n3)     = n1vec(:).';
        % Hsubmatrixdata.quadraturematrix_weights(n3,1:n3) = weightvec(:).';
        Hsubmatrixdata.quadraturematrix_pos(1:n3,3)     = n1vec(:);
        Hsubmatrixdata.quadraturematrix_weights(1:n3,n3) = weightvec(:);

    end
%             disp(['triplet no. ',int2str(nn),': n1,n2,n3 = ',int2str(n1),',',int2str(n2),',',int2str(n3)])
    shortlistnumber = Hsubmatrixdata.reftoshortlist(nn);

%             ind_start = 1 + (Hsubmatrixdata.Qvecstartnums(ii)-1) + (Hsubmatrixdata.Qvecstartnums(jj)-1)* nQ;

   if showtext >= 3
       disp(' ')
        disp(['         H-submatrix ',int2str(nn),': From edge ',int2str(edge1),' via edge ',int2str(edge2),' to edge ',int2str(edge3)])
        disp(['         of size ',int2str(n3*n2),' by ',int2str(n2*n1),'. Thinplaneboost = ',int2str(thinplaneboostvec(nn))])
        if symmetrycompression ==1
            disp(['         stored as shortlist number ',int2str(shortlistnumber)])
        end
   end

   % Build vertical [n3,n2] matrices and expand them horizontally
   % NB! n1 refers to edge 1, which is the source edge
   %     n3 refers to edge 3, which is the target edge
   % edgetriplets written in order [edge3 edge2 edge 1]

   if n1 ~= n1prev || n2 ~= n2prev || n3 ~= n3prev

%                if computer ~= 'MACI64'
       if ~strcmp(computer, 'MACI64')
           if n3 <= 255
               n3vertmat = uint8(1:n3);
           else
               n3vertmat = uint16(1:n3);                   
           end
           n3vertmat = reshape(repmat(n3vertmat,n2,1),n2*n3,1);           
           n3vertmat = n3vertmat(:,ones(1,n1*n2));

           if n2 <= 255
               n2vertmat = uint8(1:n2).';
           else
              n2vertmat = uint16(1:n2).'; 
           end
           n2vertmat = repmat(n2vertmat,n3,1);
           n2vertmat = n2vertmat(:,ones(1,n1*n2));

       % Build horizontal n2 and n1 matrices and expand them vertically
            if n1 <= 255                    
               n1hormat = uint8(1:n1);
            else
                n1hormat = uint16(1:n1);
            end
           n1hormat = repmat(n1hormat,1,n2);
               n1hormat = n1hormat(ones(n3*n2,1),:);

            if n2 <= 255
               n2hormat = uint8(1:n2);
            else
               n2hormat = uint16(1:n2); 
            end
               n2hormat = reshape(repmat(n2hormat,n1,1),1,n2*n1);           
           n2hormat = n2hormat(ones(n3*n2,1),:);

       % Now, we should use only the combinations where n2vertmat ==
       % n2hormat! Otherwise we get too many combinations: there should
       % be n1*n2*n3 unique combinations, but the Hsub-matrix has size
       % [n3*n2,n2*n1], which has n2 too many elements!
       %
       % Once we have the iv, which identifies the wanted matrix
       % elements, we can throw away n2hormat because it is identical
       % to n2vertmat in the wanted matrix elements!

           ivuse = find(n2hormat==n2vertmat);
           n3vertmat = n3vertmat(ivuse);
           n2vertmat = n2vertmat(ivuse);
           n1hormat  = n1hormat(ivuse);

       else
%                     [ivuse,n1hormat,n2vertmat,n3vertmat] = ESIE2calcedgeinteqmatrixsub2_mex(n1,n2,n3);       
            [ivuse,n1hormat,n2vertmat,n3vertmat] = EDcalcedgeinteqmatrixsub2_mex(n1,n2,n3);       

       end

       ivuse = uint64(ivuse);

   end

   % Build vectors that will be referred to by the n1,n2,n3
   % matrices

   if n1 ~= n1prev
        n1vec = full(Hsubmatrixdata.quadraturematrix_pos(n1,1:n1));
        n1vec = n1vec(:);
        weightvec = full(Hsubmatrixdata.quadraturematrix_weights(n1,1:n1));
        weightvec = weightvec(:);
        n1prev = n1;
   end       
   dzvec = weightvec*edgedata.edgelengthvec(edge1);

    if n2 ~= n2prev
        n2vec = full(Hsubmatrixdata.quadraturematrix_pos(n2,1:n2));
        n2vec = n2vec(:);
        n2prev = n2;
    end
    if n3 ~= n3prev        
        n3vec = full(Hsubmatrixdata.quadraturematrix_pos(n3,1:n3));
        n3vec = n3vec(:);
        n3prev = n3;
    end

    % Bug found and fixed 22 Sep. 2014: lines below did not always detect if another
    % edge was skewed relative to the first one, because only the
    % thetae1sho was used. Necessary to check thetae1sho and thetae2sho.

   ny = pi/(2*pi-edgedata.closwedangvec(edge2));
   thetain  = edgetoedgedata.thetae1sho(refto(edge1,edge2));
   thetain_other_end_of_edge  = edgetoedgedata.thetae2sho(refto(edge1,edge2));
%    if abs(thetain) < 1e-10 || abs(thetain-(2*pi-closwedangvec(edge2))) < 1e-10   % Old erroneous line
   if (abs(thetain) < 1e-10 || abs(thetain-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetain-thetain_other_end_of_edge) < 1e-10
       acrossface_in = 1;       
   else
       acrossface_in = 0;
   end
   thetaout = edgetoedgedata.thetae1sho(refto(edge3,edge2));
   thetaout_other_end_of_edge = edgetoedgedata.thetae2sho(refto(edge3,edge2));
%   if abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-closwedangvec(edge2))) < 1e-10    % Old erroneous line
   if (abs(thetaout) < 1e-10 || abs(thetaout-(2*pi-edgedata.closwedangvec(edge2))) < 1e-10) && abs(thetaout-thetaout_other_end_of_edge) < 1e-10
       acrossface_out = 1;
   else
       acrossface_out = 0;
   end

   if showtext >= 3
       disp(['         Acrossface-out = ',int2str(acrossface_out),' and Acrossface-in = ',int2str(acrossface_in), ' (not used)'])
   end

   ze2_re2 = n2vec*len2;
   ze3_re2 = edgetoedgedata.ze1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.ze2sho(refto(edge3,edge2))-edgetoedgedata.ze1sho(refto(edge3,edge2))   );
%        re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

   if shortlistnumber ~= shortlistprev      

       n1vec = n1vec(:).';

       ze1_re2 = edgetoedgedata.ze1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.ze2sho(refto(edge1,edge2))-edgetoedgedata.ze1sho(refto(edge1,edge2))   );

       % Bug fixed 26 Oct. 2017: the line below is wrong because the
       % r-values are not distributed linearly.
%            re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );


       if acrossface_in*acrossface_out == 0 % This case is a preparation for non-convex objects.
            % Reconvert the cylindrical coordinates of
            % edge1-re2 into cartesian! Then we can use the EDcoordtrans.
            % We define our own cartesian coord syst such that the reference
            % edge (edge2) has its starting point in [0 0 0], the z-axis along its
            % own edge, and the x-axis is along (onto) the reference plane.

            xedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
            xedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
            yedge1start = edgetoedgedata.re1sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge1,edge2)));
            yedge1end   = edgetoedgedata.re2sho(refto(edge1,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge1,edge2)));
            xe1_re2 = xedge1start + n1vec.*(xedge1end-xedge1start);
            ye1_re2 = yedge1start + n1vec.*(yedge1end-yedge1start);

            [re1_re2,thetae1_re2,~] = EDcoordtrans1([xe1_re2.' ye1_re2.' ze1_re2.'],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
            re1_re2 = re1_re2(:).';
            thetae1_re2 = thetae1_re2(:).';

            xedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
            xedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*cos(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
            yedge3start = edgetoedgedata.re1sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae1sho(refto(edge3,edge2)));
            yedge3end   = edgetoedgedata.re2sho(refto(edge3,edge2))*sin(edgetoedgedata.thetae2sho(refto(edge3,edge2)));
            xe3_re2 = xedge3start + n3vec.*(xedge3end-xedge3start);
            ye3_re2 = yedge3start + n3vec.*(yedge3end-yedge3start);

            [re3_re2,thetae3_re2,~] = EDcoordtrans1([xe3_re2 ye3_re2 ze3_re2],[0 0 0;0 0 edgedata.edgelengthvec(edge2)],[0 1 0]);
            re3_re2 = re3_re2(:);
            thetae3_re2 = thetae3_re2(:);

       else

            re1_re2 = edgetoedgedata.re1sho(refto(edge1,edge2)) + n1vec*( edgetoedgedata.re2sho(refto(edge1,edge2))-edgetoedgedata.re1sho(refto(edge1,edge2))   );
            re3_re2 = edgetoedgedata.re1sho(refto(edge3,edge2)) + n3vec*( edgetoedgedata.re2sho(refto(edge3,edge2))-edgetoedgedata.re1sho(refto(edge3,edge2))   );

       end

       ldist = (ze3_re2(n3vertmat) - ze2_re2(n2vertmat)).^2   + re3_re2(n3vertmat).^2;
       ldist = sqrt(ldist);

       mdist = (ze1_re2(n1hormat).' - ze2_re2(n2vertmat)).^2   + (re1_re2(n1hormat).').^2;
       mdist = sqrt(mdist);

       dz_expjkmoverm = exp(-1i*k*mdist)./mdist.*dzvec(n1hormat);           

       ch = ((ze2_re2(n2vertmat) - ze1_re2(n1hormat).').*(ze2_re2(n2vertmat) - ze3_re2(n3vertmat)) + mdist.*ldist)./re1_re2(n1hormat).'./re3_re2(n3vertmat);
       ch = ( real(sqrt( ch.^2-1)) + ch ).^ny;
       ch = ( ch + 1./ch)/2;

       % PS 28 March 2011: The 4 beta terms are identical!       
       if acrossface_in == 1 && acrossface_out == 1
           beta = 4*( sin(ny*(pi + thetain + thetaout  ))./(ch - cos(ny*(pi + thetain + thetaout   ))));
       else
           beta = ( sin(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                  ( sin(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi + thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                  ( sin(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' + thetae3_re2(n3vertmat)   ))+eps*10)) + ...
                  ( sin(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)  ))./(ch - cos(ny*(pi - thetae1_re2(n1hormat).' - thetae3_re2(n3vertmat)   ))+eps*10));
       end

       if ismember(edge1,edgedata.offedges)==1 || ismember(edge2,edgedata.offedges)==1 || ismember(edge3,edgedata.offedges)==1
           beta = beta*0;
       else
%           disp(['Edge1: ',int2str(edge1),' edge2: ',int2str(edge2),' edge3: ',int2str(edge3)])
       end

        Hsubdatalists{nn} = (2-acrossface_out)*Hsubconstant*ny*beta.*dz_expjkmoverm*thinplaneboostvec(nn);   
        ivuselists{nn} = ivuse;
%        if ifreq == 1
            listofactiveHsubs(Hsubcounter) = nn;
            Hsubcounter = Hsubcounter + 1;
 %       end

   end

    shortlistprev = shortlistnumber;

end
