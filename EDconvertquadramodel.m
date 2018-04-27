function [cornersout,planecornersout] = EDconvertquadramodel(cornersin,planecornersin)
% EDconvertquadramodel modifies the two matrices cornersin and
% planecornersin to give a proper EDmodel. Duplicates of corners
% definitions will be removed, and triangular planes will be written
% properly.
% 
% Input parameters:
%   cornersin       Matrix, [ncornersin,3]
%   planecornersin  Matrix, [nplanes,4]
% 
% Output parameters:
%   cornersout      Matrix, [ncornersout,4], where ncornersout < ncornersin
%                   if there were duplicates in ncornersin.
%   planecornersout Matrix, [nplanes,4]
% 
% Peter Svensson 27 Apr 2018 (peter.svensson@ntnu.no)
% 
% [cornersout,planecornersout] = EDconvertquadramodel(cornersin,planecornersin);

% 27 Apr 2018 First version


ncorners = size(corners1,1);
nplanes = size(planecorners1,1);

corners2 = [[1:ncorners].' corners1];
planecorners2 = planecorners1;
planecorners2 = reshape(planecorners2,nplanes*4,1);

% Remove corners with identical coordinates, and replace the removed
% cornernumbers in planecorners (and replace by the single "keeper").

cornercounter = 0;
while cornercounter+1 < size(corners2,1)
    cornercounter = cornercounter + 1;
    c1 = corners2(cornercounter,2:4); 
    
    distvec = EDcalcdist(c1,corners2(:,2:4));
    iv = find(abs(distvec)<1e-10);
    if length(iv) > 1
        conumber_keep    = corners2(cornercounter,1);
        conumber_replace = corners2(iv(2:end),1);
        
        for jj = 1:length(conumber_replace)
           iv2 = find(corners2(:,1) == conumber_replace(jj));
           for kk = 1:length(iv2)
              corners2(iv2(kk),:) = []; 
           end
           iv2 = find(planecorners2 == conumber_replace(jj));
           planecorners2(iv2) = conumber_keep;
        end
        
    end
end

planecorners2 = reshape(planecorners2,nplanes,4);

% Go through the planes and detect triangle definitions

for ii = 1:nplanes
    colist = [planecorners2(ii,:) planecorners2(ii,1)];
   dco = diff(colist);
   iv = find(dco==0);
   if ~isempty(iv)
        
      if iv(1) == 1
         planecorners2(ii,:) = [colist(2:4) 0]; 
      end
      if iv(1) == 3
          planecorners2(ii,:) = [colist(1:3) 0];
      end       
   end
   
end

