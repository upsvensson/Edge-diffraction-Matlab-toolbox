function [corners,planecorners] = EDconvertquadramodel(corners_input,planecorners_input)
% EDconvertquadramodel modifies the two matrices cornersin and
% planecornersin to give a proper EDmodel. Duplicates of corners
% definitions will be removed, and triangular planes will be written
% properly.
% 
% Input parameters:
%   corners_input       Matrix, [ncornersin,3]
%   planecorners_input  Matrix, [nplanes,4]
% 
% Output parameters:
%   corners      Matrix, [ncornersout,4], where ncornersout < ncornersin
%                if there were duplicates in ncornersin.
%   planecorners Matrix, [nplanes,4]
% 
% Peter Svensson 27 Apr 2018 (peter.svensson@ntnu.no)
% 
% [corners,planecorners] = EDconvertquadramodel(corners_input,planecorners_input);

% 27 Apr 2018 First version


ncorners = size(corners_input,1);
nplanes = size(planecorners_input,1);

corners = corners_input;
planecorners = planecorners_input;
planecorners = reshape(planecorners,nplanes*4,1);

% Find corners with identical coordinates, and replace the duplicate
% cornernumbers in planecorners (and replace by the single "keeper").

cornercounter = 0;
while cornercounter+1 < size(corners,1)
    cornercounter = cornercounter + 1;
    c1 = corners(cornercounter,:); 
    
    distvec = EDcalcdist(c1,corners);
    iv = find(abs(distvec)<1e-10);
    if length(iv) > 1
        conumber_keep    = cornercounter;
        conumber_replace = iv(2:end);
        
        for jj = 1:length(conumber_replace)
%            iv2 = find(corners(:,1) == conumber_replace(jj));
%            for kk = 1:length(iv2)
%               corners(iv2(kk),:) = []; 
%            end
           iv2 = find(planecorners == conumber_replace(jj));
           planecorners(iv2) = conumber_keep;
        end
        
    end
end

planecorners = reshape(planecorners,nplanes,4);

% Go through the planes and detect triangle definitions

for ii = 1:nplanes
    colist = [planecorners(ii,:) planecorners(ii,1)];
   dco = diff(colist);
   iv = find(dco==0);
   if ~isempty(iv)
        
      if iv(1) == 1
         planecorners(ii,:) = [colist(2:4) 0]; 
      end
      if iv(1) == 3
          planecorners(ii,:) = [colist(1:3) 0];
      end       
   end
   
end

