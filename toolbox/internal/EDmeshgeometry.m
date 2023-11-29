function [apertureelements,apertureelementdistances] = EDmeshgeometry(Lx,Ly,nx,ny,geotype,calcdistances)
% EDgeometry creates two structs with geometry details for a
% discretization of an aperture.
% 
% Input parameters:
%   Lx, Ly      The size of the aperture, in meters
%   nx, ny      The number of discrete points, in the x- and y-directions
%   geotype     (optional) 1 -> uniform discretization with equally sized
%               elements. 2 -> Gauss-Legendre discretization. Default: 1
%   calcdistances (optional) 1 -> apertureelementdistances is generated
%               0 -> apertureelementdistances is not generated
%               Default: 1
% 
% Output parameters:
%   apertureelements    Struct with the following fields:
%   .midcoordsx, .midcoordsy    List, size [nx*ny,1], with the x- and
%                       y-coordinates of:
%                       (geotype 1) the center point of each element
%                       (geotype 2) the nodes for the quadrature
%   .weights            List, size [nx*ny,1], with:
%                       (geotype 1) the areas of all elements
%                       (geotype 2) the weights for the quadrature 
%   .nx, .ny            Same as the input parameters nx and ny
%   .Lx, .Ly            Same as the input parameters Lx and Ly
%   .dx, .dy            The size of each element (empty, if geotype = 2)
%   apertureelementdistances   Struct with the following fields:
% For geotype = 1
%   .fromelem, .toelem  Long lists, size [(nx*ny)^2,1], of all the source
%                       element numbers and all the receiver point numbers.
%   .shortlist_distances   A shortlist, size [nunique,1], of unique distances
%   .reftoshortlist_distances   A long list, size [(nx*ny)^2,1], which for
%                       each element-to-element combination points to the
%                       place in .shortlist_distances with the correct
%                       distance.
%   .examplefromshortlist_distances    A list, size [nunique,1], which for
%                       each entry in the shortlist gives one example of
%                       element-to-element combination with the correct
%                       distance.
% For geotype = 2
%   .alldistances       Matrix, size [nx*ny,nx*ny], of all node-to-node
%                       distances
% 
% Uses the function GEOdivrect, lgwt
%
% peter.svensson@ntnu.no 6 Dec. 2021
%
% [apertureelements,apertureelementdistances] = APcreategeometry3(Lx,Ly,nx,ny,geotype,calcdistances);

% 25 Oct. 2021 First version
%  2 Nov. 2021 Expansion to include Gauss-Legendre distribution of 
% discrete points
% 18 Nov. 2021 Split up of the output parameters into two structs. Also
% made it possible to skip the distances calculation.
% 5 Dec. 2021 Made a much "cheaper" calculation to find the unique
% distances. Assuming elements of uniform size, the number of elements
% shift in x- and y-direction is studied, instead of distances.
% 6 Dec. 2021 Calculates distance squared instead of distances, to find
% unique combinations.

if nargin < 5
    geotype = 1;
    calcdistances = 1;
end
if geotype == 2
   error(['ERROR: Sorry, geotype = 2 is not implemented yet']) 
end
if nargin < 6
    calcdistances = 1;
end

c1 = [-Lx/2 -Ly/2 0];
c2 = [Lx/2 -Ly/2 0];
c3 = [Lx/2 Ly/2 0];
c4 = [-Lx/2 Ly/2 0];

apertureelements = struct;
apertureelementdistances = struct;

apertureelements.nx = nx;
apertureelements.ny = ny;
apertureelements.Lx = Lx;
apertureelements.Ly = Ly;

ntot = apertureelements.nx*apertureelements.ny;

elementaspectratio = (Lx/nx/(Ly/ny));
if round(elementaspectratio*1e5) == 1e5
    elements_are_square = 1;
else
    elements_are_square = 0;    
end

if geotype == 1
    [patchcornercoords,patchmidcoords,aperture_areas] = GEOdivrect(c1,c2,c3,c4,nx,ny);
    apertureelements.midcoordsx = patchmidcoords(:,1);
    apertureelements.midcoordsy = patchmidcoords(:,2);
    apertureelements.areas = aperture_areas;
    apertureelements.dx = Lx/nx;
    apertureelements.dy = Ly/ny;

    if calcdistances == 1
        % Generate a long list (vertical vector) of all the "fromelem" = source
        % elements. Then generate an equally long list of all the "toelem" =
        % receiver elements. These lists have the size [nx*ny,nx*ny]

        fromelem = [1:ntot];
        if ntot < 65536
            fromelem = uint16(fromelem);
        else
            fromelem = uint32(fromelem);
        end
        fromelem = fromelem(ones(ntot,1),:);
        apertureelementdistances.fromelem = reshape(fromelem,ntot*ntot,1);
        clear fromelem
%         [Ifrom,Jfrom] = ind2sub([nx,ny],apertureelementdistances.fromelem);
%         fromrow = rem(apertureelementdistances.fromelem,ny);
%         iv = find(fromrow==0);
%         fromrow(iv) = ny;        
%         fromcol = ceil(double(apertureelementdistances.fromelem)/ny);
        
        toelem = [1:ntot].';
        if ntot < 65536
            toelem = uint16(toelem);
        else
            toelem = uint32(toelem);
        end
        toelem = toelem(:,ones(1,ntot));
        apertureelementdistances.toelem = reshape(toelem,ntot*ntot,1);
        clear toelem
%         [Ito,Jto] = ind2sub([nx,ny],apertureelementdistances.toelem);
%         torow = rem(apertureelementdistances.toelem,ny);
%         iv = find(torow==0);
%         torow(iv) = ny;        
%         tocol = ceil(double(apertureelementdistances.toelem)/ny);

%         Dnx = abs(Ito-Ifrom);
%         Dny = abs(Jto-Jfrom);
%         Dnx = abs(tocol-fromcol);
%         Dny = abs(torow-fromrow);
        
        % Calculate all the distances from the "fromelem" to the "toelem"
        % Then find all the unique distance values (within 5 decimals)
        % 
        % The list shortlist_distances will contain all the unique distance values.
        % The matrix reftoshortlist_distances will have a value for each of the
        % (nx*ny)^2 values in center_to_center_distances, pointing to the
        % "shortlist" where the much fewer unique values are stored.

        % Generate a short list of midelement-to-midelement distances
        % squared
        center_to_center_distances_squared = (apertureelements.midcoordsx(apertureelementdistances.toelem)-...
            apertureelements.midcoordsx(apertureelementdistances.fromelem)).^2 + ...
            (apertureelements.midcoordsy(apertureelementdistances.toelem)-...
            apertureelements.midcoordsy(apertureelementdistances.fromelem)).^2;  
%         
% % %         Dnmatrix = [Dnx Dny];
% % %         if elements_are_square == 1
% % %             Dnmatrix = sort(Dnmatrix.').';
% % %         end
% % %         [test,examplefromshortlist_distances,reftoshortlist_distances] = unique(Dnmatrix,'rows');
% % %                     
% % %         shortlist_distances = sqrt( (apertureelements.midcoordsx(...
% % %                     apertureelementdistances.toelem(examplefromshortlist_distances))-...
% % %                     apertureelements.midcoordsx(apertureelementdistances.fromelem(examplefromshortlist_distances))).^2 + ...
% % %             (apertureelements.midcoordsy(apertureelementdistances.toelem(examplefromshortlist_distances))-...
% % %             apertureelements.midcoordsy(apertureelementdistances.fromelem(examplefromshortlist_distances))).^2 );  
% % % 
        
         [shortlist_distances,examplefromshortlist_distances,reftoshortlist_distances] = unique(round(center_to_center_distances_squared*1e5)/1e5);
         clear center_to_center_distances
         
        % To save some memory, we convert the reftoshortlist_distances, which
        % contains only integers to 16-bit or 32-bit values.
        if length(shortlist_distances) < 65536
            reftoshortlist_distances = uint16(reftoshortlist_distances);
        else
            reftoshortlist_distances = uint32(reftoshortlist_distances);
        end 
        apertureelementdistances.shortlist_distances = sqrt(shortlist_distances);
        apertureelementdistances.reftoshortlist_distances = reftoshortlist_distances;
        apertureelementdistances.examplefromshortlist_distances = examplefromshortlist_distances;

        [sortedvec,isort] = sort(apertureelementdistances.reftoshortlist_distances);
        dsort = diff(sortedvec);
        istep = find(dsort);
        listofalloccurences = cell(length(istep)+1,1);
        for jj = 1:length(istep)+1
            if jj == 1
                istart = 1;
            else
                istart = istep(jj-1)+1;
            end
            if jj == length(istep)+1
               iend = length(sortedvec);
            else
               iend = istep(jj); 
            end

            listofalloccurences{jj} = uint32(isort(istart:iend));        
        end
        apertureelementdistances.listofalloccurences = listofalloccurences;
    end
    
else
    [xvec,weightx] = lgwt(nx,-Lx/2,Lx/2);
    xvec = xvec(end:-1:1);
    weightx = weightx(end:-1:1);
    [yvec,weighty] = lgwt(ny,-Ly/2,Ly/2);
    yvec = yvec(end:-1:1);
    weighty = weighty(end:-1:1);
    [XX,YY] = meshgrid(xvec,yvec);
    XX = reshape(XX,nx*ny,1);
    apertureelements.midcoordsx = XX;
    YY = reshape(YY,nx*ny,1);
    apertureelements.midcoordsy = YY;    
    [WWX,WWY] = meshgrid(weightx,weighty);
    apertureelements.weights = WWX.*WWY;    
    apertureelements.distances = EDcalcdist([XX YY zeros(size(XX))],[XX YY zeros(size(XX))+1e-6]);    

end
