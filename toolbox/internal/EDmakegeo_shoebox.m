function [corners,planecorners] = EDmakegeo_shoebox(length_x,length_y,length_z)
% This function generates the corners and planecorners matrices on the form
% that is used in the EDtoolbox for an external shoebox. By default, the
% corners are symmetrically distributed in the x- and y-directions, but
% placed at z = 0 and z = -length_z.
%
% Input parameters:
%   length_x,length_y, length_z     The sizes in meters
%
% Output parameters:
%   corners             A matrix of size [8,3] giving the x-,y- and
%                       z-coordinates of the 8 corners.
%   planecorners        A matrix of size [6,4] giving the corners for each
%                       of the 6 planes.
% 
% Peter Svensson (peter.svensson@ntnu.no) 6 Oct. 2023
% 
% [corners,planecorners] = EDmakegeo_shoebox(length_x,length_y,length_z);

% 28 Aug. 2020 First version
% 4 Oct. 2023 Copied into the EDtoolbox
% 6 Oct. 2023 Changed name from EDmakeshoeboxgeo to EDmakegeo_shoebox

corners = [     -length_x/2   -length_y/2   0
    length_x/2   -length_y/2   0
    length_x/2    length_y/2   0
   -length_x/2    length_y/2   0
     -length_x/2   -length_y/2   -length_z
    length_x/2   -length_y/2   -length_z
    length_x/2    length_y/2   -length_z
   -length_x/2    length_y/2   -length_z   ];

planecorners = [  1 2 3 4
    1 5 6 2
    2 6 7 3 
    4 3 7 8
    1 4 8 5
    8 7 6 5];



