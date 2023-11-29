function Q = EDquadstep(a,b,fa,fc,fb,tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec)
% EDquadstep - Recursive numerical integration routine for EDwedge1st_ir.
% Modified from Matlabs QUADSTEP function. The EDquadstep version
% solves all time steps at the same time.
%
% Uses no special function 
%
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson (peter.svensson@iet.ntnu.no) 28 Jan 2018
%
% EDquadstep(a,b,fa,fc,fb,tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec)

% 21 March 2004 Functioning version in EDB toolbox
% 29 Dec. 2016 Renamed to ESIE2 without any changes
% 28 Jan 2018 Renamed from ESIE2 without any changes

% Evaluate integrand twice in interior of subinterval [a,b].
c = (a + b)/2;

nslots = length(a);

% We use the matrix x for both the x-interval values and the corresponding
% y-values

x = [(a + c)/2, (c + b)/2];
x = reshape(x,2*nslots,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instead of using a function call, we plug in the entire function here!
%% x = betaoverml_all(x,rs,rr,zr,ny,sinnyfivec,cosnyfivec);

ml = sqrt( (x-zs).^2 + rs^2 ).*sqrt( (x-zr).^2 + rr^2 );
%------------------------------------------------
% We would have liked to use the following code:
%       y = (ml + (zvec-zs).*(zvec-zr))/(rs*rr);
%       expnyeta = (sqrt(y.^2-1) + y).^ny;
%       coshnyeta = 0.5*( expnyeta + 1./expnyeta);
% but to save three vectors, zvec will be re-used for
% y and expnyeta and coshnyeta

zvec = (ml + (x-zs).*(x-zr))/(rs*rr);
zvec = (real(sqrt(zvec.^2-1)) + zvec).^ny;
zvec = 0.5*( zvec + 1./zvec);
x =     sinnyfivec(1)./(zvec - cosnyfivec(1));
x = x + sinnyfivec(2)./(zvec - cosnyfivec(2));
x = x + sinnyfivec(3)./(zvec - cosnyfivec(3));
x = x + sinnyfivec(4)./(zvec - cosnyfivec(4));         

x = x./ml;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = reshape(x,nslots,2);

fd = x(:,1);
fe = x(:,2);
ftemp = fa + fb + 2*fc;

%-----------------------------------------
% We use the output vector Q for another purpose

tempvec = (b-a)/6;

% Three point Simpson's rule.
Q1 = tempvec.*(ftemp + 2*fc);

% Five point double Simpson's rule.
%Q2 = (h/12).*(fa + 4*fd + 2*fc + 4*fe + fb);

% We would have liked to use the following code:
% % % Q2 = tempvec/2.*(ftemp + 4*fd + 4*fe);
% % % % One step of Romberg extrapolation.
% % % Q = Q2 + (Q2 - Q1)/15;
% ...... but we save one vector by re-using ftemp
% for Q2

ftemp = tempvec/2.*(ftemp + 4*fd + 4*fe);
% ftemp = tempvec/2.*(ftemp + 4*x(:,1) + 4*x(:,2));
% One step of Romberg extrapolation.
Q = ftemp + (ftemp - Q1)/15;

%-----------------------------------------
% Check if all or some of the sample intervals 
% have reached the tolerance.

tempvec = find(abs(ftemp-Q)>tol);

% Check accuracy of integral over this subinterval.

if isempty(tempvec)
   return

% Subdivide into two subintervals.
else
   Qac = EDquadstep(a(tempvec),c(tempvec),fa(tempvec),fd(tempvec),fc(tempvec),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
   Qcb = EDquadstep(c(tempvec),b(tempvec),fc(tempvec),fe(tempvec),fb(tempvec),tol,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);
   Q(tempvec) = Qac + Qcb;
end
