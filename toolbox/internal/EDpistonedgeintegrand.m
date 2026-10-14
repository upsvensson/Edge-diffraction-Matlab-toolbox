function outvalue = EDpistonedgeintegrand(x,psi2,zR2,k)
% EDpistonedgeintegrand - The integrand for the direct sound from a polygon piston.
% The direct sound from a piston is computed as a line integral along each
% piston edge (by EDmakefirstordertfs). This is the integrand
%   exp(-jkR)/(psi^2 + x^2),  R = sqrt(x^2 + psi^2 + zR^2)
%
% Input parameters:
%   x               Positions along the piston edge, in m, relative to the
%                   point on the edge (line) closest to the projection of
%                   the receiver onto the piston plane
%   psi2            The square of the distance from that projection to the
%                   edge (line)
%   zR2             The square of the distance from the receiver to the
%                   piston plane
%   k               The wavenumber
%
% Output parameters:
%   outvalue        The integrand at each x
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% outvalue = EDpistonedgeintegrand(x,psi2,zR2,k);

% R = sqrt( x.^2 + psi2 + zR2);
% outvalue = exp(-1i*k*R)./(psi2 + x.^2);
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

outvalue = exp(-1i*k*sqrt( x.^2 + psi2 + zR2))./(psi2 + x.^2);
