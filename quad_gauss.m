function [x,w] = quad_gauss(N, a, b)
%function [x,w] = quad_gauss(N, a, b)
%
%   Return the Gauss-Legendre quadrature with N points in the interval
%   [a,b].
%
%   References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'',
%       Clarendon Press, Oxford, 2004.

% $Id$

% Use OPQ routine from Gautschi
ab = r_jacobi(N, 0, 0);
[x,w] = compute_gauss(N, ab(:,1), ab(:,2));

% scale [-1,1] to [a,b]
x = (b-a)*(x+1)/2 + a;
w = w*(b-a)/2;
