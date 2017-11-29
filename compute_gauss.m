function [x,w] = compute_gauss(N, alpha, beta)
%function [x,w] = compute_gauss(N, alpha, beta)
%
%   Compute the Gaussian rule corresponding to the orthogonal polynomials
%   that are defined by the coefficients alpha and beta, according to the
%   three-term recurrence relation for orthogonal polynomials:
%
%   p[i](x) = (x-alpha[i]) * p[i-1](x) - beta[i]*p[i-2](x),    i=1..N,
%
%   with p[0](x) = 1, and p[-1](x) = 0.
%
%   In:
%       N (1x1): number of quadrature points
%       alpha, beta (N x 1): recurrence coefficients
%   Out:
%       x (N x 1): quadrature points
%       w (N x 1): quadrature weights
%
%   References:
%   W. Gautschi, ``Orthogonal Polynomials: Computation and Approximation'',
%       Clarendon Press, Oxford, 2004.

% $Id$

% Use OPQ routine from Gautschi
xw = gauss(N, [alpha beta]);
x = xw(:,1);
w = xw(:,2);
