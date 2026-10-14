function [xvec, wvec] = EDgetGLrule(n)
% EDgetGLrule returns the GL nodes and weights, using lgwt from Matlab
% Central. The values are cached which means that lgwt will not be called
% the second time this function is called with the same quad. order n.
%
% Input parameters:
%   n       The quadrature order
%
% Output parameters:
%   xvec    The array of nodes, in [0,1].
%   wvec    The array of weight values belonging to the node values.
%
% Peter Svensson (peter.svensson@ntnu.no) 15 Oct. 2026
%
% [xvec, wvec] = EDgetGLrule(n);

% 25 Sep. 2026 Firstversion, suggested by Claude
% 15 Oct. 2026 The header was updated (by Claude, prompted by PS).

persistent GLrulecache
if isempty(GLrulecache), GLrulecache = {}; end
if n > numel(GLrulecache) || isempty(GLrulecache{n})
    [xvec, wvec] = EDdistelements(n,2);
    GLrulecache{n} = [xvec, wvec];
end
xvec = GLrulecache{n}(:,1);
wvec = GLrulecache{n}(:,2);
