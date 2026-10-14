function [y,too_close] = EDquadgk_sub1(t,A,B,integralparameters)

% Transform to weaken singularities at both ends: [a,b] -> [-1,1]
tt = 0.25*(B-A)*t.*(3 - t.^2) + 0.5*(B+A);
[y,too_close] = evalFun(tt);
if ~too_close
    y = 0.75*(B-A)*y.*(1 - t.^2);
end
