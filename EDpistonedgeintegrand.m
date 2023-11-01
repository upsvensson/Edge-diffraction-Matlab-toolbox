function outvalue = EDpistonedgeintegrand(x,psi2,xp,zR2,k)

% R = sqrt( (x-xp).^2 + psi2 + zR2);
% outvalue = exp(-1i*k*R)./(psi2 + x.^2);

outvalue = exp(-1i*k*sqrt( (x-xp).^2 + psi2 + zR2))./(psi2 + x.^2);
