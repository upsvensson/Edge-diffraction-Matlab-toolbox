function [nvec,weightvec] = EDdistelements(n,method)
% EDdistelements returns two vectors: fractional positions of n elements
% ditributed from 0 to 1, and the corresponding weightvec.
%
% Input parameters:
%   n       The number of elements to distribute
%   method  0 = linear distribution
%           1 = 1,2,4,8 towards middle
%           2 = Gauss-legendre
%
% Output parameters:
%   nvec        Vector of fractional positions between 0 and 1
%   weightvec   Vector of weighting values for the sample points in nvec. 
%
% Uses function quad_gauss
%
% peter.svensson@ntnu.no       27 Nov. 2017
%
% [nvec,weightvec] = EDdistelements(n,method);

% 10 March 2011 First version
% 15 March 2012 Added Gauss-Legendre
% 27 Nov. 2017 Copied from ESIE2toolbox

if method == 0
    nvec = ([0:n-1].'+0.5)/n;
    weightvec = ones(n,1)/n;    
elseif method == 1
    sizevec = 2.^([0:floor((n-1)/2)].');

    if round(n/2)*2 == n
        sizevec = [sizevec sizevec(end:-1:1)];
    else
        sizevec = [sizevec sizevec(end-1:-1:1)];    
    end

    weightvec = sizevec/sum(sizevec);

    cumn = [0 cumsum(sizevec)];
    nvec = cumn(1:end-1).' + sizevec/2;
    nvec = nvec/sum(sizevec);
    
elseif method == 2
    if n == 10
        nvec = [   0.013046735741414
           0.067468316655508
           0.160295215850488
           0.283302302935376
           0.425562830509184;zeros(5,1)];
        nvec = [nvec(1:5);1-nvec(5:-1:1)];
        weightvec = [   0.033335672154344
           0.074725674575290
           0.109543181257991
           0.134633359654998
           0.147762112357376;zeros(5,1)];
        weightvec = [weightvec(1:5);weightvec(5:-1:1)];
    elseif n == 18
        nvec = [   0.004217415789535
           0.022088025214301
           0.053698766751222
           0.098147520513738
           0.154156478469823
           0.220114584463026
           0.294124419268579
           0.374056887154247
           0.457612493479132;zeros(9,1)];
        nvec = [nvec(1:9);1-nvec(9:-1:1)];
        weightvec = [  0.010808006763242
           0.024857274447485
           0.038212865127444
           0.050471022053144
           0.061277603355739
           0.070321457335325
           0.077342337563133
           0.082138241872916
           0.084571191481572];
        weightvec = [weightvec(1:9);weightvec(9:-1:1)];

    elseif n == 20
        nvec = [   0.003435700407453
           0.018014036361043
           0.043882785874337
           0.080441514088891
           0.126834046769925
           0.181973159636742
           0.244566499024586
           0.313146955642290
           0.386107074429177
           0.461736739433251;zeros(10,1)];
        nvec = [nvec(1:10);1-nvec(10:-1:1)];
        weightvec = [   0.008807003569576
           0.020300714900194
           0.031336024167054
           0.041638370788352
           0.050965059908620
           0.059097265980759
           0.065844319224588
           0.071048054659191
           0.074586493236302
           0.076376693565363;zeros(10,1)];
        weightvec = [weightvec(1:10);weightvec(10:-1:1)];
    elseif n == 40
        nvec = [   0.000881145144720
           0.004636880650271
           0.011370025008113
           0.021041590393104
           0.033593595860662
           0.048950596515563
           0.067020248393870
           0.087693884583344
           0.110847174286741
           0.136340872405037
           0.164021657692910
           0.193723055166010
           0.225266437452436
           0.258462099156911
           0.293110397814198
           0.329002954587121
           0.365923907496373
           0.403651209649314
           0.441957964662372
           0.480613791246975;zeros(20,1)];
        nvec = [nvec(1:20);1-nvec(20:-1:1)];
        weightvec = [
           0.002260638549266
           0.005249142265576
           0.008210529190954
           0.011122924597084
           0.013968503490012
           0.016730097641274
           0.019391083987236
           0.021935454092837
           0.024347903817536
           0.026613923491968
           0.028719884549696
           0.030653121246465
           0.032402006728300
           0.033956022907617
           0.035305823695644
           0.036443291197902
           0.037361584528984
           0.038055180950313
           0.038519909082124
           0.038752973989213;zeros(20,1)];
        weightvec = [weightvec(1:20);weightvec(20:-1:1)];
    else        
        [nvec,weightvec] = quad_gauss(n,0,1);
        nvec = nvec(:);
        weightvec = weightvec(:);
    end

end


