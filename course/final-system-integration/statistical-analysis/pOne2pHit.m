function pHit = pOne2pHit(pOne,N)
% PONE2PHIT calculates the probability of a hit given the probability of a
% single hit with N shots taken.
%   pHit = PONE2PHIT(pOne,N) calculates the probability of a single shot
%   hitting the target given the probability of a single hit when N shots
%   are taken.
%
%   M. Kutzer, 17Apr2020, USNA

%% Check inputs 
narginchk(2,2);

%% Calculate pHit
% N = log(1-pOne)/log(1-pHit)
% log(1-pHit) = log(1-pOne)/N
% 1-pHit = exp( log(1-pOne)/N )
pHit = 1-exp( log(1-pOne)/N );