function [p] = prob_nhits(n,N,Ph)
% PROB_NHITS computes probability of exactly n hits in N shots. 
% 
%   p = PROB_NHITS(n,N,Ph) finds the probability of exactly n hits
%   in N shots, assuming the probability any individual shot will hit is
%   Ph. Ph is a probability in the range [0,1]; N and n should be positive
%   integers. The value returned should be in the range [0,1]. The value
%   returned should be in the range [0,1]. 
%
%   D Evangelista (2020)
%   evangeli at usna dot edu

p = nchoosek(N,n)*Ph^n*(1-Ph)^(N-n);
end
