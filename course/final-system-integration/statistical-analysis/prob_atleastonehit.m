function [p] = prob_atleastonehit(N,Ph)
% PROB_ATLEASTONEHIT computes probability of at least 1 hit in N shots
% 
%   p = PROB_ATLEASTONEHIT(Ph,N) finds the probability of at least one hit
%   in N shots, assuming the probability any individual shot will hit is
%   Ph. Ph is a probability in the range [0,1]; N should be a positive
%   integer.
%
%   D Evangelista (2020)
%   evangeli at usna dot edu


% This is really 1-P(no hits), which is 1-(1-Ph)^N for an individual shot
% hit probability of Ph and N shots...
p = 1-(1-Ph)^N;
end
