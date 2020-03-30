function radians = counts2radians(counts)
% COUNTS2RADIANS converts encouder counts to radians.
%   radians = COUNTS2RADIANS(counts)

%radiansPERcount = 1/16; % THIS CONVERSION FACTOR NEEDS TO BE CHANGED!
radiansPERcount = 0.001963494997279; % 15 decimal places... CRAZY!
%radiansPERcount = 0.002;
radians = radiansPERcount*counts;