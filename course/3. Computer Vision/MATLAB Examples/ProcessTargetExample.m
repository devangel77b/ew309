function [TargetLocation,imBin_Seg,imBin_Tar] = ProcessTargetExample(im)
% PROCESSTARGET processes a color image given threshold values, minimum
% area, and eccentricity. 
%
%   TargetLocation = PROCESSTARGET(im) takes an RGB image "im" and returns
%   candidate target locations as an Nx2 array (where N indicates the total
%   number of targets.
%
%   [TargetLocation,imBin_Seg,imBin_Tar] = PROCESSTARGET(im) returns
%   candidate target locations, the initial image segmentation, and the
%   segmented "targets".
%
%   Unknown, --XXX----, USNA

% Updates
%   09Apr2017 - Revised to provide outputs instead of plots, M. Kutzer
%   13Mar2018 - Updated documentation, M. Kutzer
%   08Jan2019 - Updated documentation, M. Kutzer

%% Check inputs
narginchk(1,1);
if size(im,3) ~= 3
    error('Input must be a color image.');
end

%% Initialize target location
TargetLocation = []; % initialized to empty

%% Set thresholding values (YOU WILL WANT TO CHANGE THESE)
% -> Low and high red threshold values
R_low = 181;
R_high = 270;
% -> Low and high green threshold values
G_low = 60;
G_high = 162;
% -> Low and high blue threshold values
B_low = 32; 
B_high = 122;

%% Define minimum and maximum area in pixels (YOU WILL WANT TO CHANGE THESE)
MinArea = 200;
MaxArea = pi*(322)^2/4;

%% Define maximum eccentricity (YOU MAY WANT TO CHANGE THESE)
MaxEccentricity = 0.5;

%% Segment image 
imBin_Seg = im(:,:,1) >= R_low &  im(:,:,1) <= R_high & ...
            im(:,:,2) >= G_low &  im(:,:,2) <= G_high & ...
            im(:,:,3) >= B_low &  im(:,:,3) <= B_high;
    
%% Remove all segmented objects with an area less than the minimum
imBin_Seg = bwareaopen(imBin_Seg, MinArea);

%% Define connected regions and calculate region properties
cc = bwconncomp(imBin_Seg);
stats = regionprops(cc, 'Area', 'Centroid', 'Eccentricity');

%% Apply maximum eccentricity and area to define imBin2
idx = find([stats.Eccentricity]<MaxEccentricity & [stats.Area]<MaxArea );
imBin_Tar = ismember(labelmatrix(cc), idx);

%% Calculate areas and centroids
% Propeties are placed in an array
Areas = cat(1,stats(idx).Area); 
CX_CY = cat(1,stats(idx).Centroid);

%% Package target location(s) for output
TargetLocation = CX_CY;