%% SCRIPT_PeppersExample
% This example corresponds to the document:
%   PART 2 - Image Processing.docx

%% Get "Peppers" image
im = imread('Peppers.png');

%% Creating a Color Thresholding Function
if exist('createMask','file') ~= 2
    colorThresholder(im);
end

%% Color Thresholding
bin = createMask(im);

%% Filling Holes
binFill = imfill(bin,'holes');

%% Minimum Area Thresholding
minArea = 1000;
binMinArea = bwareaopen(binFill,minArea);

%% Labelling Connected Components
[lbl, n] = bwlabel(binMinArea);

%% Calculating Labelled Region Properties
stats = regionprops(lbl,'area','centroid','eccentricity');

%% Region Property Thresholding
minArea = 2000;
maxArea = 10000;
minEcc = 0.90;
maxEcc = 0.95;

idx = find(...
	[stats.Area] >= minArea & ...
	[stats.Area] <= maxArea & ... 
	[stats.Eccentricity] >= minEcc & ...
	[stats.Eccentricity] <= maxEcc);

statsTargets = stats(idx);

for i = idx
	binTargets{i} = lbl == i;
end

%% Plotting results
% Create the figure and axes
fig = figure('Parent',0);
axs = axes('Parent',fig,'NextPlot','add','XLim',[0,size(im,2)],'YLim',[0,size(im,1)]);
% Create the image object
img = imshow(im,'Parent',axs);
% Create a line object to display the binary target
thl = plot(0,0,'.c','Parent',axs);
% Create a line object to display the centroid
cnt = plot(0,0,'*m','Parent',axs,'MarkerSize',10);

% Find the x and y pixel locations for the binary of the 1st viable target
[y,x] = find(binTargets{1} == 1);
% Display the pixel locations
set(thl,'XData',x,'YData',y);

% Display the centroid location of the 1st viable target
set(cnt,'XData',statsTargets(1).Centroid(1),'YData',statsTargets(1).Centroid(2));

% Update the plot
drawnow;