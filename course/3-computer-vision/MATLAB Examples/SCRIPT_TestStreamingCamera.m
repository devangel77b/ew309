%% SCRIPT_TestStreamingCamera

%% Initialize the camera and define camera preview
[cam,prv] = initCamera;
drawnow;

%% Grab an image from the preview and hide the preview figure
im = get(prv,'CData');          % Get current image from the camera preview
% Find figure containing preview
mom = prv;
while true
    mom = get(mom,'Parent');
    switch lower( get(mom,'Type') )
        case 'figure'
            prvFig = mom;
            break;
    end
end
set(prvFig,'Visible','off');    % Hide the preview figure

%% Initialize figure, axes, and plot objects for live streaming

% Create a figure, and save the figure handle as "figObj"
figObj = figure('Name','Live Processing Camera Stream');
% Create an axes and plot an image, and save the image handle as "imgObj"
imgObj = imshow(im);
% Get the axes handle and save it as "axsObj"
axsObj = get(imgObj,'Parent');
% Adjust the axes settings
% -> Show the x/y axes and change the location of the x-axis
set(axsObj,'Visible','on','xAxisLocation','top');
% -> Add a label to the x and y axes
xlabel(axsObj,'X-Pixel Coordinates');
ylabel(axsObj,'Y-Pixel Coordinates');
% -> Adjust the position of the axis (i.e. it's location on the monitor)
pos = get(axsObj,'Position');
set(axsObj,'Position',pos - [0.00, 0.07, 0.00, 0.00]);
% -> Allow additional plots to be added to the axes
hold(axsObj,'on');

% Create plot objects to show segmentation, targets, and target locations
% -> Plot object for showing "imBin_Segmented" (segmented object)
plt_Segmented = plot(0,0,'.m');
% -> Plot object for showing "imBin_Targets" (candidate targets)
plt_Targets = plot(0,0,'.c');
% -> Plot of target location
plt_Centers = plot(0,0,'*r','Visible','off');      

% Plot crosshair across center of image
plt_hCrosshair = plot(size(im,2)/2*[1, 1],[0.5,size(im,1)+0.5],':c','LineWidth',1.5);
plt_vCrosshair = plot([0.5,size(im,2)+0.5],size(im,1)/2*[1, 1],':c','LineWidth',1.5);

% Define text objects for displaying 
%   - area and 
%   - target location (relative to image center)
txt_Area = text(0,0,'EMPTY','Visible','off');   % Text displaying target area (in pixels)
txt_Position = text(0,0,'[NaN,NaN]','BackgroundColor','w','FontSize',18,'VerticalAlignment','Top');

%% Run loop where you can process the camera feed
while true
    % Get image from camera
    im = get(prv,'CData');
    % Process image
    [TargetLocations,imBin_Segmented,imBin_Targets] = ProcessTargetExample(im);
    
    % Show the image
    set(imgObj,'CData',im);
    % Show segmentation
    [y,x] = find(imBin_Segmented);
    set(plt_Segmented,'xData',x,'yData',y);
    % Show segmented targets
    [y,x] = find(imBin_Targets);
    set(plt_Targets,'xData',x,'yData',y);
    
    if ~isempty(TargetLocations)
        % Only consider the first target
        if size(TargetLocations,1) > 1
            warning('Multiple targets found, only displaying the first.');
        end
        TargetLocation = TargetLocations(1,:);
        
        % Display target centroid
        set(plt_Centers,'xData',TargetLocation(1),'yData',TargetLocation(2));
        set(plt_Centers,'Visible','on');
        
        % Calculate area of target
        strct = regionprops(imBin_Targets, 'Area');
        Area = strct.Area;
        % Display target area
        set(txt_Area,'Position',[TargetLocation,0],'String',sprintf('%d',Area));
        set(txt_Area,'Visible','on');
        set(txt_Position,'String',sprintf('[%.2f, %.2f]',TargetLocation - [size(im,2)/2, size(im,1)/2]));
        set(txt_Position,'Visible','on');
    else
        % Hide the centroid and text if no target is found
        set(plt_Centers,'Visible','off');
        set(txt_Area,'Visible','off');
        set(txt_Position,'Visible','off');
    end
    drawnow
    
end
