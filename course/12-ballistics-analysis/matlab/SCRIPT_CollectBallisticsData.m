%% SCRIPT_CollectBallisticsData

%% Initialize Camera
[~,prv] = initWebcam;

%% Highlight center of the image
im = get(prv,'CData');
imSize = size(im);
prvAxs = get(prv,'Parent');
hold(prvAxs,'on');
prvV = plot(prvAxs,[0.5, imSize(2)-0.5],[imSize(1)/2,imSize(1)/2],':c','LineWidth',1.5);
prvH = plot(prvAxs,[imSize(2)/2,imSize(2)/2],[0.5, imSize(1)-0.5],':c','LineWidth',1.5);

%% Prompt user for test information
% Prompt for total number of tests
dlgPrompt = {'How many distances are you testing?'};
dlgTitle = 'Number of Tests';
dlgDims = [1,40];
dlgDefInput = {'3'};
answer = inputdlg(dlgPrompt,dlgTitle,dlgDims,dlgDefInput);

% Loop through the total number of tests
for i = 1:str2double(answer{1})
    % Prompt the user for test distance and units
    dlgPrompt = {'Enter Test Distance:','Enter Distance Units:'};
    title = 'Test Distance';
    dims = [1 35];
    definput = {'','inches'};
    answer = inputdlg(dlgPrompt,title,dims,definput);

    % Bring the preview to the front
    axes(prvAxs)
    drawnow
    
    % Remind the user to center the crosshair
    questdlg('Center the Crosshair in the Preview!', ...
        'Center Reminder','OK','OK');
    
    % Fire shots
    % PUT YOUR FIRE CODE HERE (BE SURE TO INITIALIZE AND OPEN YOUR SERIAL
    % PORT BEFORE THIS LOOP!
    
    % Ask user to mark the point of impact
    % Remind the user to center the crosshair
    questdlg('Clearly Mark the Points of Impact', ...
        'Impact Reminder','OK','OK');
    
    % Take image
    im = get(prv,'CData');
    fname = sprintf('BallisticsTest_%s%s.png',answer{1},answer{2});
    imwrite(im,fname);
end
