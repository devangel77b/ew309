%% SCRIPT_ProcessBallisticsImages

% Plot the image
img = imshow(im);
axs = gca;
hold(axs,'on');
pltV = plot(axs,[0.5, imSize(2)-0.5],[imSize(1)/2,imSize(1)/2],':c','LineWidth',1.5);
pltH = plot(axs,[imSize(2)/2,imSize(2)/2],[0.5, imSize(1)-0.5],':c','LineWidth',1.5);

% Define the total number of points to select
n = 5;
for i = 1:n
    % Prompt the user
    title( sprintf('Select Point %d of %d',i,n) );
    % Get the coordinate
    [x(i),y(i)] = ginput(1);
    % Plot the selected point
    plot(axs,x(i),y(i),'om');
    % Output the coordinate to the command window
    fprintf('x (pixels): %3.1f | y (pixels): %3.1f\n',x(i),y(i));
end