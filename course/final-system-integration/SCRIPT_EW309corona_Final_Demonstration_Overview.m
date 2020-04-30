%% SCRIPT_EW309corona_Final_Demonstration_Overview
% This script provides the code associated with:
%   "EW309corona Final Demonstration Overview.pptx'.
%
%   M. Kutzer, 17Apr2020, USNA
clear all
close all
clc

%% (0) Create status figure
% -> We will use this to show & save the image (optional)
fig = figure('Name','SCRIPT_EW309corona_FinalProject_Skeleton');
axs = axes('Parent',fig);

%% (1.1) Define Test Parameters
% -> Value(s) that need to be defined by the user:
% (a) Target
targetRange    = 170;           % STUDENT DEFINED VALUE IN CENTIMETERS
targetDiameter = 12.7;          % STUDENT DEFINED VALUE IN CENTIMETERS
% (b) Probability
pOne           = 0.95;          % STUDENT DEFINED VALUE BETWEEN 0 AND 1
% (c) Control
Kp  = 1.0;  % STUDENT DEFINED VALUE, Proportional Gain (Kp)
Ki  = 2.0;  % STUDENT DEFINED VALUE, Integral Gain (Ki)
Kd  = 3.0;  % STUDENT DEFINED VALUE, Derivative Gain (Kd)
tf  = 5;    % STUDENT DEFINED VALUE IN SECONDS, Turret Controller Stop Time 

%% (1.2) Calculate Target Statistics
% -> Value(s) that need to be calculated:
nShots = 10;    % Number of Shots to Fire
xBias = -2.0;   % Target X-bias (cm)
yBias = 1.5;    % Target Y-bias (cm)
% [xBias,yBias,nShots] = NerfGunStats(targetRange,targetDiameter/2,pOne);

%% (1.3) Create Simulated Target
targetSpecs = createTargetSpecs(targetDiameter,xBias,yBias);

%% (2) Get Initial Target Image
% -> Target image is generated at a random angle between -25deg & 25deg
im = getTargetImage(targetRange,[],targetSpecs);

% -> Show & save the image (optional)
img = imshow(im,'Parent',axs);
hold(axs,'on');
set(axs,'Visible','on');
drawnow
%saveas(fig,'(2) Initial Target Image.png','png');

%% (3.1) Locate Target in Pixels
% -> Value(s) that need to be calculated:
% xPixels = ???
% yPixels = ???
% [xPixels,yPixels] = STUDENTFUNCTION(im);

%% (3.2) Locate Target in Centimeters
% -> Value(s) that need to be calculated:
% x_cm = ???
% y_cm = ???
% [x_cm,y_cm] = STUDENTFUNCTION(targetRange,xPixels,yPixels);

%% (4) Calculate Desired Turret Angle
% -> Value(s) that need to be calculated:
thetaDesired = deg2rad(-10);    % Desired Angle (radians)

%% (5.1) Package Control Parameters
cParams.Kp = Kp;
cParams.Ki = Ki;
cParams.Kd = Kd;
cParams.despos = thetaDesired;

%% (5.2) Rotate the Turret
tEval = linspace(0,tf,100);
[SSE,ts,t,theta] = sendCmdtoDcMotor('closed',cParams,tEval);

%% (6) Get Updated Target Image
im = getTargetImageUpdate(theta(end));

% -> Show & save the image (optional)
set(img,'CData',im);
drawnow;
%saveas(fig,'(6) Updated Target Image.png','png');

%% (7) Repeat Steps (3) – (6) OPTIONAL


%% (8) Fire At Target
im = getShotPatternImage(nShots);

% -> Show & save the image (optional)
set(img,'CData',im);
drawnow;
%saveas(fig,'(8) Shots Fired At Target.png','png');

%% (9.1) Performance assesment
EW309coronaPerforanceEval;
%saveas(gcf,'(9.1) Performance assesment.png','png');