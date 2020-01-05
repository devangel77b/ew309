%% SCRIPT_FitTransferFunction_MultipleStepResponses
% This script is an investigation using "iddata" and "tfest" functions in 
% to fit a specified transfer function to step response data collected for
% multiple different signed PWM values from the EW309 turrets.
%
% NOTE: This example assumes that experimental data from the turret  will 
%      include:
%           PWM - signed PWM value (positive or negative) applied to 
%                 generate the step response of the system
%           t   - Nx1 array containing time data for each measured value of 
%                 theta (angular position of the turret) and theta_dot
%                 (angular velocity of the turret)
%           theta_dot - Nx1 array containing measured turret velocity data  
%           theta     - Nx1 array containing measured turret position data             
%
%   M. Kutzer, 11Apr2018, USNA

clear all
close all
clc

%% Define ideal transfer functions
% NOTE: We are using this section to create a model to *simulate* data to  
%       check if this method works. You do not need a simulated model if
%       you are collecting actual data from your actual system!
a_star = 3;
b_star = 457;
d_star = 0.0025;

s = tf('s');
% Velocity response transfer function G(s)
G_star = b_star/(s+a_star);
% Position response transfer function H(s)
H_star = b_star/((s+d_star)*(s+a_star));

fprintf('Simulation Transfer function:\n');
fprintf('\t H(s) = %f / ( (s + %f)(s + %f) ) \n',b_star,a_star,d_star);

%% Simulate data using step response applied to ideal transfer functions

% Applied PWM signals
% Define PWM applied for step response
PWM_ALL = [-1.00, -0.90, -0.80, -0.70, -0.60, 0.60, 0.70, 0.80, 0.90, 1.00];  

% Ideal (model-based) data
ideal_t = [0:0.01:10]';  % Assume we are measuring for 10 seconds on at 100Hz

for i = 1:numel(PWM_ALL)
    PWM = PWM_ALL(i);
    ideal_theta_dot = step(PWM*G_star,ideal_t); % Estimate perfect velocity data
    ideal_theta = step(PWM*H_star,ideal_t);     % Estimate perfect position data

    % Apply noise to simulate measurements
    t_ALL(:,i) = ideal_t + ... % "perfect" sampling time
        normrnd(zeros(size(ideal_t)),repmat(0.005,size(ideal_t))); % noise
    theta_dot_ALL(:,i) = ideal_theta_dot + ... % "perfect" velocity data
        normrnd(zeros(size(ideal_theta_dot)),repmat(b_star/200,size(ideal_theta_dot))); % noise
    theta_ALL(:,i) = ideal_theta + ...  % "perfect" position data 
        normrnd(zeros(size(ideal_theta)),repmat(b_star/200,size(ideal_theta))); % noise
end

%% Plot collected data
% Position step response plot
fig_pos = figure('Name','Position Step Response');
axs_pos = axes('Parent',fig_pos);
hold(axs_pos,'on');
xlabel(axs_pos,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_pos,'$\theta(t)$ (deg)','Interpreter','Latex');
title(axs_pos,'Position Step Response');
for i = 1:numel(PWM_ALL)
    color{i} = rand(1,3);
    plot(axs_pos,t_ALL(:,i),theta_ALL(:,i),'Color',color{i});
end

% Velocity step response plot
fig_vel = figure('Name','Velocity Step Response');
axs_vel = axes('Parent',fig_vel);
hold(axs_vel,'on');
xlabel(axs_vel,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_vel,'$\dot{\theta}(t)$ (deg/s)','Interpreter','Latex');
title(axs_vel,'Velocity Step Response');
for i = 1:numel(PWM_ALL)
    plot(axs_vel,t_ALL(:,i),theta_dot_ALL(:,i),'Color',color{i});
end

%% Fit transfer function to position step response
% NOTE: This is what you can use to try and fit your transfer function
%       using experimental data.

% Estimate data for uniform sample frequency ------------------------------
% -> The fitting tool "tfest" used to estimate the coefficients of the 
%    transfer function requires a uniform sample frequency.

% Create a uniform time vector
uniform_t = linspace(max(min(t_ALL)),min(max(t_ALL)),size(t_ALL,1))';
% Use linear interpolation to estimate theta for the uniform time vector
for i = 1:numel(PWM_ALL)
    uniform_theta_ALL(:,i) = interp1(t_ALL(:,i),theta_ALL(:,i),uniform_t);
end
% -------------------------------------------------------------------------

% Add equilibrium values --------------------------------------------------
% -> The fitting tool "tfest" recommends at least some equilibrium values
%    when defining the fit. These are values for the system prior to 
%    applying the step response that we will prepend as zero.

% Define input signal as an array
PWM_in = repmat(PWM_ALL,size(theta_ALL,1),1);

% Prepend equilibrium values to input and output signal
% -> Define equilibrium time values
dt = mean( diff(uniform_t) );
prepend_t = [0:dt:(max(uniform_t)/5)]';
% -> Define equilibrium input values
prepend_PWM_in = zeros(size(prepend_t,1),size(PWM_ALL,2));
% -> Define equilibrium output values
prepend_theta = zeros(size(prepend_t,1),size(theta_ALL,2));
% -> Prepend values
uniform_t = [prepend_t;...
             uniform_t - min(uniform_t) + dt + max(prepend_t)];
uniform_theta_ALL = [prepend_theta; uniform_theta_ALL];
PWM_in = [prepend_PWM_in; PWM_in];
% -------------------------------------------------------------------------

% Format collected position data
data_pos = iddata(...
    uniform_theta_ALL,...                % Measured output  
    PWM_in,...                           % Applied input
    [],...                               % Sampling frequency (empty if you supply sample times)
    'SamplingInstants',uniform_t);%,...  % Sample times
    %'InputName',{'Signed PWM'},...      % Specify input name(s)       (optional)
    %'InputUnit',{'unitless'},...        % Specify units of input(s)   (optional)
    %'OutputName',{'Turret Angle'},...   % Specify output name(s)      (optional)
    %'OutputUnit',{'deg'});              % Specify units for output(s) (optional)

% Fit trasfer function for position response
% Define number of poles
np = 2;
% Define number of zeros
nz = 0;
% Define transfer function fit options
opt = tfestOptions('InitMethod','all','InitialCondition','zero');
% Calculate plant transfer function
H = tfest(data_pos,np,nz,opt);

%% Plot fit response
t_estimate = linspace(min(uniform_t),max(uniform_t),numel(uniform_t));
for i = 1:numel(PWM_ALL)
    PWM = PWM_ALL(i);
    % Position response
    [theta_estimate,t_out] = step(PWM*H,t_estimate);
    plot(axs_pos,t_out,theta_estimate,'g','LineWidth',2);

    % Velocity response
    [theta_dot_estimate,t_out]  = step(PWM*s*H,t_estimate);
    plot(axs_vel,t_out,theta_dot_estimate,'g','LineWidth',2);
end

%% Display Transfer Function Coefficients
pvec = getpvec(H,'free');

a = pvec(2);
b = pvec(1);
d = pvec(3);

fprintf('Estimated Transfer function:\n');
fprintf('\t H(s) = %f / ( (s + %f)(s + %f) ) \n',b,a,d);
