%% SCRIPT_FitTransferFunction_SingleStepResponse
% This script is an investigation using "iddata" and "tfest" functions in 
% to fit a specified transfer function to step response data collected from
% the EW309 turrets.
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

% Simulate data using step response applied to ideal transfer functions

% Applied PWM signal
PWM_val = 0.55;     % Define PWM applied for step response

% Ideal (model-based) data
ideal_t = [0:0.01:10]';  % Assume we are measuring for 10 seconds on at 100Hz
ideal_theta_dot = step(PWM_val*G_star,ideal_t); % Estimate perfect velocity data
ideal_theta = step(PWM_val*H_star,ideal_t);     % Estimate perfect position data

% Apply noise to simulate real measurements
t = ideal_t + ...                                              % "perfect" sampling time
    normrnd(zeros(size(ideal_t)),repmat(0.005,size(ideal_t))); % noise
theta_dot = ideal_theta_dot + ...                                                   % "perfect" velocity data
    normrnd(zeros(size(ideal_theta_dot)),repmat(b_star/200,size(ideal_theta_dot))); % noise
theta = ideal_theta + ...                                                   % "perfect" position data 
    normrnd(zeros(size(ideal_theta)),repmat(b_star/200,size(ideal_theta))); % noise

% Define PWM signal array
PWM = repmat(PWM_val ,size(theta));

% Add equilibrium values --------------------------------------------------
% -> The fitting tool "tfest" recommends at least some equilibrium values
%    when defining the fit. These are values for the system prior to 
%    applying the step response that we will prepend as zero.
% -> When collecting data from your turret, all you need to do is start
%    your code and send back data, wait a second or two, then set the PWM
%    to a fixed value. This will collect equilibrium data for you.

% Prepend equilibrium values to input and output signal
% -> Define equilibrium time values
dt = mean( diff(t) );
prepend_t = [0:dt:(max(t)/5)]';
% -> Define equilibrium input values
prepend_PWM = zeros(size(prepend_t));
% -> Define equilibrium output values
prepend_theta_dot = zeros(size(prepend_t));
prepend_theta = zeros(size(prepend_t));

% -> Prepend values
t = [prepend_t;...
    ( max(prepend_t) + dt ) + ( t - min(t) )];
theta_dot = [prepend_theta_dot; theta_dot];
theta = [prepend_theta; theta];
PWM = [prepend_PWM; PWM];
% -------------------------------------------------------------------------

%% Plot collected data
% Position step response plot
fig_pos = figure('Name','Position Step Response');
axs_pos = axes('Parent',fig_pos);
hold(axs_pos,'on');
plot(axs_pos,t,theta,'b');
xlabel(axs_pos,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_pos,'$\theta(t)$ (deg)','Interpreter','Latex');
title(axs_pos,'Position Step Response');

% Velocity step response plot
fig_vel = figure('Name','Velocity Step Response');
axs_vel = axes('Parent',fig_vel);
hold(axs_vel,'on');
plot(axs_vel,t,theta_dot,'b');
xlabel(axs_vel,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_vel,'$\dot{\theta}(t)$ (deg/s)','Interpreter','Latex');
title(axs_vel,'Velocity Step Response');

%% Fit transfer function to position step response
% NOTE: This is what you can use to try and fit your transfer function
%       using experimental data.

% Estimate data for uniform sample frequency ------------------------------
% -> The fitting tool "tfest" used to estimate the coefficients of the 
%    transfer function requires a uniform sample frequency.

% Create a uniform time vector
uniform_t = linspace(min(t),max(t),numel(t))';
% Use linear interpolation to estimate PWM and theta for the uniform time 
% vector 
uniform_PWM = interp1(t,PWM,uniform_t);
uniform_theta = interp1(t,theta,uniform_t);
% -------------------------------------------------------------------------

% Format collected position data
data_pos = iddata(...
    uniform_theta,...                   % Measured output  
    uniform_PWM,...                     % Applied input
    [],...                              % Sampling frequency (empty if you supply sample times)
    'SamplingInstants',uniform_t,...    % Sample times
    'InputName',{'Signed PWM'},...      % Specify input name(s)       (optional)
    'InputUnit',{'unitless'},...        % Specify units of input(s)   (optional)
    'OutputName',{'Turret Angle'},...   % Specify output name(s)      (optional)
    'OutputUnit',{'deg'});              % Specify units for output(s) (optional)

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
% Find where step is initiated
idx = find(uniform_PWM == 0, 1, 'Last');
idx_Start = idx+1;  % Index value where step input is applied
% Define fixed PWM value of step
PWM_val = mean(uniform_PWM(idx_Start:end));

% Position response
[theta_estimate,t_out] = step(PWM_val*H,...
    uniform_t(idx_Start:end) - uniform_t(idx_Start));
% Plot response
plot(axs_pos,t_out + uniform_t(idx_Start),theta_estimate,'g','LineWidth',2);

% Velocity response
[theta_dot_estimate,t_out]  = step(PWM_val*s*H,...
    uniform_t(idx_Start:end) - uniform_t(idx_Start));
% Plot response
plot(axs_vel,t_out + uniform_t(idx_Start),theta_dot_estimate,'g','LineWidth',2);

%% Display Transfer Function Coefficients
% Parse information from the identified transfer function object
pvec = getpvec(H,'free');

% Label the coefficients from the transfer function assuming the following
% form:
%                b
% H(s) = -------------------
%         s^2 + c_1*s + c_2
b = pvec(1);
c(1) = pvec(2);
c(2) = pvec(3);

% Calculate the poles using algebraic manipulation
syms x
fcn = x^2 + c(1)*x + c(2);
factored_fcn = factor(fcn, x, 'FactorMode', 'real');
for i = 1:numel(factored_fcn)
    coeffs_fcn = coeffs(factored_fcn(i));
    p(i) = coeffs_fcn(1);
end

% Label new coefficients assuming the form of the proposed transfer
% function:
%                  b
% H(s) = ---------------------
%         (s + delta)*(s + a)
delta = p(1);
a = p(2);

% Displace transfer function
fprintf('Estimated Transfer function:\n');
fprintf('\t H(s) = %f / ( (s + %f)*(s + %f) ) \n',b,delta,a);
