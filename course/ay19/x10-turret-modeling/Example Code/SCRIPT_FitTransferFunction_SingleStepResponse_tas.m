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
PWM = 0.55;     % Define PWM applied for step response

% Ideal (model-based) data
ideal_t = [0:0.01:10]';  % Assume we are measuring for 10 seconds on at 100Hz
ideal_theta_dot = step(PWM*G_star,ideal_t); % Estimate perfect velocity data
ideal_theta = step(PWM*H_star,ideal_t);     % Estimate perfect position data

% Apply noise to simulate measurements
t = ideal_t + ... % "perfect" sampling time
    normrnd(zeros(size(ideal_t)),repmat(0.005,size(ideal_t))); % noise
theta_dot = ideal_theta_dot + ... % "perfect" velocity data
    normrnd(zeros(size(ideal_theta_dot)),repmat(b_star/200,size(ideal_theta_dot))); % noise
theta = ideal_theta + ...  % "perfect" position data 
    normrnd(zeros(size(ideal_theta)),repmat(b_star/200,size(ideal_theta))); % noise

%% Plot collected data

% Velocity step response plot
figure(1)
fig_vel = figure('Name','Velocity Step Response');
axs_vel = axes('Parent',fig_vel);
hold(axs_vel,'on');
plot(axs_vel,t,theta_dot,'b');
xlabel(axs_vel,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_vel,'$\dot{\theta}(t)$ (deg/s)','Interpreter','Latex');
title(axs_vel,'Velocity Step Response');
hold on

%If you want to fit the velocity response using fittype
fittype_measured = fittype('((0.5*b)/a)*(1-exp(-a*x))', 'dependent', {'y'}, 'independent', {'x'},'coefficients',{'a','b'});
fit_measured = fit(t,theta_dot,fittype_measured,'StartPoint',[1,1],'Lower',[0,0]);
%this will get you an estimate of a and b
co_effs = coeffvalues(fit_measured);
a_est = co_effs(1);
b_est = co_effs(2);
theta_dot_est = 0.5*b_est/a_est*(1-exp(-a_est*t));
plot(t,theta_dot_est,'r')

% Position step response plot
figure(2)
fig_pos = figure('Name','Position Step Response');
axs_pos = axes('Parent',fig_pos);
hold(axs_pos,'on');
plot(axs_pos,t,theta,'b');
xlabel(axs_pos,'$t$ (sec)','Interpreter','Latex');
ylabel(axs_pos,'$\theta(t)$ (deg)','Interpreter','Latex');
title(axs_pos,'Position Step Response');
hold on

%now try and match with delta value that moves the pole away from zero due
%to friction
d_est = 0.001;
%theta_est = 0.5*b/(a*delta)*(1 - exp(-delta*t)/(a-delta) + delta*exp(-a*t)/(a*(a-delta)));
theta_est = b_est/(2*a_est*d_est) - (b_est*exp(-d_est*t))/(2*d_est*(a_est - d_est)) + (b_est*exp(-a_est*t))/(2*a_est*(a_est - d_est));
plot(t,theta_est,'r');


%% Fit transfer function to position step response
% NOTE: This is what you can use to try and fit your transfer function
%       using experimental data.

% Estimate data for uniform sample frequency ------------------------------
% -> The fitting tool "tfest" used to estimate the coefficients of the 
%    transfer function requires a uniform sample frequency.

% Create a uniform time vector
uniform_t = linspace(min(t),max(t),numel(t))';
% Use linear interpolation to estimate theta for the uniform time vector 
uniform_theta = interp1(t,theta,uniform_t);
% -------------------------------------------------------------------------

% Add equilibrium values --------------------------------------------------
% -> The fitting tool "tfest" recommends at least some equilibrium values
%    when defining the fit. These are values for the system prior to 
%    applying the step response that we will prepend as zero.

% Define input signal as an array
PWM_in = repmat(PWM,size(theta));

% Prepend equilibrium values to input and output signal
% -> Define equilibrium time values
dt = mean( diff(uniform_t) );
prepend_t = [0:dt:(max(t)/5)]';
% -> Define equilibrium input values
prepend_PWM_in = zeros(size(prepend_t));
% -> Define equilibrium output values
prepend_theta = zeros(size(prepend_t));
% -> Prepend values
uniform_t = [prepend_t;...
             uniform_t - min(uniform_t) + dt + max(prepend_t)];
uniform_theta = [prepend_theta; uniform_theta];
PWM_in = [prepend_PWM_in; PWM_in];
% -------------------------------------------------------------------------

% Format collected position data
data_pos = iddata(...
    uniform_theta,...                   % Measured output  
    PWM_in,...                          % Applied input
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
% Position response
figure(3)
t_estimate = linspace(min(t),max(t),numel(t));
[theta_estimate,t_out] = step(PWM*H,t_estimate);
plot(axs_pos,t_out,theta_estimate,'g','LineWidth',2);

% Velocity response
figure(4)
[theta_dot_estimate,t_out]  = step(PWM*s*H,t_estimate);
plot(axs_vel,t_out,theta_dot_estimate,'g','LineWidth',2);

%% Display Transfer Function Coefficients
pvec = getpvec(H,'free');

a = pvec(2);
b = pvec(1);
d = pvec(3);

fprintf('Estimated Transfer function:\n');
fprintf('\t H(s) = %f / ( (s + %f)(s + %f) ) \n',b,a,d);
