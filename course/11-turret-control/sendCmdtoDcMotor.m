function [SSE,varargout] = sendCmdtoDcMotor(mode,control_params,varargin)
% sendCmdtoDcMotor emulates the behavior of the EW309 turret platforms
% under various user-specified operational conditions. The Primary output
% of the function is the steady-state error of the closed-loop system.
% However, optional input and output arguments can be specified to
% configure the system to simulate various open- or closed-loop
% configurations.
%
%   Inputs:
%       mode: control operational mode: 'step' or 'closed'
%               'step' mode applies a step input at t=1.0 seconds of
%               user-specified PWM duty cycle
%               'closed' mode uses a PID controller with user-specified
%               gains
%       control_params: data structure containing control input parameters
%        The following parameters are necessary when in closed-loop mode
%           control_params.despos: desired pointing angle in radians. This is the
%                                  reference input to the closed-loop controller
%           control_params.Kp: proportional gain on closed-loop controller
%           control_params.Ki: integral gain on closed-loop controller
%           control_params.Kd: derivative gain on closed loop controller
%        The following parameters are necessary when in step input mode
%           control_params.stepPWM: PWM duty cycle magnitude (0-1) of step
%           input (only applies to open-loop operational mode)
%
%   Optional inputs:
%       SSE = sendCmdtoDcMotor('closed',control_params,time)
%           time: time array to use in simulating the motor e.g. t=0:.1:3;
%       SSE = sendCmdtoDcMotor('closed',control_params,time,q0)
%           q0: initial condition of motor, as a vector of initial
%           position, initial speed, initial armature current, and initial
%           integral error, e.g.
%           theta0 = 0; dtheta=0; i0 = 0; e0 = 0; q0 = [theta0;dtheta0;i0;e0]
%
%
% Example usage (closed-loop control):
%   cntrlprms.despos = pi/4;
%   cntrlprms.Kp = 0.5;
%   cntrlprms.Ki = 0.02;
%   cntrlprms.Kd = 0.02;
%   t = 0:.05:10;
%   [SSE,t,theta,omega,eint] = sendCmdtoDcMotor('closed',cntrlprms,t);
% OR
%   SSE = sendCmdtoDcMotor('closed',cntrlprms,t);
%
%
% Example usage (open-loop step input):
%   cntrlprms.stepPWM = 0.45; % 45% duty cycle step input
%   t = 0:.05:10;
%   [SSE,t,theta,omega,eint] = sendCmdtoDcMotor('step',cntrlprms,t);
% OR
%   SSE = sendCmdtoDcMotor('step',cntrlprms,t,[0;0;0;0]);
%
%   L. DeVries, USNA, EW309, AY2020


% handle additional input quantities
if nargin == 5
    t = varargin{1};
    % initial condition
    theta0 = 0;
    dtheta0 = 0;
    i0 = 0;
    q0 = [theta0;dtheta0;i0;0];
elseif nargin == 6
    t = varargin{1};
    q0 = varargin{2};
elseif nargin > 6
    error('Too many inputs')
else
    t = 0:.01:10;
    % initial condition
    theta0 = 0;
    dtheta0 = 0;
    i0 = 0;
    q0 = [theta0;dtheta0;i0;0];
end


% Motor constants
motorParams.Ra = 5; % Armature resistance (Ohms)
motorParams.La = 0.2*10^-1; % Armature inductance (H) (~10^-3)
motorParams.Bm = .027; % coefficient of friction (Nm*s/rad)
motorParams.Km = .95; % transducer constant (Nm*s/rad) (amp*H/rad)
motorParams.J = 0.15*10^0; % moment of inertial
motorParams.friction.a0 = 0.15; % positive spin static friction (Nm)
motorParams.friction.a1 = 0.25; % positive spin coulumb friction coefficient
motorParams.friction.a2 = 1.3; % speed decay constant on coulumb friction
motorParams.friction.a3 = .36; % negative spin static friction (Nm)
motorParams.friction.a4 = 0.25; % negative spin coulumb friction coefficient
motorParams.friction.a5 = 1; % speed decay constant on coulumb friction
motorParams.friction.del = 0.05; % rad/s "linear zone" of friction
motorParams.dzone.pos = 0.25; % ten percent duty cycle on positive side 0.25 comes from trials
motorParams.dzone.neg = 0.25; % twenty percent on negative side 0.25 comes from trials


% switch operational modes (closed- or open-loop)
switch mode
    case 'closed'
        motorParams.case = 3; % closed loop control case
        % integrate EOM
        [~,Q] = ode45(@MotDynHF_sc,t,q0,[],motorParams,control_params);
        
        % steady-state error
        lng = ceil(0.05*length(Q(:,1))); % last 5% of data points
        SSE = control_params.despos - mean(Q(end-lng:end,1));
        
    case 'step'
        motorParams.case = 2; % step input case
        % integrate EOM
        [~,Q] = ode45(@MotDynHF_sc,t,q0,[],motorParams,control_params);
        
        % steady-state error (non-existent for step input)
        SSE = NaN;
        
        
end

% % % plot results
% % fig3 = figure(3); clf
% % plot(t,Q_cl(:,1))
% % hold on
% % plot(t,cntrlprms.despos*ones(size(t)),'--r')
% % legend('Closed-loop Response', 'Desired Position')
% % xlabel('Time (s)')
% % ylabel('Orientation (rad)')

% orientation
out(:,1) = t; % time vector
out(:,2) = Q(:,1); % theta
out(:,3) = Q(:,2); % omega
out(:,4) = Q(:,4); % error integral


% Optional outputs
nout = max(nargout,1)-1;
for k = 1:nout
    varargout{k} = out(:,k);
end

end


function dQ = MotDynHF_sc(t,Q,params,cntrlprms)
%MotDynHF simulates high fidelity dynamics of a DC motor. The model
%includes nonlinear input deadzone and Stribeck Friction.
%   INPUTS:
%       t: time, scalar value for current time of integration
%       Q: 4x1 dimensional state vector at time t, Q = [position; velocity;
%       current; error integral]
%       params: data structure containing motor parameters
%           params.Ra: motor armature resistance (Ohms)
%           params.La: motor armature inducatnce (H)
%           params.Bm: coefficient of linear friction (Nm*s/rad)
%           params.Km: transducer constant (Nm*s/rad) (amp*H/rad)
%           params.J: moment of inertia (Kg*m^2)
%           params.friction.a0: positive spin static friction (Nm)
%           params.friction.a1: positive spin coulumb friction coefficient (Nm)
%           params.friction.a2: speed decay constant on coulumb friction (unitless
%           params.friction.a3: negative spin static friction (Nm)
%           params.friction.a4: negative spin coulumb friction coefficient
%           params.friction.a5: speed decay constant on coulumb friction
%           params.friction.del: approximation of stiction range (rad/s)
%           params.dzone.pos: dead zone for positive inputs (duty cycle)
%           params.dzone.neg: dead zone for negative inputs (duty cycle)
%           params.case: operational case to simulate
%                       1 = sinusoidal input, amplitude 1, period 10 second
%                       2 = step input, user-specified magnitude (0-1),
%                       step input is applied at t=1.0 second
%                       3 = closed loop control
%            
%       cntrlprms: data structure containing operational mode (open- or closed-loop)
%                  and parameters of operation (PID control gains or
%                  open-loop step input duty cycle)
%           cntrlprms.mode: operational mode ('open' or 'closed')
%           cntrlprms.despos: desired position (rad)
%           cntrlprms.Kp: proportional gain
%           cntrlprms.Ki: integral gain
%           cntrlprms.Kd: derivative gain
%           cntrlprms.stepPWM: duty cycle magnitude of step input (0-1)
%   OUTPUT:
%       dQ: 4x1 derivative of state vector Q at time t
%
%   Example Usage: (open-loop sine wave input)
%       motorParams.Ra = 5; % Armature resistance (Ohms)
%       motorParams.La = 0.2*10^-1; % Armature inductance (H) (~10^-3)
%       motorParams.Bm = .027; % coefficient of friction (Nm*s/rad)
%       motorParams.Km = .95; % transducer constant (Nm*s/rad) (amp*H/rad)
%       motorParams.J = 0.15*10^0; % moment of inertial
%       motorParams.friction.a0 = 0.15; % positive spin static friction (Nm)
%       motorParams.friction.a1 = 0.25; % positive spin coulumb friction coefficient
%       motorParams.friction.a2 = 1.3; % speed decay constant on coulumb friction 
%       motorParams.friction.a3 = .36; % negative spin static friction (Nm)
%       motorParams.friction.a4 = 0.25; % negative spin coulumb friction coefficient
%       motorParams.friction.a5 = 1; % speed decay constant on coulumb friction
%       motorParams.friction.del = 0.05; % rad/s "linear zone" of friction
%       motorParams.dzone.pos = 0.25; % ten percent duty cycle on positive side 0.25 comes from trials 
%       motorParams.dzone.neg = 0.25; % twenty percent on negative side 0.25 comes from trials
%       cntrlprms.despos = 0;
%       cntrlprms.Kp = 0;
%       cntrlprms.Ki = 0;
%       cntrlprms.Kd = 0;
%       cntrlprms.stepPWM = 0.45;
% 
%       % initial conditions
%       t = 0:.01:3;
%       theta0 = 0; % position
%       dtheta0 = 0; % angular velocity
%       i0 = 0; % initial current
%       q0 = [theta0;dtheta0;i0;0]; % initial state vector
%       motorParams.case = 1; % (for testing/development) case one, sinusoidal input
%       [~,Q] = ode45(@MotDynHF,t,q0,[],motorParams,cntrlprms);
%
% L. DeVries, Ph.D., USNA
% EW309, AY2020
% Last edited 3/22/2020



% control input
switch params.case % test cases for simulation 
    case 1 % sinusoidal data, compare to experimental data
        if t<0.45
            dc = 0;
        else
            dc = sin(2*pi/10*t);
        end
        err = 0;
    case 2 % positive step input, compare to experimental data
        if t<1
            dc = 0;
        else
            dc = cntrlprms.stepPWM;
        end
        err = 0;
    case 3 % closed loop control, for students
        err = cntrlprms.despos - Q(1);
 
        % PID controller
        dc = cntrlprms.Kp*err + cntrlprms.Ki*Q(4) - cntrlprms.Kd*Q(2);
        if dc>1.0
            dc = 1.0;
        elseif dc<-1.0
            dc = -1.0;
        end
end

% introduce deadzone here! Use EW305 results to quantify!
if dc<params.dzone.pos && dc>-params.dzone.neg
    dc = 0;
end

% voltage from duty cycle
u = 12*dc;

% motor torque when assuming inductance is negligible
ia = 1/params.Ra*(u-params.Km*Q(2)); % armature current
Tm = params.Km*ia; % torque

% use this motor torque when assuming inductance is not negligible
% Tm = params.Km*Q(3);

% equations of motion (EOM)
dQ(1,1) = Q(2); % dtheta = omega
dQ(2,1) = 1/params.J*(Tm - params.Bm*Q(2) - frictionNL_sc(Q(2),params)); % domega = torques and friction

% electrical circuit EOM (assumes non-negligible inductance) (set to zero
% if ignoring inductance)
dQ(3,1) = 0;%1/params.La*(u - params.Km*Q(2) - params.Ra*Q(3)); % di/dt = circuit equation with inductance

% avoid windup
alp = (err/50).^4; % weighting function approaches zero sharply as error approaches 100 rad*s

dQ(4,1) = (1-alp)*err; % error integral term
end



function Tf = frictionNL_sc(omega,params)

Tf = (params.friction.a0 + params.friction.a1)/params.friction.del*omega;

indmd = omega <0 & omega > -params.friction.del;
Tf(indmd) = (params.friction.a3 + params.friction.a4)/params.friction.del*omega(indmd);

indp = omega>=params.friction.del;
Tf(indp) =   params.friction.a0 + params.friction.a1*exp( - (params.friction.a2*abs( omega(indp) - params.friction.del) )) ;

indm = omega <= -params.friction.del;
Tf(indm) = -(params.friction.a3 + params.friction.a4*exp( - (params.friction.a5*abs( omega(indm) + params.friction.del) )));

end
