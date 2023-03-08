function [datastruct] = InitDataStructure(type_kf)
% InitDataStructure.m is a script which is used to initialize the
% datastructure for the flight path reconstructed (FPR) variables. It general
% information on, e.g. the type of filter that is used. Furthermore it
% stores the inputs, state and measurements of the aircraft
%
% Inputs:   
%           type_kf:            The type of kalman filter that is used
%
% Outputs:  
%           datastruct:         Initializing data structure for FPR values
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Information 
% Create version information. For thesis of Marlon van den Hoek and
% Laurens van Horssen, this will be version 1.0
info.Version = '1.0';

% Obtain information regarding the gear and flaps settings. If these values
% are unknown (e.g. for the flight practical data from 2015 and earlier),
% then set the value to "Unknown".
info.Configuration.gear = ''; % [1 = Down, 0 = Up]
info.Configuration.flaps = ''; % [degrees]

% Obtain the modification date
info.Modification_date = date;

% Choose which type of filter is used. The following options are available:
% 'EKF'     : Extended Kalman Filter 
% 'IEKF'    : Iterated Extended Kalman Filter
% 'UKF'     : Unscented Kalman Filter
% 'ERTSS'  :  Extended Rauch Tung Striebel Smoother
% 'IERTSS'  : Iterated Extended Rauch Tung Striebel Smoother
% 'URTSS'   : Unscented Rauch Tung Striebel Smoother
info.FilterType = '';

% Initial condition vector of size Nx1, where N is the dimension of the
% kinematic system
info.FilterParameters.x0.data = [];
info.FilterParameters.x0.Description = 'Initial condition vector';

% Initial covariance matrix of size NxN, where N is the dimension of the
% kinematic system
info.FilterParameters.P0.data = [];
info.FilterParameters.P0.Description = 'Initial covariance matrix';

% Process/ Input noise matrix of size MxM, where M is the number of inputs
info.FilterParameters.Q.data = [];
info.FilterParameters.Q.Description = 'Process noise matrix';

% Measurement noise matrix of size LxL, where L is the number of
% measurements
info.FilterParameters.R.data = []; 
info.FilterParameters.R.Description = 'Measurement noise matrix';

if strcmp(type_kf,'UKF') == 1 || strcmp(type_kf,'URTSS') == 1 
    % Scaling parameter for the UKF. alpha determines the spread of the sigma
    % points around the mean state
    info.FilterParameters.alpha.data = [];
    info.FilterParameters.alpha.Description = 'Scaling parameter UKF';

    % Secondy scaling parameter, usually set to zero.
    info.FilterParameters.kappa.data = [];
    info.FilterParameters.kappa.Description = 'Scaling parameter UKF';

    % Scaling parameter to incorporate prior knowledge of the distribution of
    % the state. For Gaussian distribution beta = 2 is optimal.
    info.FilterParameters.beta.data = [];
    info.FilterParameters.beta.Description = 'Scaling parameter UKF';
    
elseif  strcmp(type_kf,'IEKF') == 1 || strcmp(type_kf,'IERTSS') == 1
    % Epsilon determines that maximum error between two estimated states,
    % before the iteration is ended
    info.FilterParameters.epsilon.data = [];
    info.FilterParameters.epsilon.Description = 'Maximum error before stopping iterations of the iteration';

    % maxIterations determines the maximum number of iterations that were
    % allowed in case of an iterated Kalman filter
    info.FilterParameters.maxIterations.data = [];
    info.FilterParameters.maxIterations.Description = 'Maximum number of iterations that were allowed';
elseif strcmp(type_kf,'EKF') == 1 || strcmp(type_kf,'ERTSS') == 1 
    % Do nothing, since there are no options to choose    
else
    error('Type of Kalman filter that is used is wrong. Options that could be chosen are: EKF, IEKF, UKF, ERTSS, IERTSS, URTSS')
end

% Inform the user whether or not the inputs have been prefiltered. Options
% are either: 'yes' or 'no'
info.Prefiltered.true = ''; % yes or no

% In case the inputs have been prefiltered, give the cut-off frequency of
% the filter that has been used. If no prefilter is used, leave the
% parameter empty.
info.Prefiltered.fc.data = [];
info.Prefiltered.fc.Description = 'Cut off frequency';
info.Prefiltered.fc.Unit = ''; % hz or rad

% Give information on the filter order that is used for prefiltering. If
% no prefilter is used, leave the parameter empty.
info.Prefiltered.fo.data = [];
info.Prefiltered.fo.Description = 'Filter order';
info.Prefiltered.fo.Unit = '-';

% Give information on the type of filter that is used for prefiltering the
% inputs. If no prefilter is used, leave the parameter empty.
info.Prefiltered.type.data = ''; % e.g. Butterworth, Chebyshev 
info.Prefiltered.type.Description = 'Type of filter that has been used to filter inputs';

clearvars -except info

%% Reconstruct the input
% Acceleration in x-direction
u.ax.time = [];
u.ax.raw = [];
u.ax.filtered = [];
u.ax.Description = 'Acceleration in x-direction';
u.ax.Unit = 'm/s^2';
u.ax.Sensor = ''; % 'arinc/Ahrs1/bLongAcc'

% Acceleration in y-direction
u.ay.time = [];
u.ay.raw = [];
u.ay.filtered = [];
u.ay.Description = 'Acceleration in y-direction';
u.ay.Unit = 'm/s^2';
u.ay.Sensor = ''; % 'arinc/Ahrs1/bLatAcc'

% Acceleration in z-direction
u.az.time = [];
u.az.raw = [];
u.az.filtered = [];
u.az.Description = 'Acceleration in z-direction';
u.az.Unit = 'm/s^2';
u.az.Sensor = ''; % 'arinc/Ahrs1/bNormAcc'

% Roll rate
u.p.time = [];
u.p.raw = [];
u.p.filtered = [];
u.p.Description = 'Roll rate';
u.p.Unit = 'rad/s';
u.p.Sensor = ''; % 'arinc/Ahrs1/bRollRate'

% Pitch rate
u.q.time = [];
u.q.raw = [];
u.q.filtered = [];
u.q.Description = 'Pitch rate';
u.q.Unit = 'rad/s';
u.q.Sensor = ''; % 'arinc/Ahrs1/bPitchRate'

% Yaw rate
u.r.time = [];
u.r.raw = [];
u.r.filtered = [];
u.r.Description = 'Yaw rate';
u.r.Unit = 'rad/s';
u.r.Sensor = ''; % 'arinc/Ahrs1/bYawRate'

% Elevator deflection
delta.de.time = [];
delta.de.raw = [];
delta.de.filtered = [];
delta.de.Description = 'Elevator deflection';
delta.de.Unit = 'rad';
delta.de.Sensor = ''; % 'synchro/delta_e'

% Aileron deflection
delta.da.time = [];
delta.da.raw = [];
delta.da.filtered = [];
delta.da.Description = 'Aileron deflection';
delta.da.Unit = 'rad';
delta.da.Sensor = ''; % 'synchro/delta_a'

% Rudder deflection
delta.dr.time = [];
delta.dr.raw = [];
delta.dr.filtered = [];
delta.dr.Description = 'Rudder deflection';
delta.dr.Unit = 'rad';
delta.dr.Sensor = ''; % 'synchro/delta_r'

% Trim tab deflection
delta.dte.time = [];
delta.dte.raw = [];
delta.dte.filtered = [];
delta.dte.Description = 'Elevator trim tab';
delta.dte.Unit = 'rad';
delta.dte.Sensor = ''; % ''

% Spoiler setting
delta.spoiler.time = [];
delta.spoiler.raw = [];
delta.spoiler.filtered = [];
delta.spoiler.Description = 'Spoiler deflection';
delta.spoiler.Unit = 'rad';
delta.spoiler.Sensor = ''; % ''

% Flap setting
delta.flaps.time = [];
delta.flaps.raw = [];
delta.flaps.filtered = [];
delta.flaps.Description = 'Flap deflection';
delta.flaps.Unit = 'deg';
delta.flaps.Sensor = ''; % 'analog/flaps'

% Gear in or out
delta.gear.time = [];
delta.gear.data = [];
delta.gear.Description = 'Landing gear in/out';
delta.gear.Unit = 'In/Out - (0/1)';
delta.gear.Sensor = ''; %

% Roll acceleration
RotAcc.pdot.time = [];
RotAcc.pdot.data = [];
RotAcc.pdot.Description = 'Roll acceleration';
RotAcc.pdot.Unit = 'rad/s^2';

% Pitch acceleration
RotAcc.qdot.time = [];
RotAcc.qdot.data = [];
RotAcc.qdot.Description = 'Pitch acceleration';
RotAcc.qdot.Unit = 'rad/s^2';

% Yaw acceleration
RotAcc.rdot.time = [];
RotAcc.rdot.data = [];
RotAcc.rdot.Description = 'Yaw acceleration';
RotAcc.rdot.Unit = 'rad/s^2';

% Fuel mass flow left engine
Engine.lh.FMF.time = [];
Engine.lh.FMF.data = [];
Engine.lh.FMF.Description = 'Left engine: Fuel mass flow';
Engine.lh.FMF.Unit = 'kg/s';
Engine.lh.FMF.Sensor = ''; % 'analog/lh_engine_FMF'

% Oil pressure left engine
Engine.lh.OP.time = [];
Engine.lh.OP.data = [];
Engine.lh.OP.Description = 'Left engine: Oil pressure';
Engine.lh.OP.Unit = 'psi';
Engine.lh.OP.Sensor = ''; % 'analog/lh_engine_OP'

% Inter Turbine Temperature left engine
Engine.lh.ITT.time = [];
Engine.lh.ITT.data = [];
Engine.lh.ITT.Description = 'Left engine: Inter Turbine Temperature';
Engine.lh.ITT.Unit = 'C';
Engine.lh.ITT.Sensor = ''; % 'analog/lh_engine_itt'

% Fan speed left engine
Engine.lh.N1.time = [];
Engine.lh.N1.data = [];
Engine.lh.N1.Description = 'Left engine: Fan speed (N1)';
Engine.lh.N1.Unit = '%';
Engine.lh.N1.Sensor = ''; % 'analog/lh_engine_fan_N1'

% Turbine speed left engine
Engine.lh.N2.time = [];
Engine.lh.N2.data = [];
Engine.lh.N2.Description = 'Left engine: Turbine speed (N2)';
Engine.lh.N2.Unit = '%';
Engine.lh.N2.Sensor = ''; % 'analog/lh_engine_turbine_N2'

% Fuel used left engine
Engine.lh.FU.time = [];
Engine.lh.FU.data = [];
Engine.lh.FU.Description = 'Left engine: fuel used';
Engine.lh.FU.Unit = 'kg';

% Fuel mass flow right engine
Engine.rh.FMF.time = [];
Engine.rh.FMF.data = [];
Engine.rh.FMF.Description = 'Right engine: Fuel mass flow';
Engine.rh.FMF.Unit = 'kg/s';
Engine.rh.FMF.Sensor = ''; % 'analog/rh_engine_FMF'

% Oil pressure right engine
Engine.rh.OP.time = [];
Engine.rh.OP.data = [];
Engine.rh.OP.Description = 'Right engine: Oil pressure';
Engine.rh.OP.Unit = 'psi';
Engine.rh.OP.Sensor = ''; % 'analog/rh_engine_OP'

% Inter Turbine Temperature right engine
Engine.rh.ITT.time = [];
Engine.rh.ITT.data = [];
Engine.rh.ITT.Description = 'Right engine: Inter Turbine Temperature';
Engine.rh.ITT.Unit = 'C';
Engine.rh.ITT.Sensor = ''; % 'analog/rh_engine_itt'

% Fan speed right engine
Engine.rh.N1.time = [];
Engine.rh.N1.data = [];
Engine.rh.N1.Description = 'Right engine: Fan speed (N1)';
Engine.rh.N1.Unit = '%';
Engine.rh.N1.Sensor = ''; % 'analog/rh_engine_fan_N1'

% Turbine speed right engine
Engine.rh.N2.time = [];
Engine.rh.N2.data = [];
Engine.rh.N2.Description = 'Right engine: Turbine speed (N2)';
Engine.rh.N2.Unit = '%';
Engine.rh.N2.Sensor = ''; % 'analog/rh_engine_turbine_N2'

% Fuel used right engine
Engine.rh.FU.time = [];
Engine.rh.FU.data = [];
Engine.rh.FU.Description = 'Right engine: fuel used';
Engine.rh.FU.Unit = 'kg';

% Total engine thrust
Engine.Thrust.time = [];
Engine.Thrust.data = [];
Engine.Thrust.Description = 'Total engine thrust';
Engine.Thrust.Unit = 'N';

Input = struct('u', u, 'delta', delta, 'RotAcc', RotAcc);
clearvars -except info Input Engine

%% Reconstruct the state
% x-position in the Earth Centered Earth Fixed frame
xpos.time = [];
xpos.data = [];
xpos.std = [];
xpos.Description = 'x-position in ECEF frame';
xpos.Unit = 'm';

% y-position in the Earth Centered Earth Fixed frame
ypos.time = [];
ypos.data = []; 
ypos.std = [];
ypos.Description = 'y-position in ECEF frame';
ypos.Unit = 'm';

% z-position in the Earth Centered Earth Fixed frame
zpos.time = [];
zpos.data = []; 
zpos.std = [];
zpos.Description = 'z-position in ECEF frame';
zpos.Unit = 'm';

% Body velocity in x-direction
u.time = [];
u.data = [];
u.std = [];
u.Description = 'Longitudinal body velocity';
u.Unit = 'm/s';

% Body velocity in y-direction
v.time = [];
v.data = []; 
v.std = [];
v.Description = 'Lateral body velocity';
v.Unit = 'm/s';

% Body velocity in z-direction
w.time = [];
w.data = []; 
w.std = [];
w.Description = 'Vertical body velocity';
w.Unit = 'm/s';

% Roll angle
phi.time = [];
phi.data = [];
phi.std = [];
phi.Description = 'Roll angle';
phi.Unit = 'rad';

% Pitch angle
theta.time = [];
theta.data = []; 
theta.std = [];
theta.Description = 'Pitch angle';
theta.Unit = 'rad';

% Yaw angle
psi.time = [];
psi.data = []; 
psi.std = [];
psi.Description = 'Yaw angle';
psi.Unit = 'rad';

% Wind velocity in x-direction in the Earth Centered Earth Fixed frame
Wxe.time = [];
Wxe.data = [];
Wxe.std = [];
Wxe.Description = 'Wind velocity in x-direction (ECEF frame)';
Wxe.Unit = 'm/s';

% Wind velocity in y-direction in the Earth Centered Earth Fixed frame
Wye.time = [];
Wye.data = [];
Wye.std = [];
Wye.Description = 'Wind velocity in y-direction (ECEF frame)';
Wye.Unit = 'm/s';

% Wind velocity in z-direction in the Earth Centered Earth Fixed frame
Wze.time = [];
Wze.data = [];
Wze.std = [];
Wze.Description = 'Wind velocity in z-direction (ECEF frame)';
Wze.Unit = 'm/s';

% Fuselage upwash component effect on the angle of attack
Caup.time = [];
Caup.data = [];
Caup.std = [];
Caup.Description = 'Fuselage upwash component effect on angle of attack';
Caup.Unit = '-';

% Unknown wind component in the angle of attack
Ca0.time = [];
Ca0.data = [];
Ca0.std = [];
Ca0.Description = 'Unknown wind component in angle of attack';
Ca0.Unit = 'rad';

% Fuselage sidewash component on the angle of side slip
Cbsi.time = [];
Cbsi.data = [];
Cbsi.std = [];
Cbsi.Description = 'Fuselage sidewash component effect on angle of side slip';
Cbsi.Unit = '-';

% Unknown wind component in the angle of side slip
Cb0.time = [];
Cb0.data = [];
Cb0.std = [];
Cb0.Description = 'Unknown wind component in angle of side slip';
Cb0.Unit = 'rad';

% Angle of attack as measured by the alpha-vane
alpha_v.time = [];
alpha_v.data = [];
alpha_v.std = [];
alpha_v.Description = 'Angle of attack as measured by alpha vane';
alpha_v.Unit = 'rad';

% Bias on acceleration in x-direction
lambda_Ax.time = [];
lambda_Ax.data = []; 
lambda_Ax.std = [];
lambda_Ax.Description = 'Bias on acceleration in x-direction';
lambda_Ax.Unit = 'm/s^2';

% Bias on acceleration in y-direction
lambda_Ay.time = [];
lambda_Ay.data = []; 
lambda_Ay.std = [];
lambda_Ay.Description = 'Bias on acceleration in y-direction';
lambda_Ay.Unit = 'm/s^2';

% Bias on acceleration in z-direction
lambda_Az.time = [];
lambda_Az.data = []; 
lambda_Az.std = [];
lambda_Az.Description = 'Bias on acceleration in z-direction';
lambda_Az.Unit = 'm/s^2';

% Bias on roll rate measurement
lambda_p.time = [];
lambda_p.data = []; 
lambda_p.std = [];
lambda_p.Description = 'Bias on roll rate';
lambda_p.Unit = 'rad/s';

% Bias on pitch rate measurement
lambda_q.time = [];
lambda_q.data = []; 
lambda_q.std = [];
lambda_q.Description = 'Bias on pitch rate';
lambda_q.Unit = 'rad/s';

% Bias on yaw rate measurement
lambda_r.time = [];
lambda_r.data = []; 
lambda_r.std = [];
lambda_r.Description = 'Bias on yaw rate';
lambda_r.Unit = 'rad/s';

% Body velocity in x-direction - left wing
u_L.time = [];
u_L.data = [];
u_L.std = [];
u_L.Description = 'Longitudinal body velocity - left wing';
u_L.Unit = 'm/s';

% Body velocity in x-direction - right wing
u_R.time = [];
u_R.data = [];
u_R.std = [];
u_R.Description = 'Longitudinal body velocity - right wing';
u_R.Unit = 'm/s';

% Body velocity in y-direction - left wing
v_L.time = [];
v_L.data = []; 
v_L.std = [];
v_L.Description = 'Lateral body velocity - left wing';
v_L.Unit = 'm/s';

% Body velocity in y-direction - right wing
v_R.time = [];
v_R.data = []; 
v_R.std = [];
v_R.Description = 'Lateral body velocity - right wing';
v_R.Unit = 'm/s';

% Body velocity in z-direction - left wing
w_L.time = [];
w_L.data = []; 
w_L.std = [];
w_L.Description = 'Vertical body velocity - left wing';
w_L.Unit = 'm/s';

% Body velocity in z-direction - right wing
w_R.time = [];
w_R.data = []; 
w_R.std = [];
w_R.Description = 'Vertical body velocity - right wing';
w_R.Unit = 'm/s';


State = struct('x', xpos, 'y', ypos, 'z', zpos, ...
               'u', u, 'v', v, 'w', w, ...
               'phi', phi, 'theta', theta, 'psi', psi, ...
               'Wxe', Wxe, 'Wye', Wye, 'Wze', Wze, ...
               'Caup', Caup, 'Ca0', Ca0, 'Cbsi', Cbsi, 'Cb0', Cb0, 'alpha_v', alpha_v, ...
               'lambda_Ax', lambda_Ax, 'lambda_Ay', lambda_Ay, 'lambda_Az', lambda_Az, ...
               'lambda_p', lambda_p,'lambda_q', lambda_q, 'lambda_r', lambda_r, ...
               'u_L', u_L, 'u_R', u_R, 'v_L', v_L, 'v_R', v_R, 'w_L', w_L, 'w_R', w_R);

clearvars -except info Input Engine State

%% Reconstruct the measurements\
% x-position in the Earth Centered Earth Fixed frame
xpos.time = [];
xpos.raw = []; 
xpos.reconstructed = []; 
xpos.innovation = []; 
xpos.std = [];
xpos.Description = 'x-position in ECEF frame';
xpos.Unit = 'm';
xpos.Sensor = ''; % 'arinc/Gps/xposition'

% y-position in the Earth Centered Earth Fixed frame
ypos.time = [];
ypos.raw = []; 
ypos.reconstructed = []; 
ypos.innovation = []; 
ypos.std = [];
ypos.Description = 'y-position in ECEF frame';
ypos.Unit = 'm';
ypos.Sensor = ''; % 'arinc/Gps/yposition'

% z-position in the Earth Centered Earth Fixed frame
alt.time = [];
alt.raw = []; 
alt.reconstructed = []; 
alt.innovation = []; 
alt.std = [];
alt.Description = 'Altitude';
alt.Unit = 'm';
alt.Sensor = ''; % 'arinc/Dadc1/alt'

% Velocity in x-direction in Earth Centered Earth Fixed frame
xdot.time = [];
xdot.raw = []; 
xdot.reconstructed = []; 
xdot.innovation = []; 
xdot.std = [];
xdot.Description = 'velocity in x-direction in ECEF frame';
xdot.Unit = 'm/s';
xdot.Sensor = ''; % 'arinc/Gps/nsVel'

% Velocity in y-direction in Earth Centered Earth Fixed frame
ydot.time = [];
ydot.raw = []; 
ydot.reconstructed = []; 
ydot.innovation = []; 
ydot.std = [];
ydot.Description = 'velocity in x-direction in ECEF frame';
ydot.Unit = 'm/s';
ydot.Sensor = ''; % 'arinc/Gps/ewVel'

% Altitude rate
hdot.time = [];
hdot.raw = []; 
hdot.reconstructed = []; 
hdot.innovation = []; 
hdot.std = [];
hdot.Description = 'Altitude rate';
hdot.Unit = 'm/s';
hdot.Sensor = ''; % 'arinc/Dadc1/altRate'

% Roll angle
phi.time = [];
phi.raw = []; 
phi.reconstructed = []; 
phi.innovation = []; 
phi.std = [];
phi.Description = 'Roll angle';
phi.Unit = 'rad';
phi.Sensor = ''; % 'arinc/Ahrs1/Roll'

% Pitch angle
theta.time = [];
theta.raw = []; 
theta.reconstructed = []; 
theta.innovation = []; 
theta.std = [];
theta.Description = 'Pitch angle';
theta.Unit = 'rad';
theta.Sensor = ''; % 'arinc/Ahrs1/Pitch'

% Yaw angle
psi.time = [];
psi.raw = []; 
psi.reconstructed = []; 
psi.innovation = []; 
psi.std = [];
psi.Description = 'Yaw angle';
psi.Unit = 'rad';
psi.Sensor = ''; % 'arinc/Fms1/trueHdg'

% True airspeed
Vtas.time = [];
Vtas.raw = []; 
Vtas.reconstructed = []; 
Vtas.innovation = []; 
Vtas.std = [];
Vtas.Description = 'True airspeed';
Vtas.Unit = 'm/s';
Vtas.Sensor = ''; % 'arinc/Dadc1/tas'

% Angle of attack
alpha.time = [];
alpha.raw = []; 
alpha.reconstructed = []; 
alpha.innovation = []; 
alpha.std = [];
alpha.Description = 'Angle of attack';
alpha.Unit = 'rad';
alpha.Sensor = ''; % 'analog/vane_AOA'

% Angle of side slip
beta.time = [];
beta.raw = []; 
beta.reconstructed = []; 
beta.innovation = []; 
beta.std = [];
beta.Description = 'Angle of side slip';
beta.Unit = 'rad';
beta.Sensor = ''; % 'Pseudo Beta'

% Angle of attack - left wing
alpha_L.time = [];
alpha_L.raw = []; 
alpha_L.reconstructed = []; 
alpha_L.innovation = []; 
alpha_L.std = [];
alpha_L.Description = 'Angle of attack -  left wing';
alpha_L.Unit = 'rad';
alpha_L.Sensor = ''; % 'analog/vane_AOA'

% Angle of attack - right wing
alpha_R.time = [];
alpha_R.raw = []; 
alpha_R.reconstructed = []; 
alpha_R.innovation = []; 
alpha_R.std = [];
alpha_R.Description = 'Angle of attack -  right wing';
alpha_R.Unit = 'rad';
alpha_R.Sensor = ''; % 'analog/vane_AOA'

Measurement = struct('x', xpos, 'y', ypos, 'h', alt, ...
                     'xdot', xdot, 'ydot', ydot, 'hdot', hdot, ...
                     'phi', phi, 'theta', theta, 'psi', psi, ...
                     'Vtas', Vtas, 'alpha', alpha, 'beta', beta, ...
                     'alpha_L', alpha_L, 'alpha_R', alpha_R);
                 
clearvars -except info Input Engine State Measurement

%% Put all structures together to create one final data structure
datastruct = struct('info', info, 'Input', Input, 'State', State, ...
    'Measurement', Measurement, 'Engine', Engine);
  
end