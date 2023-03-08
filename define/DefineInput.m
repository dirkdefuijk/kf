function [u] = DefineInput(flightdata)
% DefineInput.m is a script which determines the input for the Kalman
% filter automatically. The script automatically creates a 6xM matrix,
% where M is the number of data points. The 6 inputs are ordered as: 
% ax, ay, az, p, q, r.
%
% Inputs:   
%           flightdata: structure containing the measured flightdata
%
% Outputs:  
%           u:          input matrix [6xM matrix]
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Define constants
% Define the gravitational acceleration on Earth
g0 = 9.80665;

%% Obtain input matrix with correts SI units (m/s^2 & rad/s)
    
% Obtain the acceleration in x-direction
switch flightdata.arinc.Ahrs1.bLongAcc.Unit
    case 'g'
        ax = flightdata.arinc.Ahrs1.bLongAcc.data*g0;
    case 'm/s^2'
        ax = flightdata.arinc.Ahrs1.bLongAcc.data;
    otherwise
        error('Unit of the linear acceleration in x-direction is unknown')
end

% Obtain the acceleration in y-direction    
switch flightdata.arinc.Ahrs1.bLatAcc.Unit
    case 'g'
        ay = flightdata.arinc.Ahrs1.bLatAcc.data*g0;
    case 'm/s^2'
        ay = flightdata.arinc.Ahrs1.bLatAcc.data;
    otherwise
        error('Unit of the linear acceleration in y-direction is unknown')
end

% Obtain the acceleration in z-direction
switch flightdata.arinc.Ahrs1.bNormAcc.Unit
    case 'g'
        az = -flightdata.arinc.Ahrs1.bNormAcc.data*g0;
    case 'm/s^2'
        az = -flightdata.arinc.Ahrs1.bNormAcc.data;
    otherwise
        error('Unit of the linear acceleration in z-direction is unknown')
end

% Obtain the roll rate
switch flightdata.arinc.Ahrs1.bRollRate.Unit
    case 'deg/s'
        p = flightdata.arinc.Ahrs1.bRollRate.data*pi/180;
    case 'rad/s'
        p = flightdata.arinc.Ahrs1.bRollRate.data;
    otherwise
        error('Unit of the roll rate is unknown')
end

% Obtain the pitch rate    
switch flightdata.arinc.Ahrs1.bPitchRate.Unit
    case 'deg/s'
        q = flightdata.arinc.Ahrs1.bPitchRate.data*pi/180;
    case 'rad/s'
        q = flightdata.arinc.Ahrs1.bPitchRate.data;
    otherwise
        error('Unit of the pitch rate is unknown')
end

% Obtain the yaw rate
switch flightdata.arinc.Ahrs1.bYawRate.Unit
    case 'deg/s'
        r = flightdata.arinc.Ahrs1.bYawRate.data*pi/180;
    case 'rad/s'
        r = flightdata.arinc.Ahrs1.bYawRate.data;
    otherwise
        error('Unit of the yaw rate is unknown')
end
    
%% Define input matrix
u = [ax; ay; az; p; q; r];
end