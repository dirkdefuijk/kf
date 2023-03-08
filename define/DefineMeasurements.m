function [y,DataStruct] = DefineMeasurements(flightdata,DataStruct)
% DefineMeasurements.m is a script which is used to initialize the
% measurement. The measurements for the new data are:
%
% x, y, h, xdot, ydot, hdot, phi, theta, psi, Vtas, alpha (vane), pseudo beta
%
% Inputs:   
%           flightdata:         structure containing the measured flightdata
%           DataStruct:         Stores sensor information in the DataStruct
%
% Outputs:  
%           y:                  measurement matrix [matrix]
%           DataStruct:         DataStructure storing information with, e.g. the state
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Store the measurements in the observation matrix

% Obtain x-position
switch flightdata.arinc.Gps.xposition.Unit
    case 'ft'
        xpos = flightdata.arinc.Gps.xposition.data*0.3048;
    case 'm'
        xpos = flightdata.arinc.Gps.xposition.data;
    otherwise
        error('Unit of the x-position of the GPS is unknown')
end
DataStruct.Measurement.x.Sensor = 'arinc/Gps/xposition';

% Obtain y-position
switch flightdata.arinc.Gps.yposition.Unit
    case 'ft'
        ypos = flightdata.arinc.Gps.yposition.data*0.3048;
    case 'm'
        ypos = flightdata.arinc.Gps.yposition.data;
    otherwise
        error('Unit of the y-position of the GPS is unknown')
end
DataStruct.Measurement.y.Sensor = 'arinc/Gps/yposition';

% Obtain altitude
switch flightdata.arinc.Gps.alt.Unit
    case 'ft'
        alt = flightdata.arinc.Gps.alt.data*0.3048; 
    case 'm'
        alt = flightdata.arinc.Gps.alt.data;
    otherwise
        error('Unit of the altitude is unknown')
end
DataStruct.Measurement.h.Sensor = 'arinc/Gps/alt';

% Obtain xdot (North-South Velocity)
switch flightdata.arinc.Gps.nsVel.Unit
    case 'knots'
        xdot  = flightdata.arinc.Gps.nsVel.data*0.514444; 
    case 'm/s'
        xdot  = flightdata.arinc.Gps.nsVel.data;
    otherwise
        error('Unit of the North-South velocity of the GPS is unknown')
end
DataStruct.Measurement.xdot.Sensor = 'arinc/Gps/nsVel';
    
% Obtain ydot (East-West Velocity)
switch flightdata.arinc.Gps.ewVel.Unit
    case 'knots'
        ydot  = flightdata.arinc.Gps.ewVel.data*0.514444; 
    case 'm/s'
        ydot  = flightdata.arinc.Gps.ewVel.data;
    otherwise
        error('Unit of the East-West velocity of the GPS is unknown')
end
DataStruct.Measurement.ydot.Sensor = 'arinc/Gps/ewVel';

% Obtain hdot (altitude rate)
switch flightdata.arinc.Gps.vertVel.Unit
    case 'ft/min'
        hdot  = flightdata.arinc.Gps.vertVel.data*0.00508; 
    case 'm/s'
        hdot  = flightdata.arinc.Gps.vertVel.data;
    otherwise
        error('Unit of the altitude rate is unknown')
end
DataStruct.Measurement.hdot.Sensor = 'arinc/Gps/vertVel';

% Obtain phi (roll angle)
switch flightdata.arinc.Ahrs1.Roll.Unit
    case 'deg'
        phi  = flightdata.arinc.Ahrs1.Roll.data*pi/180; 
    case 'rad'
        phi  = flightdata.arinc.Ahrs1.Roll.data;
    otherwise
        error('Unit of the roll angle is unknown')
end
DataStruct.Measurement.phi.Sensor = 'arinc/Ahrs1/Roll';

% Obtain theta (pitch angle)
switch flightdata.arinc.Ahrs1.Pitch.Unit
    case 'deg'
        theta  = flightdata.arinc.Ahrs1.Pitch.data*pi/180; 
    case 'rad'
        theta  = flightdata.arinc.Ahrs1.Pitch.data;
    otherwise
        error('Unit of the pitch angle is unknown')
end
DataStruct.Measurement.theta.Sensor = 'arinc/Ahrs1/Pitch';

% Obtain psi (yaw angle)
switch flightdata.arinc.Ahrs1.magnHdg.Unit
    case 'deg'
        psi  = flightdata.arinc.Ahrs1.magnHdg.data*pi/180; 
    case 'rad'
        psi  = flightdata.arinc.Ahrs1.magnHdg.data;
    otherwise
        error('Unit of the yaw angle is unknown')
end
DataStruct.Measurement.psi.Sensor = 'arinc/Ahrs1/magnHdg';

% Obtain Vtas (true airspeed)
switch flightdata.arinc.Dadc1.tas.Unit
    case 'knots'
        Vtas  = flightdata.arinc.Dadc1.tas.data*0.514444; 
    case 'm/s'
        Vtas  = flightdata.arinc.Dadc1.tas.data;
    otherwise
        error('Unit of the true airspeed is unknown')
end
DataStruct.Measurement.Vtas.Sensor = 'arinc/Dadc1/tas';

% Obtain alpha (angle of attack)
switch flightdata.synchro.boom_alpha.Unit
    case 'deg'
        alpha   = flightdata.synchro.boom_alpha.data*pi/180; 
        alpha_L = flightdata.synchro.boom_alpha.data*pi/180;
        alpha_R = flightdata.synchro.boom_alpha.data*pi/180;
    case 'rad'
        alpha   = flightdata.synchro.boom_alpha.data;
        alpha_L = flightdata.synchro.boom_alpha.data;
        alpha_R = flightdata.synchro.boom_alpha.data;
    otherwise
        error('Unit of the angle of attack is unknown')
end
DataStruct.Measurement.alpha.Sensor     = 'synchro/boom_alpha';
DataStruct.Measurement.alpha_L.Sensor   = 'synchro/boom_alpha';
DataStruct.Measurement.alpha_R.Sensor   = 'synchro/boom_alpha';

% Obtain beta (angle of side slip)
switch flightdata.synchro.boom_beta.Unit
    case 'deg'
        beta  = flightdata.synchro.boom_beta.data*pi/180; 
    case 'rad'
        beta  = flightdata.synchro.boom_beta.data;
    otherwise
        error('Unit of the angle of sideslip is unknown')
end
DataStruct.Measurement.beta.Sensor = 'synchro/boom_beta';    
%% Define measurement matrix
y = [xpos; ypos; alt; xdot; ydot; hdot; phi; theta; psi; Vtas; alpha; beta; alpha_L; alpha_R];
end