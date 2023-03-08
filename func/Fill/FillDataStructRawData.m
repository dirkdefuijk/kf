% FillDataStructRawData.m is a script which is used fill the raw
% measurements in the datastructure. This includes all measurements such as
% the euler angle but also the control deflections and flap settings.
%
% Inputs:   
%           none
%
% Outputs:  
%           none (filled data structure with raw measurement data)
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%%
% Define the gravitational acceleration on Earth
g0 = 9.80665;

% Create exitflag in case no information is available for the flap setting
exitflag_flaps = 0;

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
        error('Unit of the pitch rate is unknown')
end

% Store information on acceleration in x-direction
DataStruct.Input.u.ax.time = t;
DataStruct.Input.u.ax.raw = ax;
DataStruct.Input.u.ax.Sensor = 'arinc/Ahrs1/bLongAcc';

% Store information on acceleration in y-direction
DataStruct.Input.u.ay.time = t;
DataStruct.Input.u.ay.raw = ay;
DataStruct.Input.u.ay.Sensor = 'arinc/Ahrs1/bLatAcc';

% Store information on acceleration in z-direction
DataStruct.Input.u.az.time = t;
DataStruct.Input.u.az.raw = az;
DataStruct.Input.u.az.Sensor = 'arinc/Ahrs1/bNormAcc';

% Store information on roll rate
DataStruct.Input.u.p.time = t;
DataStruct.Input.u.p.raw = p;
DataStruct.Input.u.p.Sensor = 'arinc/Ahrs1/bRollRate';

% Store information on pitch rate
DataStruct.Input.u.q.time = t;
DataStruct.Input.u.q.raw = q;
DataStruct.Input.u.q.Sensor = 'arinc/Ahrs1/bPitchRate';

% Store information on yaw rate
DataStruct.Input.u.r.time = t;
DataStruct.Input.u.r.raw = r;
DataStruct.Input.u.r.Sensor = 'arinc/Ahrs1/bYawRate';

%% Control deflections
% Obtain elevator deflection
switch flightdata.synchro.delta_e.Unit
    case 'deg'
        delta_e = flightdata.synchro.delta_e.data*pi/180;
    case 'rad'
        delta_e = flightdata.synchro.delta_e.data;
    otherwise
        error('Unit of the elevator deflection is unknown')
end

% Obtain aileron deflection
switch flightdata.synchro.delta_a.Unit
    case 'deg'
        delta_a = flightdata.synchro.delta_a.data*pi/180;
    case 'rad'
        delta_a = flightdata.synchro.delta_a.data;
    otherwise
        error('Unit of the aileron deflection is unknown')
end

% Obtain rudder deflection
switch flightdata.synchro.delta_r.Unit
    case 'deg'
        delta_r = flightdata.synchro.delta_r.data*pi/180;
    case 'rad'
        delta_r = flightdata.synchro.delta_r.data;
    otherwise
        error('Unit of the rudder deflection is unknown')
end

% Obtain flap deflection
switch flightdata.analog.flaps.Unit
    case 'deg'
        flaps = flightdata.analog.flaps.data*pi/180;
    case 'rad'
        flaps = flightdata.analog.flaps.data;
    otherwise
        exitflag_flaps = 1;
end

% Store information on elevator deflection
DataStruct.Input.delta.de.time = t;
DataStruct.Input.delta.de.raw = delta_e;
DataStruct.Input.delta.de.Sensor = 'synchro/delta_e';

% Store information on aileron deflection
DataStruct.Input.delta.da.time = t;
DataStruct.Input.delta.da.raw = delta_a;
DataStruct.Input.delta.da.Sensor = 'synchro/delta_a';

% Store information on rudder deflection
DataStruct.Input.delta.dr.time = t;
DataStruct.Input.delta.dr.raw = delta_r;
DataStruct.Input.delta.dr.Sensor = 'synchro/delta_r';

% Store information on flap deflection in case information is available.
% Otherwise leave the flap setting empty
if exitflag_flaps == 1
    DataStruct.Input.delta.flaps.time = [];
    DataStruct.Input.delta.flaps.raw = [];
else
    DataStruct.Input.delta.flaps.time = t;
    DataStruct.Input.delta.flaps.raw = flaps*180/pi;
    DataStruct.Input.delta.flaps.Sensor = 'analog/flaps';
end


%% Observation measurements
% Store information on the x-postion in the ECEF frame
DataStruct.Measurement.x.raw = meas(1,:);

% Store information on the y-postion in the ECEF frame
DataStruct.Measurement.y.raw = meas(2,:);

% Store information on the altitude
DataStruct.Measurement.h.raw = meas(3,:);

% Store information on the velocity in x-direction in the ECEF frame
DataStruct.Measurement.xdot.raw = meas(4,:);

% Store information on the velocity in y-direction in the ECEF frame
DataStruct.Measurement.ydot.raw = meas(5,:);

% Store information on the climb rate
DataStruct.Measurement.hdot.raw = meas(6,:);

% Store information on the roll angle
DataStruct.Measurement.phi.raw = meas(7,:);

% Store information on the pitch angle
DataStruct.Measurement.theta.raw = meas(8,:);

% Store information on the yaw angle
DataStruct.Measurement.psi.raw = meas(9,:);

% Store information on the true airspeed
DataStruct.Measurement.Vtas.raw = meas(10,:);

% Store information on the angle of attack
DataStruct.Measurement.alpha.raw = meas(11,:);

% Store information on the angle of attack - left wing
DataStruct.Measurement.alpha_L.raw = meas(11,:);

% Store information on the angle of attack - right wing
DataStruct.Measurement.alpha_R.raw = meas(11,:);

% Store information on the angle of side slip
DataStruct.Measurement.beta.raw = meas(12,:);