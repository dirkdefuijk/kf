function [u] = Correct_cg_loc(u, massmodel, dt, filter, fc, fo, filter_type)
% Correct_cg_loc.m is a script wich corrects the measurement error due to
% the measurement of the accelerations not being in the center of gravity
%
% Inputs:   
%           u:              input matrix        [6xM matrix]
%           massmodel:      mass model data structure containing c.g. location
%           filter:         string with either 'on' or 'off' containing
%                           information on whether or not the accelerations
%                           already have been pre-filtered before
%           fc:             Cut-off frequency [Hz]
%           fo:             Filter order        [-]
%           filter_type:    Type of the filter that is being used, e.g. Butterworth
%
%
% Outputs:  
%           u:              filtered input matrix [6xM matrix], with
%                           corrected acceleration measurements
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Correct acceleration measurements for offset in c.g. position and acceleration measurement
% Obtain the linear accelerations and the rotational rates
ax = u(1,:); ay = u(2,:); az = u(3,:);
p = u(4,:);  q = u(5,:);  r = u(6,:);

% Datum is 19 inch infront of the nose and 91 under the nose
% IMU system is approx. 75 behind the nose and at same height (ASSUMPTION!! NOT MEASURED!!)

% 2.3876 m = 19 + 75 inches (ASSUMPTION!! NOT MEASURED!!)
dx = massmodel.x_cg.data - 2.3876; 

% IMU is not perfectly centered, so a small offset is used. (ASSUMPTION!! NOT MEASURED!!)
dy = massmodel.y_cg.data - 0.1; 

% 2.3114 m = 91 inches (ASSUMPTION!! NOT MEASURED!!)
dz = massmodel.z_cg.data - 2.3114; 

% In case the accelerations and rotational rates are already pre-filtered
% the value of filter is 'on' otherwise it is 'off'
if strcmp(filter,'off')
    % Filter linear accelerations in case it has not been done before
%     ax = LPfilter(ax, fc, fo, dt, filter_type);
%     ay = LPfilter(ay, fc, fo, dt, filter_type);
%     az = LPfilter(az, fc, fo, dt, filter_type);

    % Filter rotational rates in case it has not been done before
    p = LPfilter(p, fc, fo, dt, filter_type);
    q = LPfilter(q, fc, fo, dt, filter_type);
    r = LPfilter(r, fc, fo, dt, filter_type);
    
    % Calculate the rotational accelerations
    pdot = diff(p)/dt; qdot = diff(q)/dt; rdot = diff(r)/dt;
    
    % Add leading zero to obtain correct vector size
    pdot = [0 pdot]; qdot = [0 qdot]; rdot = [0 rdot];
    
elseif strcmp(filter,'on')
    % Calculate the rotational accelerations
    pdot = diff(p)/dt; qdot = diff(q)/dt; rdot = diff(r)/dt;
    
    % Add leading zero to obtain correct vector size
    pdot = [0 pdot]; qdot = [0 qdot]; rdot = [0 rdot];
else
    error('Unknown input for filter, use either "on" or "off"')
end  

%% Do the correction for the offset in acceleration measurement (i.e. IMU position) and c.g. location
Ax_cg = ax + dx.*(q.^2 + r.^2) - dy.*(p.*q - rdot) - dz.*(p.*r + qdot);  %ax
Ay_cg = ay + dy.*(r.^2 + p.^2) - dz.*(q.*r - pdot) - dx.*(q.*p + rdot);  %ay
Az_cg = az + dz.*(p.^2 + q.^2) - dx.*(r.*p - qdot) - dy.*(r.*q + pdot);  %az

% Store the accelerations in the correct position in the input matrix
u(1,:) = Ax_cg;
u(2,:) = Ay_cg;
u(3,:) = Az_cg;