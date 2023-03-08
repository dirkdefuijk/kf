function [u] = FilterInputs(u, dt, fc, fo, filter_type)
% FilterInputs.m is a script which pre-filters the inputs if necessary. 
%
% Inputs:   
%           u:              input matrix        [6xM matrix]
%           dt:             sampling rate       [s]
%           fc:             Cut-off frequency   [Hz]
%           fo:             Filter order        [-]
%           filter_type:    Type of the filter that is being used, e.g. Butterworth
%
%
% Outputs:  
%           u:              filtered input matrix [6xM matrix]
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Low pass filter the linear accelerations and rotational rates
% Obtain the linear accelerations and rotational rates from the input matrix
ax = u(1,:); ay = u(2,:); az = u(3,:); 
p  = u(4,:); q  = u(5,:); r  = u(6,:); 

% Low pass filter the linear accelerations
ax = LPfilter(ax, fc, fo, dt, filter_type);
ay = LPfilter(ay, fc, fo, dt, filter_type);
az = LPfilter(az, fc, fo, dt, filter_type);

% Low pass filter the rotational rates
p = LPfilter(p, fc, fo, dt, filter_type);
q = LPfilter(q, fc, fo, dt, filter_type);
r = LPfilter(r, fc, fo, dt, filter_type);
    
% re-define input matrix
u = [ax; ay; az; p; q; r];
end