function [DataStruct] = AddGravityAcc(DataStruct)
% AddGravityAcc.m is a script which adds the gravity component to the
% accelerations in y- and z-direction. The acceleration in x-direction
% has not been corrected for the gravity and thus already contains a
% gravity component.
%
% Inputs:   
%           DataStruct: data structure containing y- and z- accelerations
%           without a gravity component
%
% Outputs:  
%           DataStruct: data structure containing the correct accelerations
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Correct the raw data accelerations

% Define the gravitational acceleration on Earth
g0 = 9.80665;

% Obtain the "incorrect" accelerations
ax = DataStruct.Input.u.ax.raw;
ay = DataStruct.Input.u.ay.raw;
az = DataStruct.Input.u.az.raw;

% Obtain the Euler angles
phi = DataStruct.State.phi.data;
theta = DataStruct.State.theta.data;
psi = DataStruct.State.psi.data;

% Add gravity component in the y- and z-directions for each data point
for i = 1:length(ax)
    % Calculate the transformation matrix from ECEF to body frame
    Tbe = TbE(phi(i),theta(i),psi(i));
    
    % Calculate gravity component in the body frame
    g = -Tbe*[0; 0; g0];
    
    % Add gravity component accelerations in y- and z-directions
    ax(i) = ax(i);  % No correction for Ax is needed
    ay(i) = ay(i);  % No correction for Ay is needed
    az(i) = az(i) - g0;
end

% Store accelerations in the data structure
DataStruct.Input.u.ax.raw = ax;
DataStruct.Input.u.ay.raw = ay;
DataStruct.Input.u.az.raw = az;

%% Correct the filtered data accelerations

% Obtain the "incorrect" accelerations
ax = DataStruct.Input.u.ax.filtered;
ay = DataStruct.Input.u.ay.filtered;
az = DataStruct.Input.u.az.filtered;

% Obtain the Euler angles
phi = DataStruct.State.phi.data;
theta = DataStruct.State.theta.data;
psi = DataStruct.State.psi.data;

% Add gravity component in the y- and z-directions for each data point
for i = 1:length(ax)
    % Calculate the transformation matrix from ECEF to body frame
    Tbe = TbE(phi(i),theta(i),psi(i));
    
    % Calculate gravity component in the body frame
    g = -Tbe*[0; 0; g0];
    
    % Add gravity component accelerations in y- and z-directions
    ax(i) = ax(i);  % No correction for Ax is needed
    ay(i) = ay(i);  % No correction for Ax is needed
    az(i) = az(i) - g0;
end

% Store accelerations in the data structure
DataStruct.Input.u.ax.filtered = ax;
DataStruct.Input.u.ay.filtered = ay;
DataStruct.Input.u.az.filtered = az;