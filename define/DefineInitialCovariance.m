function [P0] = DefineInitialCovariance()
% DefineInitialCovariance.m is a script in with which the initial
% covariance matrix can be defined. N is the number of states.
%
% Inputs:   
%           none
%
% Outputs:  
%           P0: Initial condition [NxN matrix]
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Set initial covariance matrix for the Kalman filter
P0 = diag([10000, ...           % x
           10000, ...           % y
           10000, ...           % h
           10, ...              % u
           10, ...              % v
           10, ...              % w
           10, ...              % phi
           10, ...              % theta
           10, ...              % psi
           1, ...               % Lx
           1, ...               % Ly
           1, ...               % Lz
           1, ...               % Lp
           1, ...               % Lq
           1, ...               % Lr
           10, ...              % Wxe
           10, ...              % Wye
           1, ...               % Caup
           10, ...              % u_L
           10, ...              % u_R
           10, ...              % v_L
           10, ...              % v_R
           10, ...              % w_L
           10, ...              % w_R
           ]); 
end