function [x0] = DefineInitialCondition(y)
% DefineInitialCondition.m is a script in with which the initial condition
% can be defined, where N is the number of states. The initial condition 
% can be based on the initial measurements (y) or be determined manually.
% A maximum of three initial conditions can be set, for determining the
% convergence of the filter.
%
% Inputs:   
%           y: measurement
%
% Outputs:  
%           x0: Initial condition [Nx1 vector]
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Set initial condition for the Kalman filter
x0 = ([y(1,1); ...                 % x
       y(2,1); ...                 % y
       -y(3,1); ...                % z
       y(10,1)*cos(y(11,1)); ...   % u
       0; ...                      % v
       y(10,1)*sin(y(11,1)); ...   % w
       y(7,1); ...                 % Phi
       y(8,1);...                  % Theta
       y(9,1); ...                 % Psi
       0; ...                      % Lx
       0; ...                      % Ly
       0; ...                      % Lz
       0; ...                      % Lp
       0; ...                      % Lq
       0; ...                      % Lr
       0; ...                      % Wxe
       0; ...                      % Wye
       0.2; ...                    % Caup
       y(10,1)*cos(y(11,1)); ...   % u_L
       y(10,1)*cos(y(11,1)); ...   % u_R
       0; ...                      % v_L
       0; ...                      % v_R
       y(10,1)*sin(y(11,1)); ...   % w_L
       y(10,1)*sin(y(11,1)); ...   % w_R
       ]);

%% Example of using three initial conditions
% x0 = ([y(1,1)                  , y(1,1)                  , y(1,1)              ; ...   % x
%        y(2,1)                  , y(2,1)                  , y(2,1)              ; ...   % y
%        -y(3,1)                 , -y(3,1)                 , -y(3,1)             ; ...   % z
%        y(10,1)*cos(y(11,1))    , y(10,1)*cos(y(11,1))    , y(10,1)*cos(y(11,1)); ...   % u
%        0                       , 0                       , 0                   ; ...   % v
%        y(10,1)*sin(y(11,1))    , y(10,1)*sin(y(11,1))    , y(10,1)*sin(y(11,1)); ...   % w
%        y(7,1)                  , y(7,1)                  , y(7,1)              ; ...   % Phi
%        y(8,1)                  , y(8,1)                  , y(8,1)              ; ...   % Theta
%        y(9,1)                  , y(9,1)                  , y(9,1)              ; ...   % Psi
%        0                       , 1                       , -1                  ; ...   % Lx
%        0                       , 1                       , -1                  ; ...   % Ly
%        0                       , 1                       , -1                  ; ...   % Lz
%        0                       , 1                       , -1                  ; ...   % Lp
%        0                       , 1                       , -1                  ; ...   % Lq
%        0                       , 1                       , -1                  ; ...   % Lr
%        0                       , 10                      , -10                 ; ...   % Wxe
%        0                       , 10                      , -10                 ; ...   % Wye
%        0.2                     , 1                       , -1                  ; ...   % Caup
%        ]);
end