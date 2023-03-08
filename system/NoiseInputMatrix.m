% Filename :  NoiseInputMatrix.m
%
% Constructs the noise input matrix for the state space system
%
% Inputs: 
%     xin:      state vector x:
%               x = [x, y, z, u, v, w, phi, theta, psi, Wxe, Wye, Caup,
%                    lambda_x, lambda_y, lambda_z, lambda_p, lambda_q, 
%                    lambda_r, alpha_v]^T
% Outputs:
%     G:        noise input matrix        
%
% Made by: M.A. van den Hoek and L.J. van Horssen, September 2016
function G = NoiseInputMatrix(xin)
%% Extract the state and inputs and put them into symbolic variables
x = xin(1); y = xin(2); zin = xin(3);
u = xin(4); v = xin(5); w = xin(6);
phi = xin(7); theta = xin(8); psi = xin(9);
Wxe = xin(10); Wye = xin(11);
Caup = xin(12); 
Lx = xin(13); Ly = xin(14); Lz = xin(15);
Lp = xin(16); Lq = xin(17); Lr = xin(18); 
alpha_v = xin(19);

%%
G = ...
[0,  0,  0,  0,                    0,                    0;     % x
 0,  0,  0,  0,                    0,                    0;     % y
 0,  0,  0,  0,                    0,                    0;     % z
-1,  0,  0,  0,                    w,                    v;     % u
 0, -1,  0, -w,                    0,                    u;     % v
 0,  0, -1,  v,                   -u,                    0;     % w
 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta);     % phi
 0,  0,  0,  0,            -cos(phi),             sin(phi);     % theta
 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta);     % psi
 0,  0,  0,  0,                    0,                    0;     % Wxe
 0,  0,  0,  0,                    0,                    0;     % Wye
 0,  0,  0,  0,                    0,                    0;     % Caup
 0,  0,  0,  0,                    0,                    0;     % Lx
 0,  0,  0,  0,                    0,                    0;     % Ly
 0,  0,  0,  0,                    0,                    0;     % Lz
 0,  0,  0,  0,                    0,                    0;     % Lp
 0,  0,  0,  0,                    0,                    0;     % Lq
 0,  0,  0,  0,                    0,                    0;     % Lr
 0,  0,  0,  0,                    0,                    0];    % alpha_v