function [x_k_1k_1, P_k_1k_1, z_k1_k0, Ve] = EKF(f, h, G, x_k0_k0, P_k0_k0, u, Z_k, Q, R, dt)
% EKF.m performs a single step of the Extended Kalman Filter
%
%          step 1: One-step ahead prediction
%          step 2: Calculate the jacobian and the input noise matrix
%          step 3: Discretize state transition & input matrix
%          step 4: Covariance matrix of state prediction error
%          step 5: Kalman gain calculation
%          step 6: Measurement update 
%          step 7: Covariance matrix of state estimation error
%
% Inputs:   f:          function handle for f(x)
%           h:          function handle for h(x)
%           G:          function handle for G (input matrix)
%           x_k0_k0:    "a priori" state estimate
%           P_k0_k0:    "a priori" estimated state covariance
%           u:          Input of the system
%           Z_k:        current measurement
%           Q:          process noise covariance 
%           R:          measurement noise covariance
%           dt:         rate at which the data was sampled in seconds
%
% Outputs:  x_k_1k_1:   "a posteriori" state estimate
%           P_k_1k_1:   "a posteriori" state covariance
%           z_k1_k0:    "a posteriori" measurement estimate
%           Ve:         Pe(k+1|k) (covariance matrix of innovation)
%          
%
% REFERENCES: 
% C.C. de Visser, AE4320 SystemIdentification of Aerospace Vehicles - Parameter Estimation, (2015).
%
% Made by: M.A. van den Hoek & L.J. van Horssen, August 2016 - Version 1.0
%% Run the Extended Kalman filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the number of states 
n = length(x_k0_k0);
%% Step 1: One step ahead prediction
   
% Predicted state: x(k+1|k)
[~, x_k1_k0] = rk4(f, x_k0_k0', u, zeros(6,1), [0 dt]);
   
% Predicted output: z(k+1|k)
z_k1_k0 = h(0,x_k1_k0, u, zeros(12,1));

%% Step 2: Calculate the jacobian matrix 
Fx = numjacobian(f,x_k1_k0);   
    
% Calculate the Jacobian Hx
Hx = numjacobian(h,x_k1_k0);
%% Step 3: Discretize state transition & input matrix
[Phi, Gamma] = c2d(Fx, G(x_k0_k0), dt);

%% Step 4: Covariance matrix of state prediction error 
P_k1_k0 = Phi*P_k0_k0*Phi' + Gamma*Q*Gamma'; 
        
%% Step 5: Kalman gain calculation    
% Pz(k+1|k) (covariance matrix of innovation)
Ve = (Hx*P_k1_k0 * Hx' + R);   
        
% K(k+1) (Kalman gain)
K = P_k1_k0 * Hx' / Ve;
              
%% Step 6: Measurement update   
% x(k+1|k+1) = x(k+1|k) + K*(z(k) - h(x(k+1|k),u(k)))
x_k_1k_1 = x_k1_k0' + K * (Z_k - z_k1_k0');
            
%% Step 7: Covariance matrix of state estimation error
P_k_1k_1 = (eye(n) - K*Hx) * P_k1_k0;
end

function [t,w] = rk4(fn, xin, uin, win, t)
% Calculates the output of a Runga-Kutta fourth order integration method
    a = t(1); 
    b = t(2);
    w = xin;
    N = 2;
    h = (b-a) / N;
    t = a;

    for j=1:N
        K1 = h * fn(t, w, uin, win);
        K2 = h * fn(t+h/2, w+K1/2, uin, win);
        K3 = h * fn(t+h/2, w+K2/2, uin, win);
        K4 = h * fn(t+h, w+K3, uin, win);

        w = w + (K1 + 2*K2 + 2*K3 + K4) / 6;
        t = a + j*h;
    end
end