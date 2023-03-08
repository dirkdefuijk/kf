function [x_k_1k_1, P_k_1k_1, z_k1_k0, Ve] = IEKF(f, h, G, x_k0_k0, P_k0_k0, u, Z_k, Q, R, dt, varargin)
% IEKF.m performs a single step of the Iterated Extended Kalman Filter
%
%          step 1: One-step ahead prediction
%          step 2: Calculate the jacobian and the input noise matrix
%          step 3: Discretize state transition & input matrix
%          step 4: Covariance matrix of state prediction error
%          step 5: Measurement equation Jacobian recalculation 
%          step 6: Kalman gain calculation
%          step 7: Measurement update 
%          step 8: Covariance matrix of state estimation error
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
%% Tuning parameters IEKF

% Maximum error before stopping iterations of the IEKF
epsilon = 1e-10;

% Set maximum number of iterations for the IEKF
maxIterations = 100;

% Change values of epsilon and maxIterations if additional inputs are given
switch nargin
    case 12
        if strcmp(varargin(1),'epsilon')
            epsilon = varargin{2};
        elseif strcmp(varargin(1),'maxIterations')
            maxIterations = varargin{2};
        else
            error('Unknown input type. Make sure that additional input only contain: "epsilon", "maxIterations"!')
        end
    case 14
        for i = 1:2:length(varargin)
            if strcmp(varargin(i),'epsilon')
                epsilon = varargin{i+1};
            elseif strcmp(varargin(i),'maxIterations')
                maxIterations = varargin{i+1};
            else
                error('Unknown input type. Make sure that additional input only contain: "epsilon", "maxIterations"!')
            end
        end
    otherwise
        if nargin > 14
            error('Number of inputs is too large!')
        end
end

% Determine the number of states 
n = length(x_k0_k0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Iterated Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: One step ahead prediction
   
% Predicted state: x(k+1|k)
[~, x_k1_k0] = rk4(f, x_k0_k0', u, zeros(6,1), [0 dt]);
   
% Predicted output: z(k+1|k)
z_k1_k0 = h(0,x_k1_k0, u, zeros(12,1));

%% Step 2: Calculate the jacobian matrix 
Fx = numjacobian(f,x_k1_k0);   
    
%% Step 3: Discretize state transition & input matrix
[Phi, Gamma] = c2d(Fx, G(x_k0_k0), dt);

%% Step 4: Covariance matrix of state prediction error 
P_k1_k0 = Phi*P_k0_k0*Phi' + Gamma*Q*Gamma'; 
        
%% Step 5: Measurement equation Jacobian recalculation 
    
% Initiliaze eta(i+1) with x(k+1|k)
eta2    = x_k1_k0;
    
% Define initial error
err     = 2*epsilon;
    
% Initiliaze number of iterations
iteration    = 0;
    
% Go through the iterative part while the error is larger than the set
% value of episilon
while (err > epsilon)
    if (iteration >= maxIterations)
        % Break the while loop if the number of maximum iterations is
        % reached and print in the console
        fprintf('Terminating IEKF: exceeded max iterations (%d)\n', maxIterations);
        break
    end      
        
    % Update iteration 
    iteration = iteration + 1;
        
    % Set eta(i+1) equal to eta(i)
    eta1 = eta2;      
        
    % Calculate the Jacobian Hx
    Hx = numjacobian(h,eta1);
    %% Step 6: Kalman gain calculation
              
    % Pz(k+1|k) (covariance matrix of innovation)
    Ve = (Hx*P_k1_k0 * Hx' + R);   
        
    % K(k+1) (Kalman gain)
    K = P_k1_k0 * Hx' / Ve;
        
    % Calculate new observation state
    z_p = h(0,eta1, u, zeros(12,1));
        
    %% Step 7: Measurement update
        
    % Update eta(i+1)
    eta2 = (x_k1_k0' + K * (Z_k - z_p' - Hx*(x_k1_k0 - eta1)'))';
        
    % Calculate the error between eta(i+1) and eta(i)
    err = norm((eta2 - eta1), inf) / norm(eta1, inf);
end
          
% Update the state estimate: x(k+1|k+1), with the last measured update
x_k_1k_1 = eta2';
    
% Calculate the Jacobian Hx
Hx = numjacobian(h,x_k_1k_1');
        
% Pz(k+1|k) (covariance matrix of innovation)
Ve = (Hx*P_k1_k0 * Hx' + R);   
        
% K(k+1) (Kalman gain)
K = P_k1_k0 * Hx' / Ve;
        
% Calculate new observation state
z_p = h(0,x_k_1k_1, u, zeros(12,1));
    
%% Step 8: Covariance matrix of state estimation error
P_k_1k_1 = (eye(n) - K*Hx) * P_k1_k0 * (eye(n) - K*Hx)' + K*R*K';
    
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