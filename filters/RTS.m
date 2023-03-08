function [xs, Ps, ys, Pyy] = RTS(f, h, G, x, P, u, Q, R, dt)
% RTS.m performs the Rauch-Tung-Striebel smoother
%
%          step 1: Calculate one-step ahead prediction
%          step 2: Calculate the jacobian matrices 
%          step 3: Discretize state transition & input matrix
%          step 4: Covariance matrix of state prediction error 
%          step 5: Calculate the smoother gain matrix, the smoothed state and the smoothed covariance matrix
%          step 6: Calculated smoothed measurement and smoothed measurement covariance matrix
%
% Inputs:   f:      function handle for f(x)
%           h:      function handle for h(x)
%           G:          function handle for G (input matrix)
%           x:      "a priori" state estimate
%           P:      "a priori" estimated state covariance
%           u:      Input of the system
%           Q:      process noise covariance 
%           R:      measurement noise covariance
%           dt:     rate at which the data was sampled in seconds
%
% Outputs:  xs:     smoothed state estimate
%           Ps:     smoothed state covariance matrix
%           ys:     smoothed measurement estimate
%           Pyy:    smoothed measurement covariance matrix
%
% REFERENCES: 
% Haykin, S. S. (Ed.). (2001). Kalman Filtering and Neural Networks (pp. 17). New York: Wiley.
%
% Made by: M.A. van den Hoek & L.J. van Horssen, August 2016 - Version 1.0
%% Run the Rauch-Tung-Striebel Smoother

addpath([pwd '/func/ConsoleProgressBar']) 

% Open console progress bar and set options
% Open a progress bar   
cpb = ConsoleProgressBar();
    
 
% Set progress bar parameters
cpb.setLeftMargin(4);   % progress bar left margin
cpb.setTopMargin(1);    % rows margin

cpb.setLength(40);      % progress bar length: [.....]
cpb.setMinimum(0);      % minimum value of progress range [min max]
cpb.setMaximum(length(x));      % maximum value of progress range [min max]

cpb.setElapsedTimeVisible(1);
cpb.setRemainedTimeVisible(1);

cpb.setElapsedTimePosition('left');
cpb.setRemainedTimePosition('right');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Rauch-Tung-Striebel Smoother
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the smoothed state and smoothed state covariance matrices
xs = x; % Smoothed state
Ps = P; % Smoothed state covariance matrix

% Start console progress bar
cpb.start();
        
for k = length(x)-1:-1:1
    %% Step 1: Calculate one-step ahead prediction
    % Calculate the one-step ahead prediction
    [~, m_kp1] = rk4(f, x(:,k)', u(:,k)', zeros(6,1), [0 dt]);

    %% Step 2: Calculate the jacobian matrices 
    Fx = numjacobian(f,x(:,k));   

    % Calculate the Jacobian Hx
    Hx = numjacobian(h,x(:,k));

    %% Step 3: Discretize state transition & input matrix
    [Phi, Gamma] = c2d(Fx, G(x(:,k)), dt);

    %% Step 4: Covariance matrix of state prediction error 
    P_kp1 = Phi*P(:,:,k)*Phi' + Gamma*Q*Gamma'; 

    %% Step 5: Calculate the smoother gain matrix, the smoothed state and the smoothed covariance matrix
    % Calculate the smoother gain matrix (A)
    A = P(:,:,k) * Phi' * inv(P_kp1);
    
    % Calculate the smoothed state
    xs(:,k) = x(:,k) + A*(xs(:,k+1) - m_kp1');
    
    % Calculate the smoothed state covariance matrix
    Ps(:,:,k) = P(:,:,k) + A*(Ps(:,:,k+1) - P_kp1)*A';
   
    % Print the current data point that is being filtered
    if mod(k,100) == 0
        text = sprintf('Smoothing the states: %d/%d', k, length(x));

        cpb.setValue(k);
        cpb.setText(text);
    end   
end

% Stop console progress bar
cpb.stop();
        
%% Step 6: Calculated smoothed measurement and smoothed measurement covariance matrix

% Pre-allocate memory
ys = zeros(size(R,1),length(xs));
Pyy = zeros(size(R,1),size(R,2),length(xs));

% Start console progress bar
cpb.start();

for k = 1:length(xs)
    
    % Calculate the smoothed measurements
    ys(:,k) = h(0,xs(:,k),u(:,k),zeros(12,1));
    
    % Calculate the smoothed measurement covariance matrix
    Pyy(:,:,k) = (Hx * Ps(:,:,k) * Hx' + R);  
       
    % Print the current data point that is being filtered
    if mod(k,100) == 0
        text = sprintf('Smoothing the measurements: %d/%d', k, length(x));

        cpb.setValue(k);
        cpb.setText(text);
    end   
end

% Stop console progress bar
cpb.stop();
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