function [x_kp1, P_kp1, y_k, P_ykyk] = UKF(f, h, x, P, u, y, Q, R, dt, varargin)
% UKF.m performs a single step of the Unscented Kalman Filter
%
%          step 1: Augment the system
%          step 2: Calculate Sigma points and correspondig weights
%          step 3: Transform Sigma points through system equations
%          step 4: Transform Sigma points through measurements equations
%          step 5: Measurement update equations
%
% Inputs:   f:          function handle for f(x)
%           h:          function handle for h(x)
%           x:          "a priori" state estimate
%           P:          "a priori" estimated state covariance
%           u:          Input of the system
%           y:          current measurement
%           Q:          process noise covariance 
%           R:          measurement noise covariance
%           dt:         rate at which the data was sampled in seconds
%
% Outputs:  x_kp1:      "a posteriori" state estimate
%           P_kp1:      "a posteriori" state covariance
%           y_k:        "a posteriori" measurement estimate
%           P_ykyk:     "a posteriori" measurement covariance
%
%
% REFERENCES: 
% Wan, E., & Van Der Merwe, R. (2000). The unscented Kalman filter for
% nonlinear estimation. In Adaptive Systems for Signal Processing,
% Communications, and Control Symposium 2000. AS-SPCC. The IEEE 2000 
% (pp. 153-158). IEEE.
%
% Learning the Unscented Kalman Filter - Yi Cao (12 Dec 2010)
% link: http://www.mathwo(ks.com/matlabcentral/fileexchange/18217-learning-the-unscented-kalman-filter/content/ukf.m
% Accessed at 27 October 2015
%
% Made by: M.A. van den Hoek & L.J. van Horssen, August 2016 - Version 1.0
%% Tuning parameters UKF
% Set scaling parameter eta (in article named alpha), which determines the spread of the sigma
% points around the state x, default: 1e-3
alpha = 1e-3;
 
% Set the secondary scaling parameter kappa, default: 0
kappa = 0;
 
% Set scaling factor (in article named beta) to incorportate prior knowledge of the distrubtion of
% x (for Gaussian distributions), default: 2
beta = 2;

% Change values of alpha, beta and/or kappa if additional inputs are given
switch nargin
    case 11
        if strcmp(varargin(1),'alpha')
            alpha = varargin{2};
        elseif strcmp(varargin(1),'kappa')
            kappa = varargin{2};
        elseif strcmp(varargin(1),'beta')
            beta = varargin{2};
        else
            error('Unknown input type. Make sure that additional input only contain: "alpha", "beta" or "kappa"!')
        end
    case 13
        for i = 1:2:length(varargin)
            if strcmp(varargin(i),'alpha')
                alpha = varargin{i+1};
            elseif strcmp(varargin(i),'kappa')
                kappa = varargin{i+1};
            elseif strcmp(varargin(i),'beta')
                beta = varargin{i+1};
            else
                error('Unknown input type. Make sure that additional input only contain: "alpha", "beta" or "kappa"!')
            end
        end
    case 15
        for i = 1:2:length(varargin)
            if strcmp(varargin(i),'alpha')
                alpha = varargin{i+1};
            elseif strcmp(varargin(i),'kappa')
                kappa = varargin{i+1};
            elseif strcmp(varargin(i),'beta')
                beta = varargin{i+1};
            else
                error('Unknown input type. Make sure that additional input only contain: "alpha", "beta" or "kappa"!')
            end
        end
    otherwise
        if nargin > 15
            error('Number of inputs is too large!')
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the the Unscented Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Augment the system
x_a = [x', zeros(1,length(Q)), zeros(1,length(R))]';

P_a = [P                           , zeros(length(P),length(Q)) , zeros(length(P),length(R)); ...
        zeros(length(Q), length(P)), Q                           , zeros(length(Q),length(R)); ...
        zeros(length(R), length(P)), zeros(length(R),length(Q)) , R];
   
%% Determine L, lambda and the dimensions of the input noise and measurement noise
% Set the value for L, which is the dimension of the state x
L = length(x_a);

% Calculate lambda, which is a scaling parameter
lambda = alpha^2 * (L + kappa) - L;

% Dimension of the input noise
Nv = length(Q);

% Dimension of the measurement noise
Nn = length(R);
    
%% Step 2: Calculate Sigma points and correspondig weights
% Calculate the value of L + lambda
gamma = L + lambda;

% Calculate the sigma points
X_a = sigmas(x_a,P_a,sqrt(gamma));

% Calculate the correspondig weights to the sigma points
Wm=[lambda/gamma 0.5/gamma+zeros(1,2*length(x_a))];  %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);                         %weights for covariance

% Select sigma points 
Xx = X_a(1:L-Nv-Nn,1:2*length(x_a)+1)'; %Sigma points correspondig to state estimate
Xv = X_a(L - Nv - Nn + 1:L - Nn,1:2*length(x_a)+1)'; %Sigma points corresponding to process noise
Xn = X_a(L - Nn + 1:L,1:2*length(x_a)+1)'; %Sigma points corresponding to measurement noise

%% Step 3: Transform Sigma points through system equations
x_k = 0;
X_k = zeros(2*length(x_a)+1,length(x));
for j = 1:2*length(x_a)+1
    % Calculate transformed sampling points
    [~, X_k(j,:)] = rk4(f, Xx(j,:), u, Xv(j,:), [0 dt]);
end

% Calculate transformed mean
x_k = Wm*X_k;

% Calculate transformed deviations
X1 = X_k-x_k(ones(1,2*length(x_a)+1),:);

% Calculate transformed covariance
P_k = X1'*diag(Wc)*X1;

%% Step 4: Transform Sigma points through measurements equations
y_k = 0;
Y_k = zeros(2*length(x_a)+1,length(y));
for j = 1:2*length(x_a)+1
    % Calculate transformed sampling points
    Y_k(j,:) = h(0,X_k(j,:), u, Xn(j,:)); 

    % Calculate transformed mean
    y_k = y_k + Wm(j)*Y_k(j,:);
end

% Calculate transformed deviations
Y1 = Y_k-y_k(ones(1,2*length(x_a)+1),:);

%% Step 5: Measurement update equations
% Calculate the transformed covariance of the output
P_ykyk  = Y1'*diag(Wc)*Y1; 

% Calculate the transformed cross-covariance 
P_xkyk  = X1'*diag(Wc)*Y1; 

% Calculate the Kalman Gain
K = P_xkyk/P_ykyk;

% Calculate the one step ahead prediction
x_kp1 = x_k'+K*(y-y_k');

% Calculate the one step ahead covariance matrix
P_kp1 = P_k - K*P_ykyk*K';
end

function X=sigmas(x,P,c)
% Sigmas calculates the sigma points around a reference point
%
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];    
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