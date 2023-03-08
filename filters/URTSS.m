function [xs, Ps, ys, Pyys] = URTSS(f, h, x, P, u, Q, R, dt, varargin)
% URTSS.m performs the Unscented Rauch Tung Striebel Smoother
%
%          step 1: Augment the system
%          step 2: Calculate Sigma points and correspondig weights
%          step 3: Propagate Sigma points through system equations
%          step 4: Compute the predicted mean, covariance and cross-covariance
%          step 5: Compute the smoother Gain, smoothed mean and covariance
%
% Inputs:   f:      function handle for f(x)
%           h:      function handle for h(x)
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
%           Pyys:   smoothed measurement covariance matrix
%
%
% REFERENCES: 
%
% Unscented Rauch-Tung-Striebel Smoother - Simo Särkka
% link: https://users.aalto.fi/~ssarkka/pub/uks-preprint.pdf
% Accessed at 05 September 2016
%
% Made by: M.A. van den Hoek & L.J. van Horssen, August 2016 - Version 1.0
%% Open a console progress bar to show current data point being smoothed
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
    case 10
        if strcmp(varargin(1),'alpha')
            alpha = varargin{2};
        elseif strcmp(varargin(1),'kappa')
            kappa = varargin{2};
        elseif strcmp(varargin(1),'beta')
            beta = varargin{2};
        else
            error('Unknown input type. Make sure that additional input only contain: "alpha", "beta" or "kappa"!')
        end
    case 12
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
    case 14
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
% Run the Unscented Rauch Tung Striebel Smoother
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the smoothed mean and covariance
xs = x;
Ps = P;
ys(:,length(x)) = h(0,x(:,length(x)), u, zeros(12,1)); 

% Start console progress bar
cpb.start();

for k = length(x)-1:-1:1
%% Step 1: Augment the system      
x_a = [x(:,k)', zeros(1,length(Q)), zeros(1,length(R))]';

P_a = [P(:,:,k)                           , zeros(length(P(:,:,k)),length(Q)) , zeros(length(P(:,:,k)),length(R)); ...
        zeros(length(Q), length(P(:,:,k))), Q                                 , zeros(length(Q),length(R)); ...
        zeros(length(R), length(P(:,:,k))), zeros(length(R),length(Q))        , R];
    
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
    
%% Step 3: Propagate Sigma points through system equations
X_k = zeros(2*length(x_a)+1,length(x(:,k)));
for j = 1:2*length(x_a)+1
    % Calculate transformed sampling points
    [~, X_k(j,:)] = rk4(f, Xx(j,:), u(:,k), Xv(j,:), [0 dt]);  
end

%% Step 4: Compute the predicted mean, covariance and cross-covariance
% Calculate transformed mean
x_k = Wm*X_k;

% Calculate transformed deviations
X1 = X_k-x_k(ones(1,2*length(x_a)+1),:);
X2 = Xx - (x(:,k)*ones(1,2*length(x_a)+1))';

% Calculate transformed covariance
P_k = X1'*diag(Wc)*X1;
    
% Calculate cross covariance
C_k = X2'*diag(Wc)*X1;

%% Step 5: Compute the smoother Gain, smoothed mean and covariance
% Calculate smoother gain
D_k = C_k/P_k;
  
% Calculate smoothed mean
xs(:,k) = x(:,k) + D_k*(xs(:,k+1) - x_k');
    
% Calculate smoothed covariance
Ps(:,:,k) = P(:,:,k) + D_k * (Ps(:,:,k+1) - P_k) * D_k';
    
Y_k = zeros(2*length(x_a)+1,length(R));
for j = 1:2*length(x_a)+1
    % Calculate transformed sampling points
    Y_k(j,:) = h(0,X_k(j,:), u, Xn(j,:)); 
end

% Calculate transformed mean
y_k = Wm*Y_k;
    
ys(:,k) = y_k;
    
% Calculate transformed deviations
Y1 = Y_k-y_k(ones(1,2*length(x_a)+1),:);
    
% Calculate the transformed covariance of the output
Pyys(:,:,k)  = Y1'*diag(Wc)*Y1 + R; 
    
    % Print the current data point that is being filtered
    if mod(k,100) == 0
        text = sprintf('Running RTS: %d/%d', k, length(x));

        cpb.setValue(k);
        cpb.setText(text);
    end 
end

% Stop console progress bar
cpb.stop();

% Add an extra value for the covariance matrix of the measurement to obtain
% the correct length
Pyys(:,:,length(Pyys)+1) = Pyys(:,:,length(Pyys));

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