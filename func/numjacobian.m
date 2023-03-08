function J = numjacobian(f,x)
% Numerical approximation of the Jacobian of a non-linear system
% in the form f(x(t),u(t),t) + G w(t), using five-point method
% or complex step arithmetic for improved accuracy. Note that this function 
% only works for function handles with four variables: t, x, u and w. Where
% t is the time instant (used for one step-ahead prediction), x the state
% vector, u the input vector and w the noise input vector.
%
% Input
%           f = function handle with, respectively, inputs t, x, u and w
%               representing the: time instant, state vector, input vector
%               and noise input vector
%           x = state vector
%
% REFERENCES:
% Lai, Kok-Lam, et al. 
% "New complex-step derivative approximations with application to second-order kalman filtering." 
% AIAA Guidance, Navigation and Control Conference, San Francisco, California. 2005.
%
% Made by: M.A. van den Hoek and L.J. van Horssen, September 2016
%% Numerical Calculation of the Jacobian
% Ensure correct input size, otherwise transpose the state vector
    if size(x,1) > size(x,2)
        x = x';
    end

    % set dimensions and step size
    u = zeros(1,20);                 % set input to zero, quick fix works up to 20 inputs
    n = size(x,2);                   % set size of state vector
    m = size(f(0,x,u,zeros(1,n)),2); % number of functions in system
    step = 1e-8;                     % difference step size

    % preallocate matrix for speed
    J = zeros(m,n);
    
    % Switch (1) complex arithmetic or (0) 5-point method
    % Note, if using atan2 in ackinematics its not possible to use complex
    % step method to calculate the Jacobian matrix
    complex = 0;
    
    if complex
        % Alternatively, use complex step arithmetic for higher accuracy
        h = n * eps;                            % differentiation step size, eps = machine epsilon
        for k = 1:n                             % loop for each independent variable 
            x1 = x;                             % reference point
            x1(k) = x1(k) + h*1i;               % increment in kth independent variable
            J(:,k) = imag(f(x1,u))/h;           % complex step differentiation
        end  
    else 
        
        % Create a check to see if a random walk is used in the navigation
        % system. This is done by evaluating the function two times. When
        % the results are not equal it is an indication that the
        % differential equation is contains a random walk. In that case the
        % calculated is invalid for these equations and will be set to zero
        evaluation1 = feval(f,0,x,u,zeros(1,n));
        evaluation2 = feval(f,0,x,u,zeros(1,n));
        
        % Check if the values are the same or not
        test = evaluation1 == evaluation2;
        
        % Find the indices at which the functions uses a random walk
        index = test == 0;
       
        for i=1:n
           % Make step vector for i-th index by only setting variable xi
           % non-zero to find df(x)/dxi
           h = zeros(1,n); h(i) = step;

           % Five-point method, higher order first derivative approximation
           xstep2 = x + 2*h;
           xstep1 = x + h;
           xstepm2 = x - 2*h;
           xstepm1 = x - h;

           J(:,i) = (-feval(f,0,xstep2,u,zeros(1,n))  + 8*feval(f,0,xstep1,u,zeros(1,n)) - ...
                    8*feval(f,0,xstepm1,u,zeros(1,n)) + feval(f,0,xstepm2,u,zeros(1,n))) ...
                    /(12*step);
        end
        
        % Set the value of the Jacobian to zero in case the differential
        % equations uses a random number/ random walk.
        J(index,:) = 0;
    end
end