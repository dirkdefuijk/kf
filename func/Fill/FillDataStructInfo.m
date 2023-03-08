% FillDataStructInfo.m is a script which is used fill the FPR datastructure
% with the general information on things such as: which type of filter is
% used, was the data pre-filtered, what are the initial conditions etc.
%
% Inputs:   
%           none
%
% Outputs:  
%           none (filled data structure with general information)
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%%
% Store the type of filter that is used
DataStruct.info.FilterType = type_kf;

% Store the initial condition vector
DataStruct.info.FilterParameters.x0.data = x0(:,1);

% Store the initial covariance matrix
DataStruct.info.FilterParameters.P0.data = P0;

% Store the process noise covariance matrix
DataStruct.info.FilterParameters.Q.data = Q;

% Store the measurement noise covariance matrix
DataStruct.info.FilterParameters.R.data = R;

if strcmp(type_kf,'UKF') == 1 || strcmp(type_kf,'URTSS') == 1 
    % Scaling parameter for the UKF. alpha determines the spread of the sigma
    % points around the mean state
    DataStruct.info.FilterParameters.alpha.data = ukf_alpha;

    % Secondy scaling parameter, usually set to zero.
    DataStruct.info.FilterParameters.kappa.data = ukf_kappa;

    % Scaling parameter to incorporate prior knowledge of the distribution of
    % the state. For Gaussian distribution beta = 2 is optimal.
    DataStruct.info.FilterParameters.beta.data = ukf_beta;
    
elseif  strcmp(type_kf,'IEKF') == 1 || strcmp(type_kf,'IERTSS') == 1
    % Epsilon determines that maximum error between two estimated states,
    % before the iteration is ended
    DataStruct.info.FilterParameters.epsilon.data = iekf_epsilon;

    % maxIterations determines the maximum number of iterations that were
    % allowed in case of an iterated Kalman filter
    DataStruct.info.FilterParameters.maxIterations.data = iekf_maxIterations;

elseif strcmp(type_kf,'EKF') == 1 || strcmp(type_kf,'ERTSS') == 1 
    % Do nothing, since there are no options to choose   
else
    error('Type of Kalman filter that is used is wrong. Options that could be chosen are: EKF, IEKF, UKF, ERTSS, IERTSS, URTSS')
end

% Store information on the pre-filter
if strcmp(filter_inputs,'on')
    DataStruct.info.Prefiltered.true = 'yes';
    DataStruct.info.Prefiltered.fc.ahrs = filter_fc_ahrs;
    DataStruct.info.Prefiltered.fc.synchro = filter_fc_synchro;
    DataStruct.info.Prefiltered.fc.Unit = 'Hz';
    DataStruct.info.Prefiltered.fo.data = filter_fo;
    DataStruct.info.Prefiltered.type.data = filter_type;
elseif strcmp(filter_inputs,'off')
    DataStruct.info.Prefiltered.true = 'no';
else
    error('Unknown input for variable "filter_inputs" ')
end

% Store data on the gear and flap settings
DataStruct.info.Configuration.gear = 'Unknown'; % [1 = Down, 0 = Up]
DataStruct.info.Configuration.flaps = 'Unknown'; % [degrees]
