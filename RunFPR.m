% Main script for running the flight path reconstruction (FPR) on data
% gathered by the Cessna Citation II including a boom signal.
%
% Original by:  Laurens van Horssen & Marlon van den Hoek
% Update by:    Joost van Ingen
% Last update:  30 March 2017

close all; clear; clc
addpath([pwd '/func']) 

% Select the files which need to be filtered
% OPTION 1: Choose the files manually using the uigetfile window, with this
%           all selected files must be located in the same folder.
[FileName,PathName] = uigetfile('*.mat','Load Data File', ...
                             '../../data/processed','MultiSelect','on');
% OPTION 2: Use the make_list_of_data_paths function to automatically
%           select all data sets (can be from multiple folders).
% [FileName,PathName] = make_list_of_data_paths();

% Save figures regarding innovation, measurements, state and inputs. 
% Options: either 'on' or 'off'
savefig = 'on';

% Type of filter, the following options are available: 'EKF', 'IEKF',
% 'UKF', 'ERTSS', 'IERTSS', 'URTSS'
type_kf = 'UKF';

% Choose whether to pre-filter the inputs (ax, ay, az, p, q and r) with a
% lowpass filter, to reduce high frequency noise. Options: either
% 'on' or 'off'
filter_inputs = 'on';

% Choose the cut-off frequencies of the different pre-filter procedures
filter_fc_ahrs    = 1.5; % [Hz] cutoff freq for ahrs signals
filter_fc_synchro = 4.0; % [Hz] cufoff freq for control surf. + boom signals

% Choose the filter order
filter_fo = 4; % [-]

% Choose the filter type: 'Butterworth'
filter_type = 'Butterworth';

% Determines the spread of the sigma points around the mean state, a larger
% value of alpha means a larger spread. Value must be between 0 < alpha < 1
ukf_alpha = 0.3; 

% Secondy scaling parameter, usually set to zero.
ukf_kappa = 0;

% Scaling parameter to incorporate prior knowledge of the distribution of
% the state. For Gaussian distribution beta = 2 is optimal.
ukf_beta = 2;

% Maximum error before stopping iterations of the IEKF
iekf_epsilon = 1e-10;

% Set maximum number of iterations in case the IEKF is chosen
iekf_maxIterations = 500;

% Process noise covariance matrix   "default":
Q = diag([ 0.02, ...        % ax    0.01
           0.02, ...        % ay    0.03
           0.02, ...        % az    0.06
           0.001, ...       % p     0.001
           0.001, ...       % q     0.001
           0.001]).^2;      % r     0.004
       
% Measurement noise covariance matrix           "default":
R = diag([ 1.0, ...                  % x        1.0
           1.0, ...                  % y        1.0
           1.0, ...                  % h        0.3
           0.1, ...                  % xdot     0.1
           0.1, ...                  % ydot     0.1
           0.1, ...                  % hdot     0.1
           0.1*pi/180, ...           % phi      0.1*pi/180
           0.1*pi/180, ...           % theta    0.1*pi/180
           0.15*pi/180, ...          % psi      0.15*pi/180
           0.1, ...                  % Vtas     0.1
           0.1*pi/180, ...           % alpha    0.1*pi/180
           0.1*pi/180, ...           % beta     0.1*pi/180
           0.1*pi/180, ...           % alpha_L  0.1*pi/180
           0.1*pi/180, ...           % alpha_R  0.1*pi/180
           ]).^2;

% Add folders to run the Kalman filter
addpath([pwd '/system'])
addpath([pwd '/func/Engine']) 
addpath([pwd '/filters']) 
addpath([pwd '/func/ConsoleProgressBar']) 
addpath([pwd '/func/Fill']) 
addpath([pwd '/func/Plotting']) 
addpath([pwd '/define/']) 
addpath([pwd '/../output']) 

% Add data folder
if iscell(PathName)
  for i = 1 : length(PathName)
      addpath(PathName{i})
  end
else
  addpath(PathName)
end

%% Determine how many runs are needed for the Kalman Filter
if iscell(FileName)
    M = length(FileName);
else
    M = 1;   
end

%% Run the Kalman Filter
for i = 1:M
    % Load data set
    if iscell(FileName)
        load(FileName{i})
        disp(['Dataset (',num2str(i),'/',num2str(M),') is being reconstruced: ', FileName{i}]);
    else
        load(FileName)   
        disp(['Dataset is being reconstruced: ', FileName]);
    end
       
    % Create time vector. Since data must be resampled it does not really
    % matter which time signal is used.
    t = flightdata.analog.vane_AOA.time;

    % Obtain sample time
    dt = mean(diff(t));

    % Initialize final data structure
    DataStruct = InitDataStructure(type_kf);

    % Load the system equations in the function handle f
    f = @NavigationSystem;

    % Load the measurements equations in the function handle h
    h = @ObservationModel;
    
    % Load the noise input matrix (only needed for the IEKF)
    G = @NoiseInputMatrix;
      
    % Load the inputs for the system equations
    u = DefineInput(flightdata);

    % Pre-filter the input in case it is chosen do to so
    if strcmp(filter_inputs,'on')
        u = FilterInputs(u, dt, filter_fc_ahrs, filter_fo, filter_type);
    elseif strcmp(filter_inputs,'off')
        % Do nothing
    else
        error('Unknown input for filter, use either "on" or "off", line 14 in RunFPR.m!')
    end
    
    % Correct acceleration measurements for c.g. location if the new
    % dataset is used, since the measurement is not in the c.g. The
    % accelerations in the old data set are measured close to the c.g. and
    % it is thus not necessary to compensate for the c.g. location
    u = Correct_cg_loc(u, massmodel, dt, filter_inputs, ...
                       filter_fc_ahrs, filter_fo, filter_type);

    % Load the measurements
    [meas,DataStruct] = DefineMeasurements(flightdata,DataStruct);

    % Filter the boom resonant oscillations
    if strcmp(filter_inputs,'on')
        [meas] = FilterBoomSignals(meas, dt, filter_fc_synchro, ...
                                    filter_fo, filter_type);
    elseif strcmp(filter_inputs,'off')
        % Do nothing
    else
        error('Unknown input for filter, use either "on" or "off", line 14 in RunFPR.m!')
    end

    
    % Initial state estimate
    x0 = DefineInitialCondition(meas);

    % Initial covariance matrix
    P0 = DefineInitialCovariance();

    % Store data struct information
    FillDataStructInfo

    % Fill raw data structure of the input u
    FillDataStructRawData

    % Add engine data
    FillDataStructEngine

    for l = 1:size(x0,2)
        %% Start Kalman Filter

        % Number of data points to be filtered
        N = length(meas);

        % allocate memory
        xV = zeros(length(x0),N);                % State estimate
        yV = zeros(length(R),N);                 % Measurements estimate
        Pyy = zeros(length(R),length(R),N);      % Measurement covariance matrix
        Pxx = zeros(length(x0),length(x0),N);    % State covariance matrix

        % Set the initial conditions
        x = x0(:,l);
        P = P0;

        % Offset distances between c.g. and boom vanes, declared as
        % global variables so the Kalman filter funcs can access them
        global xva zva xvb zvb xva_L xva_R yva_L yva_R
        
        % Boom vane positions w.r.t. the datum line. Datum line is 
        % 0.38 meter in front of, and 2.31 meter below the nose.
        inch2m = 0.0254; 
        xpos_va = -27.5*inch2m; % [m] !! VEHICLE FRAME* !! 
        zpos_va = 128.5*inch2m; % [m] !! VEHICLE FRAME* !!
        xpos_vb = -16.3*inch2m; % [m] !! VEHICLE FRAME* !!
        zpos_vb = 128.5*inch2m; % [m] !! VEHICLE FRAME* !!
        
        % Open console progress bar and set options
        SetCPBoptions

        % Start console progress bar
        cpb.start();

        for m = 1:N
            
            % Center of gravity location w.r.t the datum line. Due to fuel
            % burn, c.g. can move during flight.
            xpos_cg = massmodel.x_cg.data(m); % [m] !! VEHICLE FRAME* !! 
            zpos_cg = massmodel.z_cg.data(m); % [m] !! VEHICLE FRAME* !! 
            
            % *note:    left-handed reference frame, with pos directions:
            %           X -> nose to tail, Y -> left wing, Z -> up
            xva = xpos_cg - xpos_va; % [m] !! BODY FRAME !! 
            zva = zpos_cg - zpos_va; % [m] !! BODY FRAME !!
            xvb = xpos_cg - xpos_vb; % [m] !! BODY FRAME !! 
            zvb = zpos_cg - zpos_vb; % [m] !! BODY FRAME !! 

            xva_L = xva;
            xva_R = xva;
            yva_L = 15./2;
            yva_R = -15./2;
            
            % Select the type of Kalman filter 
            switch type_kf
                case 'EKF'
                    [x, P, y_k, Pyy_t] = EKF(f, h, G, x, P, u(:,m), meas(:,m), Q, R, dt);
                case 'IEKF'
                    [x, P, y_k, Pyy_t] = IEKF(f, h, G, x, P, u(:,m), meas(:,m), Q, R, dt, ...
                        'epsilon', iekf_epsilon, 'maxIterations', iekf_maxIterations);
                case 'UKF'
                    [x, P, y_k, Pyy_t] = UKF(f, h, x, P, u(:,m), meas(:,m), Q, R, dt, ...
                        'alpha',ukf_alpha,'kappa',ukf_kappa,'beta',ukf_beta);
                case 'ERTSS'
                    [x, P, y_k, Pyy_t] = EKF(f, h, G, x, P, u(:,m), meas(:,m), Q, R, dt);
                case 'IERTSS'
                    [x, P, y_k, Pyy_t] = IEKF(f, h, G, x, P, u(:,m), meas(:,m), Q, R, dt, ...
                        'epsilon', iekf_epsilon, 'maxIterations', iekf_maxIterations);
                case 'URTSS'
                    [x, P, y_k, Pyy_t] = UKF(f, h, x, P, u(:,m), meas(:,m), Q, R, dt, ...
                        'alpha',ukf_alpha,'kappa',ukf_kappa,'beta',ukf_beta);
                otherwise
                    error('Type of Kalman filter is unknown!')
            end

            % Store the filtered data
            xV(:,m) = x;
            yV(:,m) = y_k;
            Pyy(:,:,m) = Pyy_t;
            Pxx(:,:,m) = P;

            % Print the current data point that is being filtered
            if mod(m,100) == 0
                text = sprintf('Progress: %d/%d', m, N);

                cpb.setValue(m);
                cpb.setText(text);
            end   
        end

        % Stop console progress bar
        cpb.stop();
        
        % In case it is needed, apply the Rauch-Tung-Striebel smoother
        if strcmp(type_kf,'URTSS')
            [xV, Pxx, yV, Pyy] = URTSS(f, h, xV, Pxx, u, Q, R, dt, ...
                        'alpha',ukf_alpha,'kappa',ukf_kappa,'beta',ukf_beta);
        elseif strcmp(type_kf,'IERTSS') || strcmp(type_kf,'ERTSS')
            [xV, Pxx, yV, Pyy] = RTS(f, h, G, xV, Pxx, u, Q, R, dt);
        end
        

        %% Save data
        FilteredData{l} = FillFPRdata(xV, Pxx, yV, Pyy, u, t, DataStruct);

    end

    % Calculate the rotational accelerations and put them in the data structure
    [FilteredData{1}] = CalcRotAcc(FilteredData{1});

    % Add gravitational component of the earth to the accelerations
    [FilteredData{1}] = AddGravityAcc(FilteredData{1});

    % Add engine thrust
    [FilteredData{1}] = AddEngineThrust(FilteredData{1});

    % Filter control deflections
    [FilteredData{1}] = FilterControlDeflection(FilteredData{1}, ...
                        filter_fc_synchro, filter_fo, dt, filter_type);

    % Store the datastructure into FPR
    FPR = FilteredData{1};

    %% Save the filtered file  
    % If the operating system is unix/ mac then use forward slash
    if isunix || ismac
        SaveDir = [pwd,'/../output/'];
    else
        SaveDir = [pwd,'\..\output\'];
    end
    
    % First check existence of directory, if not existing create new
    if exist(SaveDir,'dir') == 0
        mkdir(SaveDir);
    end
    
    % Make a folder in which the figurs can be saved
    % If the operating system is unix/ mac then use forward slash
    if isunix || ismac
        if iscell(FileName)
            FigureFolder = [SaveDir,'Figures/',FileName{i}(1:end-4)];
        else
            FigureFolder = [SaveDir,'Figures/',FileName(1:end-4)];  
        end
    else
        if iscell(FileName)
            FigureFolder = [SaveDir,'Figures\',FileName{i}(1:end-4)];
        else
            FigureFolder = [SaveDir,'Figures\',FileName(1:end-4)];  
        end
    end
   
    % First check existence of directory, if not existing create new
    if exist(FigureFolder,'dir') == 0
        mkdir(FigureFolder);
    end
    
    % Save filtered flightdata
    if iscell(FileName)
        matfile = fullfile(SaveDir, ['FPRCV_',FileName{i}]);
    else
        matfile = fullfile(SaveDir, ['FPRCV_',FileName]);
    end
    
    save(matfile,'FPR','flightdata', 'massmodel');

    %% Show flight path reconstruction results
    if strcmp(savefig,'on') == 1
        close all
        ShowFPRresults(FilteredData)

        movefile('Innovation.eps',FigureFolder)
        movefile('Input.eps',FigureFolder)
        movefile('Measurements.eps',FigureFolder)
        movefile('State1.eps',FigureFolder)
        movefile('State2.eps',FigureFolder)
%         movefile('Convergence.eps',FigureFolder)
    end
    
%     clearvars FilteredData flightdata massmodel FPR
    
    [~,name,ext] = fileparts(matfile);
    disp(['Data is saved to: ' SaveDir]) 
    disp([' under the file name: ', name, ext])
end

disp('Kalman Filter is finished!')