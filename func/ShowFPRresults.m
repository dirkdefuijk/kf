function [] = ShowFPRresults(FilteredData)
% ShowFPRresults.m is a script which acts as a parent for the figure
% plotting
%
% Inputs:   
%           FilteredData:  Cell containing up to 3 data structures
%
% Outputs:  
%           none
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Check if more than three initial conditions were used
N = length(FilteredData);

if length(N) > 3
    error(['Not able to plot more than three signals per time. ', ...
      'Rerun the Kalman filter with a maximum of three initial ', ...
      'conditions or initial covariance matrices at a time!'])
end

%% Plot innovations
plotInnovation(FilteredData)

%% Plot measurements
plotResults(FilteredData)

%% Plot the state
plotState(FilteredData)

%% Plot the input
plotInput(FilteredData)

%% Plot the convergence
% length(FilteredData)
% if length(FilteredData) > 1
%     plotConvergence(FilteredData)
% end
end