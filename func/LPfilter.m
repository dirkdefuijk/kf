function [y] = LPfilter(x, fc, fo, dt, type)
% LPfilter.m is a script which low pass filters the signal. For now the
% only type that is available is the Butterworth filter
%
% Inputs:   
%           x:              input signal 
%           fc:             Cut-off frequency   [Hz]
%           fo:             Filter order        [-]      
%           dt:             sampling rate       [s]
%           type:           Type of the filter that is being used, e.g. Butterworth   
%
% Outputs:  
%           y:              filtered signal
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Low pass filter

% Calculate the sample and Nyquist frequencies
fs = 1/dt;                  % Sample frequency
fn = fs/2;                  % Nyquist frequency

% Select the correct type of filter
   switch type
       case 'Butterworth'
          % Obtain the coefficients for a Butterworth Filter
          [num,den] = butter(fo,fc/fn);      % IIR Filter coefficients
          
          % Low pass filter the input signal
          y = filtfilt(num,den,x);           
       otherwise
           error('Unknown filter type, currently only "Butterworth" is available')
   end   