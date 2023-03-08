function [DataStruct] = FilterControlDeflection(DataStruct, fc, fo, dt, type)
% FilterControlDeflection.m is a script which filters the control
% deflections given a certain cut-off frequency.
%
% Inputs:   
%           DataStruct: data structure without filtered control deflections
%           fc:             Cut-off frequency   [Hz]
%           dt:             sampling rate       [s]
%           fo:             Filter order        [-]      
%           type:    Type of the filter that is being used, e.g. Butterworth   
%
% Outputs:  
%           DataStruct: data structure containing the filtered control deflections
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0   
%% Filter the control deflections

% Obtain the sample and Nyquist frequency
fs = 1/dt;                  % Sample frequency
fn = fs/2;                  % Nyquist frequency

% Select the correct type of filter
   switch type
       case 'Butterworth'
          [num,den] = butter(fo,fc/fn);      % IIR Filter coefficients
          
          % Filterd the elevator deflection
          DataStruct.Input.delta.de.filtered = filtfilt(num,den,DataStruct.Input.delta.de.raw);
          
          % Filterd the aileron deflection
          DataStruct.Input.delta.da.filtered = filtfilt(num,den,DataStruct.Input.delta.da.raw);  
          
          % Filterd the rudder deflection
          DataStruct.Input.delta.dr.filtered = filtfilt(num,den,DataStruct.Input.delta.dr.raw); 
          
          % If the flaps angle is empty, then do nothing. Otherwise filter
          % the flap angle as well.
          if isempty(DataStruct.Input.delta.flaps.raw)
              DataStruct.Input.delta.flaps.filtered = [];
          else
              DataStruct.Input.delta.flaps.filtered = filtfilt(num,den,DataStruct.Input.delta.flaps.raw);  
          end
          
       otherwise
           error('Unknown filter type, currently only "Butterworth" is available')
   end