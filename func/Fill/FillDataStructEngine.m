% FillDataStructEngine.m is a script which is used fill the engine
% parameters in the structure. This includes all measured engine settings.
%
% Inputs:   
%           none
%
% Outputs:  
%           none (filled data structure with raw measurement data)
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Left engine parameters
% Obtain and store the fuel mass flow of the left engine
switch flightdata.analog.lh_engine_FMF.Unit
                case 'lbs/hr'
                    DataStruct.Engine.lh.FMF.time = flightdata.analog.lh_engine_FMF.time;
                    DataStruct.Engine.lh.FMF.data = flightdata.analog.lh_engine_FMF.data*0.000125998;
                    DataStruct.Engine.lh.FMF.Sensor = 'analog/lh_engine_FMF';
                case 'kg/s'
                    DataStruct.Engine.lh.FMF.time = flightdata.analog.lh_engine_FMF.time;
                    DataStruct.Engine.lh.FMF.data = flightdata.analog.lh_engine_FMF.data;
                    DataStruct.Engine.lh.FMF.Sensor = 'analog/lh_engine_FMF';
    otherwise
        error('Unit of the fuel mass flow is unknown')
end

% Obtain and store the oil pressure of the left engine
switch flightdata.analog.lh_engine_OP.Unit
                case 'N/m^2'
                    DataStruct.Engine.lh.OP.time = flightdata.analog.lh_engine_OP.time;
                    DataStruct.Engine.lh.OP.data = flightdata.analog.lh_engine_OP.data*0.000145038; 
                    DataStruct.Engine.lh.OP.Sensor = 'analog/lh_engine_OP';
                case 'pa'
                    DataStruct.Engine.lh.OP.time = flightdata.analog.lh_engine_OP.time;
                    DataStruct.Engine.lh.OP.data = flightdata.analog.lh_engine_OP.data*0.000145038;
                    DataStruct.Engine.lh.OP.Sensor = 'analog/lh_engine_OP';
                case 'psi'
                    DataStruct.Engine.lh.OP.time = flightdata.analog.lh_engine_OP.time;
                    DataStruct.Engine.lh.OP.data = flightdata.analog.lh_engine_OP.data;
                    DataStruct.Engine.lh.OP.Sensor = 'analog/lh_engine_OP';
    otherwise
        error('Unit of the oil pressure is unknown')
end

% Obtain and store the inter turbine temperature of the left engine
switch flightdata.analog.lh_engine_itt.Unit
                case 'K'
                    DataStruct.Engine.lh.ITT.time = flightdata.analog.lh_engine_itt.time;
                    DataStruct.Engine.lh.ITT.data = flightdata.analog.lh_engine_itt.data - 273.15;
                    DataStruct.Engine.lh.ITT.Sensor = 'analog/lh_engine_itt';
                case 'deg C'
                    DataStruct.Engine.lh.ITT.time = flightdata.analog.lh_engine_itt.time;
                    DataStruct.Engine.lh.ITT.data = flightdata.analog.lh_engine_itt.data;
                    DataStruct.Engine.lh.ITT.Sensor = 'analog/lh_engine_itt';
    otherwise
        error('Unit of the inter turbine temperature is unknown')
end

% Obtain and store the fan speed of the left engine
switch flightdata.interval.lh_engine_fan_N1.Unit
                case '%'
                    DataStruct.Engine.lh.N1.time = flightdata.interval.lh_engine_fan_N1.time;
                    DataStruct.Engine.lh.N1.data = flightdata.interval.lh_engine_fan_N1.data;
                    DataStruct.Engine.lh.N1.Sensor = 'analog/lh_engine_fan_N1';
    otherwise
        error('Unit of the engine fan speed (N1) is unknown')
end

% Obtain and store the turbine speed of the left engine
switch flightdata.interval.lh_engine_turbine_N2.Unit
                case '%'
                    DataStruct.Engine.lh.N2.time = flightdata.interval.lh_engine_turbine_N2.time;
                    DataStruct.Engine.lh.N2.data = flightdata.interval.lh_engine_turbine_N2.data;
                    DataStruct.Engine.lh.N2.Sensor = 'analog/lh_engine_turbine_N2';
    otherwise
        error('Unit of the turbine fan speed (N2) is unknown')
end

% Obtain and store the fuel used of the left engine
switch flightdata.analog.lh_engine_FU.Unit
                case 'lbs'
                    DataStruct.Engine.lh.FU.time = flightdata.analog.lh_engine_FU.time;
                    DataStruct.Engine.lh.FU.data = flightdata.analog.lh_engine_FU.data*0.453592;
                case 'kg'
                    DataStruct.Engine.lh.FU.time = flightdata.analog.lh_engine_FU.time;
                    DataStruct.Engine.lh.FU.data = flightdata.analog.lh_engine_FU.data;
    otherwise
        error('Unit of the fuel used is unknown')
end

%% Right engine parameters
% Obtain and store the fuel mass flow of the right engine
switch flightdata.analog.rh_engine_FMF.Unit
                case 'lbs/hr'
                    DataStruct.Engine.rh.FMF.time = flightdata.analog.rh_engine_FMF.time;
                    DataStruct.Engine.rh.FMF.data = flightdata.analog.rh_engine_FMF.data*0.000125998;
                    DataStruct.Engine.rh.FMF.Sensor = 'analog/rh_engine_FMF';
                case 'kg/s'
                    DataStruct.Engine.rh.FMF.time = flightdata.analog.rh_engine_FMF.time;
                    DataStruct.Engine.rh.FMF.data = flightdata.analog.rh_engine_FMF.data;
                    DataStruct.Engine.rh.FMF.Sensor = 'analog/rh_engine_FMF';
    otherwise
        error('Unit of the fuel mass flow is unknown')
end

% Obtain and store the oil pressure of the right engine
switch flightdata.analog.rh_engine_OP.Unit
                case 'N/m^2'
                    DataStruct.Engine.rh.OP.time = flightdata.analog.rh_engine_OP.time;
                    DataStruct.Engine.rh.OP.data = flightdata.analog.rh_engine_OP.data*0.000145038;  
                    DataStruct.Engine.rh.OP.Sensor = 'analog/rh_engine_OP';
                    
                case 'pa'
                    DataStruct.Engine.rh.OP.time = flightdata.analog.rh_engine_OP.time;
                    DataStruct.Engine.rh.OP.data = flightdata.analog.rh_engine_OP.data*0.000145038;
                    DataStruct.Engine.rh.OP.Sensor = 'analog/rh_engine_OP';
                case 'psi'
                    DataStruct.Engine.rh.OP.time = flightdata.analog.rh_engine_OP.time;
                    DataStruct.Engine.rh.OP.data = flightdata.analog.rh_engine_OP.data;
                    DataStruct.Engine.rh.OP.Sensor = 'analog/rh_engine_OP';
    otherwise
        error('Unit of the oil pressure is unknown')
end

% Obtain and store the inter turbine temperature of the right engine
switch flightdata.analog.rh_engine_itt.Unit
                case 'K'
                    DataStruct.Engine.rh.ITT.time = flightdata.analog.rh_engine_itt.time;
                    DataStruct.Engine.rh.ITT.data = flightdata.analog.rh_engine_itt.data -273.15;
                    DataStruct.Engine.rh.ITT.Sensor = 'analog/rh_engine_itt';
                case 'deg C'
                    DataStruct.Engine.rh.ITT.time = flightdata.analog.rh_engine_itt.time;
                    DataStruct.Engine.rh.ITT.data = flightdata.analog.rh_engine_itt.data;
                    DataStruct.Engine.rh.ITT.Sensor = 'analog/rh_engine_itt';
    otherwise
        error('Unit of the inter turbine temperature is unknown')
end

% Obtain and store the fan speed of the right engine
switch flightdata.interval.rh_engine_fan_N1.Unit
                case '%'
                    DataStruct.Engine.rh.N1.time = flightdata.interval.rh_engine_fan_N1.time;
                    DataStruct.Engine.rh.N1.data = flightdata.interval.rh_engine_fan_N1.data;
                    DataStruct.Engine.rh.N1.Sensor = 'analog/rh_engine_fan_N1';
    otherwise
        error('Unit of the engine fan speed (N1) is unknown')
end

% Obtain and store the turbine speed of the right engine
switch flightdata.interval.rh_engine_turbine_N2.Unit
                case '%'
                    DataStruct.Engine.rh.N2.time = flightdata.interval.rh_engine_turbine_N2.time;
                    DataStruct.Engine.rh.N2.data = flightdata.interval.rh_engine_turbine_N2.data;
                    DataStruct.Engine.rh.N2.Sensor = 'analog/rh_engine_turbine_N2';
    otherwise
        error('Unit of the turbine fan speed (N2) is unknown')
end

% Obtain and store the fuel used of the right engine
switch flightdata.analog.rh_engine_FU.Unit
                case 'lbs'
                    DataStruct.Engine.rh.FU.time = flightdata.analog.rh_engine_FU.time;
                    DataStruct.Engine.rh.FU.data = flightdata.analog.rh_engine_FU.data*0.453592;
                case 'kg'
                    DataStruct.Engine.rh.FU.time = flightdata.analog.rh_engine_FU.time;
                    DataStruct.Engine.rh.FU.data = flightdata.analog.rh_engine_FU.data;
    otherwise
        error('Unit of the fuel used is unknown')
end

%% Set time vector the engine thrust
DataStruct.Engine.Thrust.time = flightdata.analog.rh_engine_FU.time;