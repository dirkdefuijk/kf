function [DataStruct] = CalcRotAcc(DataStruct)
% CalcRotAcc.m is a script which calculates the angular accelerations
% by numerically differentiating the angular rates.
%
% Inputs:   
%           DataStruct:     Data structure without angular accelerations
%
% Outputs:  
%           DataStruct:     Data structure containing the angular accelerations
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Calculate the angular accelerations, add a leading zero to get correct vector size
dt = diff(DataStruct.Input.u.p.time);
pdot = diff(DataStruct.Input.u.p.filtered)./dt; pdot = [0 pdot];

dt = diff(DataStruct.Input.u.q.time);
qdot = diff(DataStruct.Input.u.q.filtered)./dt; qdot = [0 qdot];

dt = diff(DataStruct.Input.u.r.time);
rdot = diff(DataStruct.Input.u.r.filtered)./dt; rdot = [0 rdot];

%% Store the angular accelerations in the data structure
DataStruct.Input.RotAcc.pdot.time = DataStruct.Input.u.p.time;
DataStruct.Input.RotAcc.pdot.data = pdot;

DataStruct.Input.RotAcc.qdot.time = DataStruct.Input.u.q.time;
DataStruct.Input.RotAcc.qdot.data = qdot;

DataStruct.Input.RotAcc.rdot.time = DataStruct.Input.u.r.time;
DataStruct.Input.RotAcc.rdot.data = rdot;    