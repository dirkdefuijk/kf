function [DataStruct] = FillFPRdata(x, Pxx, y, Pyy, u, t, DataStruct)
% FillFPRdata.m is a script which stores the flight path reconstructed
% variables into the main data structure. It stores the input (corrected
% with the bias term), the state and the measurement information.
%
% Inputs:   
%           x:          reconstructed state variables, consisting of:
%                       x = [x, y, z, u, v, w, phi, theta, psi, Wxe, Wye, Caup,
%                            lambda_x, lambda_y, lambda_z, lambda_p, lambda_q, 
%                            lambda_r, alpha_v]^T
%           Pxx:        Covariance matrix corresponding to the reconstructed state variables
%           y:          the filtered measurements, consisting of:
%                       y = [x, y, h, xdot, ydot, hdot, phi, theta, psi, 
%                       Vtas, alpha_v, pseudo_beta]^T     
%           Pyy:        Covariance matrix corresponding to the filtered measurements
%           t:          time vector
%           DataStruct: data structure containing all important information
%
% Outputs:  
%           DataStruct: re-defined data structure with FPR variables
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016 - Version 1.0
%% Reconstruct the input
DataStruct.Input.u.ax.filtered = u(1,:) - x(10,:);
DataStruct.Input.u.ay.filtered = u(2,:) - x(11,:);
DataStruct.Input.u.az.filtered = u(3,:) - x(12,:);
DataStruct.Input.u.p.filtered = u(4,:) - x(13,:);
DataStruct.Input.u.q.filtered = u(5,:) - x(14,:);
DataStruct.Input.u.r.filtered = u(6,:) - x(15,:);

%% Reconstruct the state
DataStruct.State.x.time = t;
DataStruct.State.x.data = x(1,:);
DataStruct.State.x.std = sqrt(squeeze(Pxx(1,1,:))');

DataStruct.State.y.time = t;
DataStruct.State.y.data = x(2,:);
DataStruct.State.y.std = sqrt(squeeze(Pxx(2,2,:))');

DataStruct.State.z.time = t;
DataStruct.State.z.data = x(3,:);
DataStruct.State.z.std = sqrt(squeeze(Pxx(3,3,:))');

DataStruct.State.u.time = t;
DataStruct.State.u.data = x(4,:);
DataStruct.State.u.std = sqrt(squeeze(Pxx(4,4,:))');

DataStruct.State.v.time = t;
DataStruct.State.v.data = x(5,:);
DataStruct.State.v.std = sqrt(squeeze(Pxx(5,5,:))');

DataStruct.State.w.time = t;
DataStruct.State.w.data = x(6,:);
DataStruct.State.w.std = sqrt(squeeze(Pxx(6,6,:))');

DataStruct.State.phi.time = t;
DataStruct.State.phi.data = x(7,:);
DataStruct.State.phi.std = sqrt(squeeze(Pxx(7,7,:))');

DataStruct.State.theta.time = t;
DataStruct.State.theta.data = x(8,:);
DataStruct.State.theta.std = sqrt(squeeze(Pxx(8,8,:))');

DataStruct.State.psi.time = t;
DataStruct.State.psi.data = x(9,:);
DataStruct.State.psi.std = sqrt(squeeze(Pxx(9,9,:))');

DataStruct.State.lambda_Ax.time = t;
DataStruct.State.lambda_Ax.data = x(10,:);
DataStruct.State.lambda_Ax.std = sqrt(squeeze(Pxx(10,10,:))');

DataStruct.State.lambda_Ay.time = t;
DataStruct.State.lambda_Ay.data = x(11,:);
DataStruct.State.lambda_Ay.std = sqrt(squeeze(Pxx(11,11,:))');

DataStruct.State.lambda_Az.time = t;
DataStruct.State.lambda_Az.data = x(12,:);
DataStruct.State.lambda_Az.std = sqrt(squeeze(Pxx(12,12,:))');

DataStruct.State.lambda_p.time = t;
DataStruct.State.lambda_p.data = x(13,:);
DataStruct.State.lambda_p.std = sqrt(squeeze(Pxx(13,13,:))');

DataStruct.State.lambda_q.time = t;
DataStruct.State.lambda_q.data = x(14,:);
DataStruct.State.lambda_q.std = sqrt(squeeze(Pxx(14,14,:))');

DataStruct.State.lambda_r.time = t;
DataStruct.State.lambda_r.data = x(15,:);
DataStruct.State.lambda_r.std = sqrt(squeeze(Pxx(15,15,:))');

DataStruct.State.Wxe.time = t;
DataStruct.State.Wxe.data = x(16,:);
DataStruct.State.Wxe.std = sqrt(squeeze(Pxx(16,16,:))');

DataStruct.State.Wye.time = t;
DataStruct.State.Wye.data = x(17,:);
DataStruct.State.Wye.std = sqrt(squeeze(Pxx(17,17,:))');

DataStruct.State.Caup.time = t;
DataStruct.State.Caup.data = x(18,:);
DataStruct.State.Caup.std = sqrt(squeeze(Pxx(18,18,:))');



%
DataStruct.State.u_L.time = t;
DataStruct.State.u_L.data = x(19,:);
DataStruct.State.u_L.std = sqrt(squeeze(Pxx(19,19,:))');

DataStruct.State.u_R.time = t;
DataStruct.State.u_R.data = x(20,:);
DataStruct.State.u_R.std = sqrt(squeeze(Pxx(20,20,:))');

DataStruct.State.v_L.time = t;
DataStruct.State.v_L.data = x(21,:);
DataStruct.State.v_L.std = sqrt(squeeze(Pxx(21,21,:))');

DataStruct.State.v_R.time = t;
DataStruct.State.v_R.data = x(22,:);
DataStruct.State.v_R.std = sqrt(squeeze(Pxx(22,22,:))');

DataStruct.State.w_L.time = t;
DataStruct.State.w_L.data = x(23,:);
DataStruct.State.w_L.std = sqrt(squeeze(Pxx(23,23,:))');

DataStruct.State.w_R.time = t;
DataStruct.State.w_R.data = x(24,:);
DataStruct.State.w_R.std = sqrt(squeeze(Pxx(24,24,:))');


%% Reconstruct the measurements
u = DataStruct.State.u.data;
v = DataStruct.State.v.data;
w = DataStruct.State.w.data;
phi = DataStruct.State.phi.data;
theta = DataStruct.State.theta.data;
psi = DataStruct.State.psi.data;
Wxe = DataStruct.State.Wxe.data;
Wye = DataStruct.State.Wye.data;

%
u_L = DataStruct.State.u_L.data;
u_R = DataStruct.State.u_R.data;
v_L = DataStruct.State.v_L.data;
v_R = DataStruct.State.v_R.data;
w_L = DataStruct.State.w_L.data;
w_R = DataStruct.State.w_R.data;

DataStruct.Measurement.x.time = t;
DataStruct.Measurement.x.reconstructed = x(1,:);
DataStruct.Measurement.x.innovation = DataStruct.Measurement.x.raw - y(1,:);
DataStruct.Measurement.x.std = sqrt(squeeze(Pyy(1,1,:))');

DataStruct.Measurement.y.time = t;
DataStruct.Measurement.y.reconstructed = x(2,:);
DataStruct.Measurement.y.innovation = DataStruct.Measurement.y.raw - y(2,:);
DataStruct.Measurement.y.std = sqrt(squeeze(Pyy(2,2,:))');

DataStruct.Measurement.h.time = t;
DataStruct.Measurement.h.reconstructed = -x(3,:);
DataStruct.Measurement.h.innovation = DataStruct.Measurement.h.raw - y(3,:);
DataStruct.Measurement.h.std = sqrt(squeeze(Pyy(3,3,:))');

DataStruct.Measurement.xdot.time = t;
DataStruct.Measurement.xdot.reconstructed = (u.*cos(theta) + (v.*sin(phi) + ...
    w.*cos(phi)).*sin(theta)).*cos(psi) - (v.*cos(phi) - w.*sin(phi)).*sin(psi) + Wxe;
DataStruct.Measurement.xdot.innovation = DataStruct.Measurement.xdot.raw - y(4,:);
DataStruct.Measurement.xdot.std = sqrt(squeeze(Pyy(4,4,:))');

DataStruct.Measurement.ydot.time = t;
DataStruct.Measurement.ydot.reconstructed = (u.*cos(theta) + (v.*sin(phi) + ...
    w.*cos(phi)).*sin(theta)).*sin(psi) + (v.*cos(phi) - w.*sin(phi)).*cos(psi) + Wye;
DataStruct.Measurement.ydot.innovation = DataStruct.Measurement.ydot.raw - y(5,:);
DataStruct.Measurement.ydot.std = sqrt(squeeze(Pyy(5,5,:))');

DataStruct.Measurement.hdot.time = t;
DataStruct.Measurement.hdot.reconstructed = u.*sin(theta) - (v.*sin(phi) + ...
                                            w.*cos(phi)).*cos(theta);
DataStruct.Measurement.hdot.innovation = DataStruct.Measurement.hdot.raw - y(6,:);
DataStruct.Measurement.hdot.std = sqrt(squeeze(Pyy(6,6,:))');

DataStruct.Measurement.phi.time = t;
DataStruct.Measurement.phi.reconstructed = x(7,:);
DataStruct.Measurement.phi.innovation = DataStruct.Measurement.phi.raw - y(7,:);
DataStruct.Measurement.phi.std = sqrt(squeeze(Pyy(7,7,:))');

DataStruct.Measurement.theta.time = t;
DataStruct.Measurement.theta.reconstructed = x(8,:);
DataStruct.Measurement.theta.innovation = DataStruct.Measurement.theta.raw - y(8,:);
DataStruct.Measurement.theta.std = sqrt(squeeze(Pyy(8,8,:))');

DataStruct.Measurement.psi.time = t;
DataStruct.Measurement.psi.reconstructed = x(9,:);
DataStruct.Measurement.psi.innovation = DataStruct.Measurement.psi.raw - y(9,:);
DataStruct.Measurement.psi.std = sqrt(squeeze(Pyy(9,9,:))'); 

DataStruct.Measurement.Vtas.time = t;
DataStruct.Measurement.Vtas.reconstructed = sqrt(u.^2 + v.^2 + w.^2);
DataStruct.Measurement.Vtas.innovation = DataStruct.Measurement.Vtas.raw - y(10,:);
DataStruct.Measurement.Vtas.std = sqrt(squeeze(Pyy(10,10,:))'); 

DataStruct.Measurement.alpha.time = t;
DataStruct.Measurement.alpha.reconstructed = atan(w./u);
DataStruct.Measurement.alpha.innovation = DataStruct.Measurement.alpha.raw - y(11,:);
DataStruct.Measurement.alpha.std = sqrt(squeeze(Pyy(11,11,:))'); 

DataStruct.Measurement.beta.time = t;
DataStruct.Measurement.beta.reconstructed = atan(v./(sqrt(u.^2 + w.^2)));
DataStruct.Measurement.beta.innovation = DataStruct.Measurement.beta.raw - y(12,:);
DataStruct.Measurement.beta.std = sqrt(squeeze(Pyy(12,12,:))'); 

DataStruct.Measurement.alpha_L.time = t;
DataStruct.Measurement.alpha_L.reconstructed = atan(w_L./u_L);
DataStruct.Measurement.alpha_L.innovation = DataStruct.Measurement.alpha_L.raw - y(13,:);
DataStruct.Measurement.alpha_L.std = sqrt(squeeze(Pyy(13,13,:))'); 

DataStruct.Measurement.alpha_R.time = t;
DataStruct.Measurement.alpha_R.reconstructed = atan(w_R./u_R);
DataStruct.Measurement.alpha_R.innovation = DataStruct.Measurement.alpha_R.raw - y(14,:);
DataStruct.Measurement.alpha_R.std = sqrt(squeeze(Pyy(14,14,:))'); 
end