% Filename :  ObservationModel.m
%
% Constructs the observation model and calculates the measurements using
% the system dynamic equations
%
% Inputs: 
%     xin:     state vector x:
%               x = [x, y, z, u, v, w,  
%                    phi, theta, psi, Wxe, Wye, Caup,
%                    lambda_x, lambda_y, lambda_z, lambda_p, lambda_q, 
%                    lambda_r, alpha_v]^T
%
%     uin:     input vector u:
%               u = [Ax, Ay, Az, p, q, r]^T
%
%     vin:      noise vector v:
%               v = [vx, vy, vh, vxdot, vydot, vhdot, vphi, vtheta, vpsi,
%                    vVtas, valpha, vbeta]^T 
%     
% Outputs:
%     z:        output vector z:
%               z = [x, y, h, xdot, ydot, hdot, phi, theta, psi, 
%                    Vtas, alpha_v, pseudo_beta]^T                
%
% Made by: M.A. van den Hoek and L.J. van Horssen, September 2016
% Update by: J.B. van Ingen, November 2016
function z = ObservationModel(t,xin,uin,vin)
% Extract the state and inputs and put them into variables
x = xin(1);     y = xin(2);     zin = xin(3);
u = xin(4);     v = xin(5);     w = xin(6);
phi = xin(7);   theta = xin(8); psi = xin(9);
Lx = xin(10);   Ly = xin(11);   Lz = xin(12);
Lp = xin(13);   Lq = xin(14);   Lr = xin(15); 
Wxe = xin(16);  Wye = xin(17);
Caup = xin(18); 

u_L = xin(19); u_R = xin(20);
v_L = xin(21); v_R = xin(22);
w_L = xin(23); w_R = xin(24);

Ax = uin(1); Ay = uin(2); Az = uin(3);
p = uin(4); q = uin(5); r = uin(6);

vx = vin(1);     vy = vin(2);       vh = vin(3);
vxdot = vin(4);  vydot = vin(5);    vhdot = vin(6); 
vphi = vin(7);   vtheta = vin(8);   vpsi = vin(9);
vVtas = vin(10); valpha = vin(11);  vbeta = vin(12);

valpha_L = vin(13); valpha_R = vin(14);

% Distance between the AoA vane and the center of gravity
global xva zva xvb zvb

% Distance between the AoA vane and the left wing location
global xva_L yva_L

% Distance between the AoA vane and the right wing location
global xva_R yva_R

% Define observation model
z(1) = x                                                        + vx; 
z(2) = y                                                        + vy;
z(3) = -zin                                                     + vh;
z(4) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*cos(psi) ...
        - (v*cos(phi) - w*sin(phi))*sin(psi) + Wxe              + vxdot; 
z(5) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*sin(psi) ...
        + (v*cos(phi) - w*sin(phi))*cos(psi) + Wye              + vydot; 
z(6) = u*sin(theta) - (v*sin(phi) + w*cos(phi))*cos(theta)      + vhdot;
z(7) = phi                                                      + vphi;
z(8) = theta                                                    + vtheta;
z(9) = psi                                                      + vpsi;
z(10) = sqrt(u^2 + v^2 + w^2)                                   + vVtas;
z(11) = (1+Caup)*atan(w/u) - xva*(q-Lq)/u                       + valpha;
z(12) = atan(v/u) + xvb*(r-Lr)/u - zvb*(p-Lp)/u                 + vbeta;
z(13) = (1+Caup)*atan(w_L/u_L) - (xva_L*(q-Lq)-yva_L*(p-Lp)) / (u + yva_L*(r-Lr))  + valpha_L;
z(14) = (1+Caup)*atan(w_R/u_R) - (xva_R*(q-Lq)-yva_R*(p-Lp)) / (u + yva_R*(r-Lr))  + valpha_R;
end