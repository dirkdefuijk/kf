% Filename :  NavigationSystem.m
%
% Constructs the navigation system and calculates the system dynamics
% equation f(x,u,t)
%
% Inputs: 
%     t:        time instant
%     xin:      state vector x:
%               x = [x, y, z, u, v, w, phi, theta, psi, Wxe, Wye, Caup,
%                    lambda_x, lambda_y, lambda_z, lambda_p, lambda_q, 
%                    lambda_r, alpha_v]^T
%
%     uin:      input vector u:
%               u = [Ax, Ay, Az, p, q, r]^T
%   
%     win:      noise vector w:
%               w = [wx, wy, wz, wp, wq, wr]^T
%     
% Outputs:
%     xdot:     time derivative of the states        
%
% Made by: M.A. van den Hoek and L.J. van Horssen, September 2016
% Update by: J.B. van Ingen, November 2016
function xdot = NavigationSystem(t,xin, uin, win)
% Define the gravitational acceleration on Earth
g = 9.80665;

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

wx = win(1); wy = win(2); wz = win(3);
wp = win(4); wq = win(5); wr = win(6);

% Define navigation system
xdot(1) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*cos(psi) ...
            - (v*cos(phi) - w*sin(phi))*sin(psi) + Wxe; 
xdot(2) = (u*cos(theta) + (v*sin(phi) + w*cos(phi))*sin(theta))*sin(psi) ...
            + (v*cos(phi) - w*sin(phi))*cos(psi) + Wye; 
xdot(3) = -u*sin(theta) + (v*sin(phi) + w*cos(phi))*cos(theta);
xdot(4) = (Ax-Lx-wx) - g*sin(theta) + (r-Lr-wr)*v - (q-Lq-wq)*w;
xdot(5) = (Ay-Ly-wy) + g*cos(theta)*sin(phi) + (p-Lp-wp)*w - (r-Lr-wr)*u;
xdot(6) = (Az-Lz-wz) - g + g*cos(theta)*cos(phi) + (q-Lq-wq)*u - (p-Lp-wp)*v;
xdot(7) = (p-Lp-wp) + (q-Lq-wq)*sin(phi)*tan(theta) ...
            + (r-Lr-wr)*cos(phi)*tan(theta);
xdot(8) = (q-Lq-wq)*cos(phi) - (r-Lr-wr)*sin(phi);
xdot(9) = (q-Lq-wq)*sin(phi)/cos(theta) + (r-Lr-wr)*cos(phi)/cos(theta);
xdot(10) = 0;
xdot(11) = 0;
xdot(12) = 0;
xdot(13) = 0;
xdot(14) = 0;
xdot(15) = 0;
xdot(16) = randn*0.01;
xdot(17) = randn*0.01;
xdot(18) = 0.01*randn*pi/180;

xdot(19) = (Ax-Lx-wx) - g*sin(theta) + (r-Lr-wr)*v_L - (q-Lq-wq)*w_L;
xdot(20) = (Ax-Lx-wx) - g*sin(theta) + (r-Lr-wr)*v_R - (q-Lq-wq)*w_R;

xdot(21) = (Ay-Ly-wy) + g*cos(theta)*sin(phi) + (p-Lp-wp)*w_L - (r-Lr-wr)*u_L;
xdot(22) = (Ay-Ly-wy) + g*cos(theta)*sin(phi) + (p-Lp-wp)*w_R - (r-Lr-wr)*u_R;

xdot(23) = (Az-Lz-wz) - g + g*cos(theta)*cos(phi) + (q-Lq-wq)*u_L - (p-Lp-wp)*v_L;
xdot(24) = (Az-Lz-wz) - g + g*cos(theta)*cos(phi) + (q-Lq-wq)*u_R - (p-Lp-wp)*v_R;
end