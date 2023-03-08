% Filename :  TbE.m
%
% Transforms from the navigation (i.e. Earth) frame to the body frame
%
% Inputs: 
%     phi:        Roll angle in radians
%     theta:      Pitch angle in radians
%     psi:        Yaw angle in radians
%
% Outputs:
%     Tbe:        transformation matrix [3x3]
%
% Example:  [Tbe] = TbE(0.4, 0.65, 0.25)
%
% Made by: M.A. van den Hoek & L.J. van Horssen, September 2016
%% Transformation matrix from ECEF to body-frame
function [Tbe] = TbE(phi,theta,psi)
% Transformation matrix
Tbe = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
      sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);
      cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)];
end
  

