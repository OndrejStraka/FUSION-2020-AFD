function [xf,Pxxf,yp,Pyyp] = kff(y,xp,Pxxp,C,Sigmav)
%KFF Filtering step of the Kalman filter
%   [xf,Pxxf,yp,Pyyp] = kff(y,xp,Pxxp,C,Sigmav) computes filtering mean
%   value 'xf', filtering covariance matrix 'Pxxf', predictive mean value
%   'yp', and predictive covariance matrix 'Pyyp' based on the measurement
%   'y', predictive mean value 'xp', predictive covariance matrix 'Pxxp',
%   measurement gain matrix 'C' and noise covariance matrix 'Sigmav'

% Compute the predictive mean value and covariance matrix of output
yp = C*xp;
Pyyp = C*Pxxp*C' + Sigmav;

% Compute the Kalman gain
K = Pxxp*C'/Pyyp;

% Compute the filtering mean value and covariance matrix
xf = xp + K*(y-yp);
Pxxf = Pxxp - K*Pyyp*K'; % helps keep symmetry

% Basic - computationally efficient but suffers from numerical problems
% nx = size(xp,1);
% Pxxf = (eye(nx)-K*C)*Pxxp;

% Joseph form - helps keep symmetry and is valid for nonoptimal gain
% nx = size(xp,1);
% tmp = eye(nx) - K*C;
% Pxxf = tmp*Pxxp*tmp' + K*Sigmav*K';


end