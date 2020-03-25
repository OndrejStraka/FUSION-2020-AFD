function [w,y,Pyy] = mynoise_GPB1_1D(xi,u,model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE PREDICTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xiPred_means = model.As'.*xi(1) +  model.Bs'*u;
xiPred_means_replicated = xiPred_means([1,1,2,2]);
xiPred_vars = model.As'.^2.*xi(2) + model.Gs'.^2;
%xiPred_weights = [model.P(:,1)*xi(5);model.P(:,2)*(1-xi(5))];
A =  model.P *[1 0;0 -1];
w = (A(:)*xi(3)+[0 0 model.P(1,2) model.P(2,2)]')';
% ybars: nModels^2 x length(y)
%ybars = y(ones(4,1),:)- model.Css'.*xiPred_means([1,1,2,2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASUREMENT PREDICTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = (model.Css'.*xiPred_means_replicated)';
% Pyys: nModels^2 x 1
Pyys = model.Css'.^2.*xiPred_vars([1 1 2 2]) + model.Hss'.^2;
Pyy = reshape(Pyys,1,1,4);
end
