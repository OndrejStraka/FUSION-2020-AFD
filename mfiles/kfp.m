function [xp,Pxxp] = kfp(u,xf,Pxxf,A,B,Sigmaw)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(u)
    xp = A*xf + B*u;
else
    xp = A*xf;
end
Pxxp = A*Pxxf*A' + Sigmaw;


end