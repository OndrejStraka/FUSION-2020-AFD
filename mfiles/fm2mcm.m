function [mfm,Pfm] = fm2mcm(w,m,P)
%FM2MCM Compute mean and covariance matrix of finite mixture
%   Detailed explanation goes here

% Get sizes of input arguments
sizw = size(w);
sizm = size(m);
sizP = size(P);

if numel(sizw)~=2 || sizw(1)~=1
    error('w must be a row')
end

if numel(sizm)~=2 || sizm(2)~=sizw(2)
    error('m must contain means as columns')
end

if numel(sizP)~=3 || sizP(1)~=sizP(2) || sizP(1)~=sizm(1) || sizP(3)~=sizw(2)
    error('P has inconsistent dimension')
end

% Compute the mean value of the finite mixture
mfm = sum(bsxfun(@times,m,w),2);

% Compute differences between individual means and the mean of the finite mixture
tmp = bsxfun(@minus,m,mfm);

% Compute the covariance matrix of the finite mixture
Pfm = zeros(sizP(1));
for i = 1:length(w)
    Pfm = Pfm + w(i)*(P(:,:,i) + tmp(:,i)*tmp(:,i)');
end

% Loopless computation of the the covariance matrix - it is slower than
% loop for finite mixture with a few terms
%y = bsxfun(@times,w,tmp);
%Pfm = sum(bsxfun(@times,P,reshape(w,1,1,sizw(2))),3) + y*tmp'


end