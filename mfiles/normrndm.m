function x = normrndm(mx,Px,m)
% NORMRNDM Random sample from normal distribution
%
% X = NORMRNDM(MX,PX,M) generates random vectors from normal distribution
% with mean value MX and positive definite covariance matrix PX.

if nargin < 3
    m = 1;
end

% Dimension of x
nx = length(mx);

% Compute Cholesky decomposition
sPx = chol(Px)';

% Generate sample vectors
x = sPx*randn(nx,m) + mx(:,ones(1,m));