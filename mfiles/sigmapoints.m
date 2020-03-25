function [X,w,zmX] = sigmapoints(mx,Pxx,varargin)
%
%   sigmapoints(mx,Pxx,kappa,sqrtPxx)

% The dimension of the random variable
nx = size(mx,1);

% Accept only two optional inputs at most
numvarargs = length(varargin);
if numvarargs > 2
    error('mytoolbox:sigmapoints:TooManyInputs', ...
        'The fuction  requires at most 2 optional inputs: the scaling parameter kappa and square root matrix');
end

% Set defaults for optional inputs
optargs = {max(0,3 - nx)  []};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[kappa, sqrtPxx] = optargs{:};

% Is empty, compute the value for parameter sqrtPxx
if isempty(sqrtPxx)
    sqrtPxx = chol(Pxx,'lower');
end

% Compute the number of sigma points
nSigmaPoints = 2*nx + 1;

% Compute the coefficient c
c = sqrt(nx + kappa);

% Create zero mean sigma points
zmSigmaPoints = [zeros(nx,1) sqrtPxx -sqrtPxx];

% Compute sigma points
zmX = c*zmSigmaPoints;
X = mx(:,ones(1,nSigmaPoints)) + zmX;

% Compute weights of sigma points
w = [kappa 0.5*ones(1,nSigmaPoints-1)]/(nx + kappa);