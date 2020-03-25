function f = normpdfm(x,m,Sigma)
%NORMPDFM Normal probability density function for a vector argument
%   F = NORMPDFM(X,M,SIGMA) returns the pdf of normal distribution with
%   mean value M and covariance matrix SIGMA, evaluated at the values in X.
%   Either several vectors X are given, or several mean values M are given.
%   Only one covariance matrix SIGMA is expected.

% The dimension of x and the number of points to be evaluated
sizx = size(x);
sizm = size(m);
sizSigma = size(Sigma);

% Check dimensions
if numel(sizx)~=2 || numel(sizm)~=2 || numel(sizSigma)~=2 || ...
        sizx(1)~=sizm(1) || sizx(1)~=sizSigma(1) || sizSigma(1)~=sizSigma(2) || ...
        (sizx(2)>1 && sizm(2)>1)
    exception = MException('MyToolbox:normpdfm:wrongDimInArg',...
        'The input arguments have inconsistent dimensions');
    throw(exception);
end

cx = bsxfun(@minus,x,m);
md = -0.5*sum(cx.*(Sigma\cx),1);
f = (1/(sqrt(det(Sigma)*(2*pi)^sizx(1))))*exp(md);


end