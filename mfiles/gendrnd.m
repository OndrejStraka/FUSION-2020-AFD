function x = gendrnd(p,varargin)
% GENDRND  Generate random array using user defined discrete distribution.
%
% X = GENDRND(p,m,n,values,cdfTol) returns an M-by-N matrix containing
% pseudorandom values drawn from the set VALUES (default set {1, 2, 3, ...,
% K-1, K}) with user defined probability ditribution P over the set {1, 2,
% 3, ..., K-1, K}.

sizp = size(p);
if numel(sizp)>2 || sizp(2)>1 || ~isreal(p)
    exception = MException('MyToolbox:gendrnd:wrongReqInArgs',...
        'Argument P has to be a column vector of real numbers.');
    throw(exception);
end

% Check the maximum number of optional input arguments
maxnumvarargs = 4;
numvarargs = length(varargin);
if numvarargs > maxnumvarargs
    exception = MException('MyToolbox:gendrnd:tooManyOptInArgs',...
        'The function accepts at most %d optional input arguments.',maxnumvarargs);
    throw(exception);
end

% Skip any new inputs if they are empty
newVals = cellfun(@(x) isnotempty(x), varargin);

% Set defaults for optional input arguments
optargs = {1,1,(1:length(p))',sqrt(eps)};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);

% Place optional input arguments in memorable variable names
[m,n,values,cdfTol] = optargs{:};

% Compute cumulative distribution
cdfp = [0;cumsum(p)];

% Check that the conditions for p to be a probability vector are satisfied
if any(p < 0) || abs(cdfp(end)-1) > cdfTol
    exception = MException('MyToolbox:gendrnd:notValidDistribution',...
        'Argument P has to contain nonnegative real numbers that sum up to 1 to be a valid pmf.');
    throw(exception);
end

% Generates random numbers from the open interval (0,1)
randomNumbers = rand(m*n,1);

% Find the indices
%[~,~,idx] = histcounts(randomNumbers,cdfp);
idx = discretize(randomNumbers,cdfp); % Discretize seems to be slightly faster than histcounts

% Index the values and reshape them
x = values(idx);
x = reshape(x,m,n);


end