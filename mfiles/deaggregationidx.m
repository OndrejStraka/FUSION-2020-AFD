function x = deaggregationidx(xqtlvl,subidx)
% DEAGGREGATIONIDX creates vectors from quantization levels
%   X = DEAGGREGATION(XQTLVL,SUBIDX) creates 2D array X of column vectors using
%   quantization levels in row cell array XQTLVL and indices SUBIDX

% The dimension and the number of the vectors to be created
sizsubidx = size(subidx);

if numel(sizsubidx) > 2
exception = MException('MyToolbox:deaggregationidx:inconsistentDimension',...
        'The SUBIDX must be at most two dimensional array.');
    throw(exception);
end

% Check that dimensions are consistent
if sizsubidx(1) ~= length(xqtlvl)
    exception = MException('MyToolbox:deaggregationidx:inconsistentDimension',...
        'The dimensions of state and subscript must be the same.');
    throw(exception);
end

% Preallocate array
x = zeros(sizsubidx);

% Create states using subscripts
for i = 1:sizsubidx(1)
    x(i,:) = xqtlvl{i}(subidx(i,:));
end