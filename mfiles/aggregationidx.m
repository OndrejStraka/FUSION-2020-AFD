function subidx = aggregationidx(x,xqtie)
% AGGREGATIONIDX projects nD vectors onto a regular grid in nD space
%
% SUBIDX = AGGREGATIONIDX(x,innerEdges) projects vectors x onto a regular
% grid with inner edges innerEdges and returns subscript indexes

% The dimension and the number of the state vectors
sizx = size(x);

if numel(sizx) > 2 || any(sizx == 0)
    exception = MException('MyToolbox:aggregationidx:InconsistentDimensions',...
        'The input X must be at most a two dimensional array');
    throw(exception);
end

% Check that the dimensions of the state and inner edges are the same
if  sizx(1) ~= length(xqtie)
    exception = MException('MyToolbox:aggregationidx:InconsistentDimensions',...
        'The dimension of the state and inner edges must be the same');
    throw(exception);
end

% Preallocate array for subscript indexes
subidx = zeros(sizx);

% Computes subscript indexes for each dimension
% NOTE: HISTC uses comparison edges(k)<= x < edge(k+1). For more details
% see help to HISTC
% for i = 1:nx
%     [~,subidx(i,:)] = histc(x(i,:),[-inf xqtie{i} inf]);
% end

% Since HISTC is not recommended in new versions of MATLAB, function
% DISCRETIZE is used instead
for i = 1:sizx(1)
   subidx(i,:) = discretize(x(i,:),[-inf xqtie{i} inf]);
end