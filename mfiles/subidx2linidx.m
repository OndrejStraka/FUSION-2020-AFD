function linidx = subidx2linidx(siz,subidx,checkRanges)
%SUBIDX2LINIDX Computes linear indices to given subscript indices
%   It is a vectorized version of the Matlab function sub2ind. The first
%   input argument SIZ should be a row vector of matrix size. The second
%   input argument SUBIDX is a matrix whose columns are subscript indices.
%   The output is a row vector of linear indices that corresponds to the
%   subscript indices. The relationship btween linear and subscript indices
%   is: linidx = subidx(1,:) + [1 siz(1) siz(1)*siz(2) siz(1)*siz(2)*siz(3)
%   ... ] * (subidx(2:end,:) - 1).
%
%   Example
%   linidx = subidx2linidx([2 2],[1 1 2 2; 1 2 1 2])
%
%   See also SUB2IND, LINIDX2SUBIDX.
%
%   Copyright 2018 IDM (Ivo Punèocháø)

% The third input argument is optional
if nargin<3
    checkRanges = true;
end


% Get sizes of the input arguments
sizSiz = size(siz);
sizSubidx = size(subidx);

% Check the consistency of the input arguments
if numel(sizSiz)>2 || sizSiz(1)~=1
    exception = MException('MyToolbox:subidx2linidx:wrongDim',...
        'The first argument must be a row vector');
    throw(exception);
end

if numel(sizSubidx)>2 || sizSiz(2)~=sizSubidx(1)
    exception = MException('MyToolbox:subidx2linidx:wrongDim',...
        'The second argument must be a matrix that has the same number of rows as the first argument columns');
    throw(exception);
end

if checkRanges && (any(any(subidx<1)) || any(any(bsxfun(@lt,siz',subidx)))) % FIXME this check takes a lot of time
    exception = MException('MyToolbox:subidx2linidx:outOfRange',...
        'A subindex in the second input argument is out of range');
    throw(exception);
end

% Compute multiples for individual dimensions
dimMultiples = [1 cumprod(siz(1:end-1))];

% Compute the linear index
linidx = 1 + dimMultiples*(subidx-1);


end