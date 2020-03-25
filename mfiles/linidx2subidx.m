function subidx = linidx2subidx(siz,linidx,checkRanges)
%LINIDX2SUBIDX Compute subscript indices for linear indices
%   SUBIDX = LINIDX2SUBIDX(SIZ,LINIDX,CHECKRANGES) compute the subscript
%   indices SUBIDX that correspond to linear indices LINIDX for given
%   matrix size SIZ. It is a vectorized version of the Matlab function
%   IND2SUB where one matrix size in row vector SIZ is assumed but several
%   linear indices in a row vector LINIDX can be evaluated at once. The
%   output  SUBIDX is matrix whose columns contain corresponding subscript
%   indices. The individual rows corresponds to individual dimensions. The
%   relationship btween linear and subscript indices is: linidx =
%   subidx(1,:) + [1 siz(1) siz(1)*siz(2) siz(1)*siz(2)*siz(3) ... ] *
%   (subidx(2:end,:) - 1). Note that quite time demaning check of linear
%   indices is also performed, it can be switched off using optional
%   logical input argument CHECKRANGES.
%
%   Example
%   linidx = subidx2linidx([2 2],[1 1 2 2; 1 2 1 2])
%
%   See also IND2SUB, SUBIDX2LINIDX.
%
%   Copyright 2018 IDM (Ivo Punèocháø)

% The third input argument is optional
if nargin<3
    checkRanges = true;
end

% Get sizes of the input arguments
sizSiz = size(siz);
sizLinidx = size(linidx);

% Check the consistency of the input arguments
if numel(sizSiz)~=2 || sizSiz(1)~=1
    exception = MException('MyToolbox:subidx2linidx:wrongDim',...
        'The first argument must be a row vector with matrix size');
    throw(exception);
end

if numel(sizLinidx)~=2 || sizLinidx(1)~=1
    exception = MException('MyToolbox:subidx2linidx:wrongDim',...
        'The second argument must be a row vector a linear indices');
    throw(exception);
end

if checkRanges && (any(linidx<1) || any(linidx>prod(siz))) % FIXME this check takes a lot of time
    exception = MException('MyToolbox:subidx2linidx:outOfRange',...
        'A linear index in the second input argument is out of range');
    throw(exception);
end

% Prealocate array for subscript indices
subidx = zeros(sizSiz(2),sizLinidx(2));

% Compute multiples for individual dimensions
dimMultiples = [1 cumprod(siz(1:end-1))];

linidx = linidx - 1;
for i = sizSiz(2):-1:1
    subidx(i,:) =  fix(linidx/dimMultiples(i));
    linidx = linidx - dimMultiples(i)*subidx(i,:);
end
subidx = subidx + 1;


end