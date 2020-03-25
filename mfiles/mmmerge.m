function [pmuseqmerged,xseqmerged,Pxxseqmerged] = mmmerge(pmuseq,xseq,Pxxseq,nXmerged)
%MMMERGE merge model sequences
%   Detailed explanation goes here

% The dimension of the state  and the number of model sequences
[nx,nX] = size(xseq);

if nXmerged < nX
    % Compute probabilities of model sequences after merging
    pmuseqmerged = ps2pfs(pmuseq,nXmerged);
    
    % Compute the mean values of the merged sequences
    xseqmerged = reshape(sum(reshape(bsxfun(@times,xseq,pmuseq'),nXmerged*nx,[]),2),nx,nXmerged);
    xseqmerged = bsxfun(@rdivide,xseqmerged,pmuseqmerged');
    
    % Replace nan that results from zero division by zeros
    xseqmerged(isnan(xseqmerged)) = 0;
    
    
    % Compute the covariance matrices of the merged sequences
    % Preallocate aray
    Pw = zeros(nx,nx,nX);
    % Preallocate aray
    Pxxseqmerged = zeros(nx,nx,nXmerged);
    
    dxseq = xseq - repmat(xseqmerged,1,nX/nXmerged);
    for i = 1:nX
        Pw(:,:,i) = (Pxxseq(:,:,i) + dxseq(:,i)*dxseq(:,i)')*pmuseq(i);
    end
    
    for i = 1:nXmerged
        if pmuseqmerged(i) ~= 0
            Pxxseqmerged(:,:,i) = sum(Pw(:,:,i:nXmerged:end),3)/pmuseqmerged(i);
        else
            Pxxseqmerged(:,:,i) = eye(nx,nx);
        end
    end
else
    pmuseqmerged = pmuseq;
    xseqmerged = xseq;
    Pxxseqmerged = Pxxseq;
end

end