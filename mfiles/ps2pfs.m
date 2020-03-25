function [pmuseqnew, idx] = ps2pfs(pmuseq,n)
% PS2PFS - computes the probabilities of sequences that have common final
% model sequence
%
% [pseqnew, idx_max] = ps2pfs(pmuseq,n)
%
% PSEQ is row vector of probabilities of individiual sequences
% NSEQREQ is required number of resulting sequences


% Number of model sequences
nseq = length(pmuseq);

% Check that the number of required sequences is the same or lower than
% current number of sequnces
if n < nseq
    % Reshape the column vector of probabilities to get matrix P, where each column of
    % P contains probabilities of model sequences that have the same final
    % sequence
    pmuseqReshaped = reshape(pmuseq,n,[]);
    
    % The probabilities of model sequences that have the same final sequence
    pmuseqnew = sum(pmuseqReshaped,2);
    
    % Index of model sequence that has maximum probability within group of sequences with common final model sequence 
    [~,idx] = max(pmuseqReshaped,[],2);
    idx = (idx-1)*n + (1:n)';
else
    pmuseqnew = pmuseq;
    [~,idx] = max(pmuseq,[],2);
end

