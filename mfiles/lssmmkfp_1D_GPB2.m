function [estimate,estimateLSS] = lssmmkfp_1D_GPB2(u,estimate,estimateLSS,model,LSS,nodeStates)
  % nM means nM vars nM weights for each node -> nM means nM vars nM^2 weights for each node

  nNode = length(model);

  % Perform prediction of local nodes
  for iNode = 1:nNode
    % Get the number of models of node
    nM = size(model(iNode).P,1);
    % Mean
    estimate{iNode}.xpseq = [model(iNode).M(1).A*estimate{iNode}.xfseqfused(:,1) model(iNode).M(2).A*estimate{iNode}.xfseqfused(:,2)]  + [model(iNode).M(1).B model(iNode).M(2).B]*u{iNode};
    % Variance
    tmp = [model(iNode).M(1).A*estimate{iNode}.Pxxfseqfused(:,:,1)*model(iNode).M(1).A' model(iNode).M(2).A*estimate{iNode}.Pxxfseqfused(:,:,2)*model(iNode).M(2).A'] + [model(iNode).M(1).G*model(iNode).M(1).G' model(iNode).M(2).G*model(iNode).M(2).G'];
    estimate{iNode}.Pxxpseq = reshape(tmp,1,1,nM);
    % probabilities
    tmp = model(iNode).P .* estimate{iNode}.pmufseq';
    estimate{iNode}.pmupseq = tmp(:);
  end
  % Perform prediction of central node
  estimateLSS{1}.pmupseq_full = LSS.P(:).* kron(estimateLSS{1}.pmufseq,ones(4,1)); % predict 16 probabilities
end
