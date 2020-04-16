function estimate = lssmmkfp_1D_GPB2(u,estimate,model,nodeStates)
  % nM means nM vars nM weights for each node -> nM means nM vars nM^2 weights for each node

  nNode = length(model);

  % Perform prediction for each node
  for iNode = 1:nNode
    % Get the number of models for node
    nM = size(model(iNode).P,1);
    estimate{iNode}.xpseq = [model(iNode).M(1).A*estimate{iNode}.xfseqfused(:,1) model(iNode).M(2).A*estimate{iNode}.xfseqfused(:,2)]  + [model(iNode).M(1).B model(iNode).M(2).B]*u{iNode};
    tmp = [model(iNode).M(1).A*estimate{iNode}.Pxxfseqfused(:,:,1)*model(iNode).M(1).A' model(iNode).M(2).A*estimate{iNode}.Pxxfseqfused(:,:,2)*model(iNode).M(2).A'] + [model(iNode).M(1).G*model(iNode).M(1).G' model(iNode).M(2).G*model(iNode).M(2).G'];
    estimate{iNode}.Pxxpseq = reshape(tmp,1,1,nM);
    tmp = model(iNode).P .* estimate{iNode}.pmufseq';
    estimate{iNode}.pmupseq = tmp(:);
  end
end
