function estimate = lssmmkfp_1D_GPB1(u,estimate,model,nodeStates)

  nNode = length(model);

  % Perform prediction for each node
  for iNode = 1:nNode
    % Get the number of models for node
    nM = size(model(iNode).P,1);
    % Perform prediction for each model
    estimate{iNode}.xpseq = [model(iNode).M(1).A*estimate{iNode}.xfseqfused model(iNode).M(2).A*estimate{iNode}.xfseqfused]  + [model(iNode).M(1).B model(iNode).M(2).B]*u{iNode};
    tmp = [model(iNode).M(1).A*estimate{iNode}.Pxxfseqfused*model(iNode).M(1).A' model(iNode).M(2).A*estimate{iNode}.Pxxfseqfused*model(iNode).M(2).A'] + [model(iNode).M(1).G*model(iNode).M(1).G' model(iNode).M(2).G*model(iNode).M(2).G'];
    estimate{iNode}.Pxxpseq = reshape(tmp,1,1,nM);
    estimate{iNode}.pmupseq =model(iNode).P*estimate{iNode}.pmufseq;
    % Perform merging to a single term
    tmp = estimate{iNode}.xpseq*estimate{iNode}.pmufseq;
    estimate{iNode}.Pxxpseq = (estimate{iNode}.Pxxpseq(:,:,1) + (estimate{iNode}.xpseq(:,1)-tmp)*(estimate{iNode}.xpseq(:,1)-tmp)')*estimate{iNode}.pmufseq(1)+...
                              (estimate{iNode}.Pxxpseq(:,:,2) + (estimate{iNode}.xpseq(:,2)-tmp)*(estimate{iNode}.xpseq(:,2)-tmp)')*estimate{iNode}.pmufseq(2);
    estimate{iNode}.xpseq = tmp;
  end
end
