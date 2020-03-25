function estimate = lssmmkff_GPB1(y,estimate,model,nodeStates)
  % Get the number of nodes
  nNode = length(y);

  % Get the dimension of the state
  nx = size(nodeStates,1);

  % Perform filtering step for each node individualy
  for iNode = 1:nNode
    if length(estimate{iNode}.xpseq) == 1
      estimate{iNode} = smmkff_GPB2(y{iNode},estimate{iNode},model(iNode));
    else
      estimate{iNode} = smmkff_GPB1(y{iNode},estimate{iNode},model(iNode));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXCHANGE
  %
  % Nodes exchange their local estimates and each computes the same set of
  % means and convariance matrices using covariance intersection algorithm
  % The code is written for two nodes only

  % Get the number of filtering estimates for individual nodes
  nModel = [size(estimate{1}.xfseq,2);size(estimate{2}.xfseq,2)];

  %Pci = zeros(nx,nx,nModel(1),nModel(2));
  %xci = zeros(nx,nModel(1),nModel(2));
  % Compute covariance intersection for every model combination
  for i = 1:nModel(1)
    for j = 1:nModel(2)
      % Create state
      xci(nodeStates(:,1),i,j) = estimate{1}.xfseq(:,i);
      xci(nodeStates(:,2),i,j) = estimate{2}.xfseq(:,j);

      Pci(nodeStates(:,1),nodeStates(:,1),i,j) = estimate{1}.Pxxfseq(:,:,i);
      Pci(nodeStates(:,2),nodeStates(:,2),i,j) = estimate{2}.Pxxfseq(:,:,j);
    end
  end

  % Multiply the covariance matrices to account for unknown correlation
  Pci = 2*Pci;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUSION


  % Gaussian sum to single gaussian for node 1 for each model
  estimate{1}.xfseqfused = zeros(nx,nModel(1));
  estimate{1}.Pxxfseqfused = zeros(nx,nx,nModel(1));
  for iModel = 1:nModel(1)
    [estimate{1}.xfseqfused(:,iModel),estimate{1}.Pxxfseqfused(:,:,iModel)] = fm2mcm(estimate{2}.pmufseq',squeeze(xci(:,iModel,:)),squeeze(Pci(:,:,iModel,:)));
  end

  % Gaussian sum to single gaussian for node 2 for each model
  estimate{2}.xfseqfused = zeros(nx,nModel(2));
  estimate{2}.Pxxfseqfused = zeros(nx,nx,nModel(2));
  for jModel = 1:nModel(2)
    [estimate{2}.xfseqfused(:,jModel),estimate{2}.Pxxfseqfused(:,:,jModel)] = fm2mcm(estimate{1}.pmufseq',squeeze(xci(:,:,jModel)),squeeze(Pci(:,:,:,jModel)));
  end


  % Gaussian sum to single gaussian for each node
  for iNode = 1:nNode
    [estimate{iNode}.xfseqfused,estimate{iNode}.Pxxfseqfused] = fm2mcm(estimate{iNode}.pmufseq',squeeze(estimate{iNode}.xfseqfused),squeeze(estimate{iNode}.Pxxfseqfused));
  end
end
