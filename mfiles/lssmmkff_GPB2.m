function [estimate,estimateLSS] = lssmmkff_GPB2(y,estimate,estimateLSS,model,LSS,nodeStates)
  % nM means nM vars nM^2 weights for each node -> nM means nM vars nM weights for each node

  % Get the number of nodes
  nNode = length(y);

  % Get the dimension of the state
  nx = size(nodeStates,1);
  
  % Detect first step according to 1st model
  first_step = (size(estimate{1}.xpseq,2)==1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERING (FOR EACH NODE INDIVIDUALLY)
  for iNode = 1:nNode
    estimate{iNode} = smmkff_GPB2(y{iNode},estimate{iNode},model(iNode));
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTRAL NODE PROCESSING
  if first_step
    estimateLSS{1}.likelihood = kron(estimate{1}.likelihood,estimate{2}.likelihood); % combine local likelihoods
    estimateLSS{1}.pmufseq = estimateLSS{1}.pmupseq .* estimateLSS{1}.likelihood; % multiply by likelihood
    estimateLSS{1}.pmufseq = estimateLSS{1}.pmufseq/sum(estimateLSS{1}.pmufseq); % normalize
  else
    estimateLSS{1}.likelihood_full = [kron(estimate{1}.likelihood(1:2),estimate{2}.likelihood(1:2));kron(estimate{1}.likelihood(1:2),estimate{2}.likelihood(3:4));...
      kron(estimate{1}.likelihood(3:4),estimate{2}.likelihood(1:2));kron(estimate{1}.likelihood(3:4),estimate{2}.likelihood(3:4))]; % combine local likelihoods
    estimateLSS{1}.pmufseq_full = estimateLSS{1}.pmupseq_full.*estimateLSS{1}.likelihood_full; % multiply by likelihood
    estimateLSS{1}.pmufseq_full = estimateLSS{1}.pmufseq_full/sum(estimateLSS{1}.pmufseq_full); % normalize
    estimateLSS{1}.pmufseq = kron(ones(1,4),eye(4))*estimateLSS{1}.pmufseq_full; % merging
  end
  % calculate probabilties of local nodes
  estimate{1}.pmufseq_central = kron(eye(2),ones(1,2))*estimateLSS{1}.pmufseq; %[1 1 0 0;0 0 1 1] * [(1,1) (1,2) [2,1) (2,2)]'
  estimate{2}.pmufseq_central = kron(ones(1,2),eye(2))*estimateLSS{1}.pmufseq; %[1 0 1 0;0 1 0 1] * [(1,1) (1,2) [2,1) (2,2)]'

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MERGING (FOR EACH NODE INDIVIDUALLY)
  if first_step % in the first step do not merge, only feed back if necessary
    if LSS.returnPmuToLocalnodes % if the probabilities are fed back
      for iNode = 1:nNode
        estimate{iNode}.pmufseq = estimate{iNode}.pmufseq_central;
      end
    end
  else
    if LSS.returnPmuToLocalnodes % if the probabilities are fed back
      estimate{1}.pmufseq = kron(eye(2),kron(ones(1,2),kron(eye(2),ones(1,2))))*estimateLSS{1}.pmufseq_full; % obtain local weights of node 1
      estimate{2}.pmufseq = kron(ones(1,2),kron(eye(2),kron(ones(1,2),eye(2))))*estimateLSS{1}.pmufseq_full; % obtain local weights of node 2
    end
    for iNode = 1:nNode
      [estimate{iNode}.pmufseq,estimate{iNode}.xfseq,estimate{iNode}.Pxxfseq] = mmmerge(estimate{iNode}.pmufseq,estimate{iNode}.xfseq,estimate{iNode}.Pxxfseq,size(model(iNode).P,1));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXCHANGE
  %
  % Nodes exchange their local estimates and each computes the same set of
  % means and convariance matrices using covariance intersection algorithm
  % The code is written for two nodes only

  % Get the number of filtering estimates for individual nodes
  nM = [size(estimate{1}.xfseq,2);size(estimate{2}.xfseq,2)];

  Pci = zeros(nx,nx,nM(1),nM(2));
  xci = zeros(nx,nM(1),nM(2));
  % Compute covariance intersection for every model combination
  for iModel = 1:nM(1)
    for jModel = 1:nM(2)
      % Create state
      xci(nodeStates(:,1),iModel,jModel) = estimate{1}.xfseq(:,iModel);
      xci(nodeStates(:,2),iModel,jModel) = estimate{2}.xfseq(:,jModel);

      Pci(nodeStates(:,1),nodeStates(:,1),iModel,jModel) = estimate{1}.Pxxfseq(:,:,iModel);
      Pci(nodeStates(:,2),nodeStates(:,2),iModel,jModel) = estimate{2}.Pxxfseq(:,:,jModel);
    end
  end

  % Multiply the covariance matrices to account for unknown correlation
  Pci = 2*Pci;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUSION

  % Gaussian sum to single gaussian for node 1
  estimate{1}.xfseqfused = zeros(nx,nM(1));
  estimate{1}.Pxxfseqfused = zeros(nx,nx,nM(1));
  for iModel = 1:nM(1)
    [estimate{1}.xfseqfused(:,iModel),estimate{1}.Pxxfseqfused(:,:,iModel)] = fm2mcm(estimate{2}.pmufseq',squeeze(xci(:,iModel,:)),squeeze(Pci(:,:,iModel,:)));
  end

  % Gaussian sum to single gaussian for node 2
  estimate{2}.xfseqfused = zeros(nx,nM(2));
  estimate{2}.Pxxfseqfused = zeros(nx,nx,nM(2));
  for jModel = 1:nM(2)
    [estimate{2}.xfseqfused(:,jModel),estimate{2}.Pxxfseqfused(:,:,jModel)] = fm2mcm(estimate{1}.pmufseq',squeeze(xci(:,:,jModel)),squeeze(Pci(:,:,:,jModel)));
  end
