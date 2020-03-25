function estimate = smmkfp(u,estimate,model)
  % SMMKFP prediction step of Switching Multiple Model estimation
  % algorithm
  %
  % ESTIMATE = SMMFILTER_P(U,ESTIMATE,MODELS) returns predictive mean values,
  % covariance matrices and sequence probabilities using filtering values in
  % ESTIMATE, model parameters in MODEL and input U 

  % Dimension of state and number of filtering estimates
  [nx,nXf] = size(estimate.xfseq);

  % Number of modes (models)
  nModels = length(model.M);

  % Set prediction values to zeros
  estimate.xpseq = zeros(nx,nXf);
  estimate.Pxxpseq = zeros(nx,nx,nXf);
  estimate.pmupseq = zeros(nModels*nXf,1);

  % Compute predictions
  idxModel = 0;
  for i = 1:nXf 
    idxModel = idxModel + 1;
    if idxModel > nModels
      idxModel = 1;
    end

    % Compute predictive mean values and covariances for model sequences
    [estimate.xpseq(:,i),estimate.Pxxpseq(:,:,i)] = ...
      kfp(u,estimate.xfseq(:,i),estimate.Pxxfseq(:,:,i),...
      model.M(idxModel).A,model.M(idxModel).B,model.M(idxModel).G*model.M(idxModel).G');

    % Compute the predictive probabilities of model sequences
    estimate.pmupseq((i-1)*nModels+1:i*nModels) = model.P(:,idxModel)*estimate.pmufseq(i);
  end
end
