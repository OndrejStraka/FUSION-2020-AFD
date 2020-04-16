function estimate = smmkff_GPB2(y,estimate,model,computeOnlyPredictiveOutput)
  % SMMKFF filtering step of Switching Multiple Model estimation
  % algorithm
  %
  % ESTIMATE = SMMFILTER_F(Y,ESTIMATE,MODEL) returns filtering mean values,
  % covariance matrices and sequence probabilities computed using predictive
  % values in ESTIMATE, model parameters in MODEL and measuremnt Y

  % !!!!!!!!!!!!!!!!
  % Square root version should be used in future - needs to be done
  % !!!!!!!!!!!!!!!!!!!!!!!!


  if nargin<4
    computeOnlyPredictiveOutput = false;
  end

  % Dimension of state and number of prediction estimates
  [nx,nXp] = size(estimate.xpseq);

  % Dimension of output
  ny = size(model.M(1).C,1);

  % Number of modes (models)
  nModels = length(model.M);

  % Determined the number of filtering estimates
  nXf = nModels*nXp;

  % Set filtering values to zeros
  estimate.xfseq = zeros(nx,nXf);
  estimate.Pxxfseq = zeros(nx,nx,nXf);
  estimate.pmufseq = zeros(nXf,1);
  estimate.likelihood = zeros(nXf,1);

  % Set predictive output values to zeros
  estimate.ypseq = zeros(ny,nXf);
  estimate.Pyypseq = zeros(ny,ny,nXf);

  for i = 1:nXp
    for j = 1:nModels
      % Precompute index to shorten notation and speed up code
      idx = nModels*(i-1)+j;
      [estimate.xfseq(:,idx),estimate.Pxxfseq(:,:,idx),...
        estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx)] = ...
        kff(y,estimate.xpseq(:,i),estimate.Pxxpseq(:,:,i),model.M(j).C,model.M(j).H*model.M(j).H');

      % Compute non-normalized filtering probability of sequence
      % it's not numerically robust :(
      % This is computed only for the case when the central node does not feed the probabilities back
      estimate.pmufseq(idx) = normpdfm(y,estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx)) * estimate.pmupseq(idx);
      estimate.likelihood(idx) = normpdfm(y,estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx));
    end
  end

  if computeOnlyPredictiveOutput
    return
  end

  % The weights that are less or equal to EPS(max_weight) are set to zero
  estimate.pmufseq(estimate.pmufseq<=eps(max(estimate.pmufseq))) = 0;
  % Normalize
  sumpmufseq = sum(estimate.pmufseq);
  estimate.pmufseq = estimate.pmufseq/sumpmufseq;
