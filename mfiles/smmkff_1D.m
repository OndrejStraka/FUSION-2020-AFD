function estimate = smmkff_nx1nz1(y,estimate,model,computeOnlyPredictiveOutput)
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

  % Set predictive output values to zeros
  estimate.ypseq = zeros(ny,nXf);
  estimate.Pyypseq = zeros(ny,ny,nXf);

  myp = zeros(nXf,1);
  srdetPyp = zeros(nXf,1);

  for i = 1:nXp
    for j = 1:nModels
      % Precompute index to shorten notation and speed up code
      idx = nModels*(i-1)+j;
      [estimate.xfseq(:,idx),estimate.Pxxfseq(:,:,idx),...
        estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx)] = ...
        kff(y,estimate.xpseq(:,i),estimate.Pxxpseq(:,:,i),model.M(j).C,model.M(j).H*model.M(j).H');

      % Compute non-normalized filtering probability of sequence
      % it's not numerically robust :(
      estimate.pmufseq(idx) = ...
        normpdfm(y,estimate.ypseq(:,idx),estimate.Pyypseq(:,:,idx))*estimate.pmupseq(idx);

      % Prepare some values used in numerically robust computation of
      % filtering sequence probabilities
      e = y-estimate.ypseq(:,idx);
      myp(idx) = e'*(estimate.Pyypseq(:,:,idx)\e);
      srdetPyp(idx) = sqrt(det(estimate.Pyypseq(:,:,idx))); 
    end
  end
  Pxxps = squeeze(estimate.Pxxpseq)';
  tmp = [model.M(1).C model.M(2).C]' * estimate.xpseq;
  estimate.ypseq = tmp(:)';
  ypbars = y - estimate.ypseq;
  Pyyps = [model.M(1).C model.M(2).C]'.^2 .* Pxxps + [model.M(1).H model.M(2).H]'.^2 ;
  estimate.Pyypseq = reshape(Pyyps(:),1,1,nXf);
  Pxy = [model.M(1).C model.M(2).C]' .* Pxxps;
  KGs = Pxy(:)'./Pyyps(:)';
  replicate = ones(nModels,1)* [1:nXp];
  estimate.xfseq = estimate.xpseq(:,replicate(:)) + KGs.*ypbars;
  estimate.Pxxfseq = reshape(Pxxps(:,replicate(:)) - KGs.*Pxy(:)',1,1,nXf);
  estimate.pmufseq = 1./sqrt(2*pi*Pyyps(:)).*exp(-0.5*ypbars'.^2./Pyyps(:)).*estimate.pmupseq;



  if computeOnlyPredictiveOutput
    return
  end

  % Perform weights normalization to obtain probabilities
  % The weights are normalized to sum up 1
  % The weights that are less or equal to EPS(max_weight) are set to zero
  estimate.pmufseq(estimate.pmufseq<=eps(max(estimate.pmufseq))) = 0;

  % Normalize
  sumpmufseq = sum(estimate.pmufseq);
  estimate.pmufseq = estimate.pmufseq/sumpmufseq;
