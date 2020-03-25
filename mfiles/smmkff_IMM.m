function estimate = smmkff_IMM(y,estimate,model,computeOnlyPredictiveOutput)
  % SMMKFF filtering step of Switching Multiple Model estimation
  % algorithm with IMM
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
  nXf = nXp;

  % Set filtering values to zeros
  estimate.xfseq = zeros(nx,nXf);
  estimate.Pxxfseq = zeros(nx,nx,nXf);
  estimate.pmufseq = zeros(nXf,1);

  % Set predictive output values to zeros
  estimate.ypseq = zeros(ny,nXf);
  estimate.Pyypseq = zeros(ny,ny,nXf);

  %myp = zeros(nXf,1);
  %srdetPyp = zeros(nXf,1);

  for i = 1:nXp
    % INDICES of probabilities for the i-th mode P(mu(k+1)=*,mu(k)=i|I)
    % idxs = [(i-1)*nModels+1:i*nModels];
    % CALCULATING NORMALIZATION P(mu(k+1)|I^k)
    % INDICES of probabilities for the i-th mode P(mu(k+1)=i,mu(k)=*|I)
    idxs_new = [0:nModels:nModels^2-1]+i;
    % MIXING PROBABILITIES for i-th mode
    pmumseqNORM = sum(estimate.pmupseq(idxs_new)); % sum of predictive probabilities for i-th model 
    if pmumseqNORM == 0
      pmumseqTMP = zeros(nModels,1);
    else
      pmumseqTMP = estimate.pmupseq(idxs_new) / pmumseqNORM; % THIS GIVES TRANSIENT PROBABILTIES OF THE MODEL
    end

    % MIXING MODES
    [xmseq(:,i),Pxxmseq(:,:,i)] = fm2mcm(pmumseqTMP',estimate.xpseq,estimate.Pxxpseq);

    % FILTERING
    [estimate.xfseq(:,i),estimate.Pxxfseq(:,:,i),...
      estimate.ypseq(:,i),estimate.Pyypseq(:,:,i)] = ...
      kff(y,xmseq(:,i),Pxxmseq(:,:,i),model.M(i).C,model.M(i).H*model.M(i).H');
    % Compute non-normalized filtering probability of sequence
    % it's not numerically robust :(
    estimate.pmufseq(i) = normpdfm(y,estimate.ypseq(:,i),estimate.Pyypseq(:,:,i))*pmumseqNORM;

    % Prepare some values used in numerically robust computation of
    % filtering sequence probabilities
    %e = y-estimate.ypseq(:,idx);
    %myp(idx) = e'*(estimate.Pyypseq(:,:,idx)\e);
    %srdetPyp(idx) = sqrt(det(estimate.Pyypseq(:,:,idx))); 
  end

  if computeOnlyPredictiveOutput
    return
  end

  % Perform weights normalization to obtain probabilities
  % The weights are normalized to sum up 1
  sumpmufseq = sum(estimate.pmufseq);

  if sumpmufseq == 0
    disp('Sequence probabilities were computed in an alternative way.')
    keyboard
  end

  % The weights that are less or equal to EPS(max_weight) are set to zero
  estimate.pmufseq(estimate.pmufseq<=eps(max(estimate.pmufseq))) = 0;

  % Normalize
  sumpmufseq = sum(estimate.pmufseq);
  if sumpmufseq == 0
    estimate.pmufseq = zeros(size(estimate.pmufseq));
  else
    estimate.pmufseq = estimate.pmufseq/sumpmufseq;
  end
