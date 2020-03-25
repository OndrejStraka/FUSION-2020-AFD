function xiNext = myphi_IMM_1D(xi,u,y,model)
  %UNTITLED Summary of this function goes here
  %   Detailed explanation goes here
  % xi(1) - mean of the 1st model, x(2) - mean of the 2nd model
  % xi(3) - var of the 1st model, x(4) - var of the 2nd model
  % xi(5) - probability of the 1st model

  ly  = size(y,2);
  % mu(k+1)=1 & mu(k)=1, mu(k+1)=2 & mu(k)=1, mu(k+1)=1 & mu(k)=2, mu(k+1)=2 & mu(k)=2  
  merge = [1 0 1 0;...
           0 1 0 1]; % which models should be merged 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PREDICTION: (nM means, vars and weights) -> (nM means and vars, nM^2 weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiPred_means = model.As'.*xi(1:2) + model.Bs'*u;
  xiPred_vars = model.As'.^2.*xi(3:4) + model.Gs'.^2;
  %xiPred_weights = [model.P(:,1)*xi(5);model.P(:,2)*(1-xi(5))];
  A =  model.P * [1 0;0 -1];
  xiPred_weights = A(:)*xi(5)+[0 0 model.P(1,2) model.P(2,2)]';
  %[mu(k+1)=1&mu(k)=1 mu(k+1)=2&mu(k)=1 mu(k+1)=1&mu(k)=2 mu(k+1)=2&mu(k)=2] 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MIXING:  (nM means and vars, nM^2 weights) -> (nM means, vars and weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INDICES of probabilities for the i-th mode P(mu(k+1)=i,mu(k)=*|I)
  idxs = [1 2;3 4];
  xiMix_probabilities = xiPred_weights(idxs); % columns for mu(k+1)=i
  %[mu(k+1)=1&mu(k)=1 mu(k+1)=2&mu(k)=1;mu(k+1)=1&mu(k)=2 mu(k+1)=2&mu(k)=2] 
  % should be same as (model.P*([1;1]*[xi(5) 1-xi(5)]))'
  xiMix_weights = sum(xiMix_probabilities)'; % [mu(k+1)=1;mu(k+1)=2] (sums of each column)
  normalizer = 1./xiMix_weights'; %[1/mu(k+1)=1 1/mu(k+1)=2]
  normalizer(isinf(normalizer)) = 0; 
  xiMix_probabilities = xiMix_probabilities.*normalizer; % normalized mixing probabilities (columns)
  %[mu(k+1)=1&mu(k)=1/mu(k+1)=1 mu(k+1)=2&mu(k)=1/mu(k+1)=2;
  % mu(k+1)=1&mu(k)=2/mu(k+1)=1 mu(k+1)=2&mu(k)=2/mu(k+1)=2] 
  xiMix_means = xiMix_probabilities'*xiPred_means;
  %[mu(k+1)=1&mu(k)=1/mu(k+1)=1 mu(k+1)=1&mu(k)=2/mu(k+1)=1;
  % mu(k+1)=2&mu(k)=1/mu(k+1)=2 mu(k+1)=2&mu(k)=2/mu(k+1)=2]*[m(1);
  %                                                           m(2)] 
  xiMix_vars = [sum(xiMix_probabilities(:,1).*(xiPred_vars + (xiPred_means - xiMix_means(1)).^2));... 
                sum(xiMix_probabilities(:,2).*(xiPred_vars + (xiPred_means - xiMix_means(2)).^2))];
  %sum([mu(k+1)=1&mu(k)=1/mu(k+1)=1;
  %     mu(k+1)=1&mu(k)=2/mu(k+1)=1].*[v(1);
  %                                    v(2)]) - xiMixmeans(1)^2
  %sum([mu(k+1)=2&mu(k)=1/mu(k+1)=2;
  %     mu(k+1)=2&mu(k)=2/mu(k+1)=2].*[v(1);
  %                                    v(2)]) - xiMixmeans(2)^2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FILTERING%  (nM means and vars and weights) -> (nM means, vars and weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ybars: nM x length(y)
  ybars = y- model.Cs'.*xiMix_means;
  % Pyys: nM x 1
  Pyys = model.Cs'.^2.*xiMix_vars + model.Hs'.^2;
  % KGs: nM x 1
  KGs = model.Cs'.*xiMix_vars./Pyys;
  % weights: nM x 1 ! WORKS ONLY IN NEW MATLAB 2019 !
  weights = 1./sqrt(2*pi*Pyys).*exp(-0.5*ybars.^2./Pyys);
  % xiFilt_means: nM x length(y) (after the addition) 
  xiFilt_means = xiMix_means + KGs.*(y - model.Cs'.*xiMix_means); 
  % xiFilt_vars: nM x 1
  xiFilt_vars = (model.Hs'.^2./Pyys).*xiMix_vars;
  % xiFilt_weights: nM x length(y)
  xiFilt_weights = weights .* xiMix_weights;
  %xiFilt_weights(xiFilt_weights<=eps(max(xiFilt_weights))) = 0
  xiFilt_weights = xiFilt_weights./sum(xiFilt_weights);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PASSING OUTPUT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiNext = zeros(5,ly);
  xiNext(1:2,:) = xiFilt_means;
  xiNext(3:4,:) = xiFilt_vars(:,ones(1,ly));
  xiNext(5,:) = xiFilt_weights(1,:);
end
