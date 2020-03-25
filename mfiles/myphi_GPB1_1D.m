function xiNext = myphi_GPB1_1D(xi,u,y,model)
  %UNTITLED evaluate phi for multiple measurements y
  %   Detailed explanation goes here
  ly  = size(y,2);
  % xi(1) - mean of the 1st model, x(2) - mean of the 2nd model
  % xi(3) - var of the 1st model,  x(4) - var of the 2nd model
  % xi(5) - probability of the 1st model
  xiFilt_mean = xi(1);
  xiFilt_var = xi(2);
  xiFilt_weights = [xi(3);1-xi(3)];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PREDICTION  (1 mean, var and weights) -> (nM means, vars and weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiPred_means = model.As'.*xiFilt_mean + model.Bs'*u;
  xiPred_vars = model.As'.^2.*xiFilt_var + model.Gs'.^2;
  xiPred_weights = model.P*xiFilt_weights;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MERGING To SINGLE TERM (nM means and vars) -> (1 mean and var)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiPredMerged_mean = sum(xiPred_means.*xiFilt_weights);
  xiPredMerged_var = sum((xiPred_vars + (xiPred_means-xiPredMerged_mean([1 1],:)).^2).*xiFilt_weights);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FILTERING (1 mean and var) -> (nM means and vars)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ybars: nM x length(y)
  ybars = y- model.Cs'.*xiPredMerged_mean;
  % Pyys: nM x 1
  Pyys = model.Cs'.^2.*xiPredMerged_var + model.Hs'.^2;
  % KGs: nM x 1
  KGs = model.Cs'.*xiPredMerged_var./Pyys;
  % weights: nM x 1 ! WORKS ONLY IN NEW MATLAB 2019 !
  %weights = 1./sqrt(2*pi*Pyys).*exp(-0.5*ybars.^2./Pyys(:,ones(1,ly)));
  weights = 1./sqrt(2*pi*Pyys).*exp(-0.5*ybars.^2./Pyys);
  % xiFilt_means: nM x length(y) (after the addition) 
  xiFilt_means = xiPredMerged_mean + KGs.*ybars; 
  % xiFilt_vars: nM x 1
  xiFilt_vars = (model.Hs'.^2./Pyys).*xiPredMerged_var;
  % xiFilt_weights: nM x length(y)
  xiFilt_weights = weights .* xiPred_weights;
  %xiFilt_weights(xiFilt_weights<=eps(max(xiFilt_weights))) = 0
  xiFilt_weights = xiFilt_weights./sum(xiFilt_weights);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MERGING To SINGLE TERM (nM means, vars and weights) -> (1 mean, var and weight=1)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiFiltMerged_mean = sum(xiFilt_means.*xiFilt_weights);
  xiFiltMerged_var = sum((xiFilt_vars + (xiFilt_means-xiFiltMerged_mean([1 1],:)).^2).*xiFilt_weights);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PASSING OUTPUT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiNext = zeros(3,ly);
  xiNext(1,:) = xiFiltMerged_mean;
  xiNext(2,:) = xiFiltMerged_var(:,ones(1,ly)); % same variances for all estimates
  xiNext(3,:) = xiFilt_weights(1,:);
end
