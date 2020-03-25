function [xiNext,varargout] = myphi_GPB2_1D(xi,u,y,model)
  %UNTITLED Summary of this function goes here
  %   Detailed explanation goes here
  % xi(1) - mean of the 1st model, x(2) - mean of the 2nd model
  % xi(3) - var of the 1st model, x(4) - var of the 2nd model
  % xi(5) - probability of the 1st model

  ly  = size(y,2);
  % mu(k+1)=1 & mu(k)=1, mu(k+1)=2 & mu(k)=1, mu(k+1)=1 & mu(k)=2, mu(k+1)=2 & mu(k)=2  
  merge = [1 0 1 0;...
           0 1 0 1]; % which models should be merged  (kron(ones(1,2),eye(2))
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  PREDICTION: (nM means, vars and weights) -> (nM means and vars, nM^2 weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiPred_means = model.As'.*xi(1:2) + model.Bs'*u;
  xiPred_means_replicated = xiPred_means([1,1,2,2]);
  xiPred_vars = model.As'.^2.*xi(3:4) + model.Gs'.^2;
  %xiPred_weights = [model.P(:,1)*xi(5);model.P(:,2)*(1-xi(5))];
  A =  model.P * [1 0;0 -1];
  xiPred_weights = A(:)*xi(5)+[0 0 model.P(1,2) model.P(2,2)]';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FILTERING%  (nM means and vars, nM^2 weights) -> (nM^2 means, vars and weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ybars: nM^2 x length(y)
  %ybars = y(ones(4,1),:)- model.Css'.*xiPred_means([1,1,2,2]);
  ybars = y- model.Css'.*xiPred_means_replicated;
  % Pyys: nM^2 x 1
  Pyys = model.Css'.^2.*xiPred_vars([1 1 2 2]) + model.Hss'.^2;
  % KGs: nM^2 x 1
  KGs = model.Css'.*xiPred_vars([1 1 2 2])./Pyys;
  % weights: nM^2 x 1 ! WORKS ONLY IN NEW MATLAB 2019 !
  %weights = 1./sqrt(2*pi*Pyys).*exp(-0.5*ybars.^2./Pyys(:,ones(1,ly)));
  weights = 1./sqrt(2*pi*Pyys).*exp(-0.5*ybars.^2./Pyys);
  % xiFilt_means: nM^2 x length(y) (after the addition) 
  %xiFilt_means = (ones(4,1)-KGs.*model.Css').*(xiPred_means([1 1 2 2])) + KGs*y; 
  xiFilt_means = xiPred_means_replicated + KGs.*(y - model.Css'.*xiPred_means_replicated); 
  % xiFilt_vars: nM^2 x 1
  xiFilt_vars = (model.Hss'.^2./Pyys).*(xiPred_vars([1 1 2 2]));
  % xiFilt_weights: nM^2 x length(y)
  %xiFilt_weights = weights .* xiPred_weights(:,ones(1,ly));
  xiFilt_weights = weights .* xiPred_weights;
  %xiFilt_weights(xiFilt_weights<=eps(max(xiFilt_weights))) = 0
  xiFilt_weights = xiFilt_weights./sum(xiFilt_weights);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % MERGING (nM^2 means, vars and weights) -> (nM means, vars and weights)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % The following five lines commented and modified by Ivo
  %xiMerged_means = merge*(xiFilt_means.*xiFilt_weights);
  %xiMerged_vars = merge*((xiFilt_vars + (xiFilt_means-xiMerged_means([1 1 2 2],:)).^2).*xiFilt_weights);
  %xiMerged_vars = xiMerged_vars./(merge*xiFilt_weights);
  %xiMerged_vars(isnan(xiMerged_vars))=1;
  %xiMerged_weights = [1 0]*merge*xiFilt_weights;
  % Modification by Ivo
  xiMerged_weights = merge*xiFilt_weights;
  normalizer = 1./xiMerged_weights;
  normalizer(isnan(normalizer)) = 1;
  xiMerged_means = merge*(xiFilt_means.*xiFilt_weights);
  xiMerged_means = xiMerged_means.*normalizer;
  xiMerged_vars = merge*((xiFilt_vars + (xiFilt_means-xiMerged_means([1 2 1 2],:)).^2).*xiFilt_weights);
  xiMerged_vars = xiMerged_vars.*normalizer;
  xiMerged_weights = [1 0]*xiMerged_weights;
  % End of modification by Ivo
  %xiMerged = [xiMerged_means; xiMerged_vars; xiMerged_weights; xi(6:7,ones(1,ly))];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PASSING OUTPUT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  xiNext = zeros(5,ly);
  xiNext(1:2,:) = xiMerged_means;
  xiNext(3:4,:) = xiMerged_vars;
  xiNext(5,:) = xiMerged_weights;
  if nargout == 3
    varargout{1} = weights;
    tmp.xiFilt_means = xiFilt_means;
    tmp.xiFilt_vars = xiFilt_vars;
    tmp.xiFilt_weights = xiFilt_weights;
    varargout{2} = tmp;
  end
end
