function xiNext = myphi_GPB2(xi,u,y,model)
  %UNTITLED Summary of this function goes here
  %   Detailed explanation goes here

  % Construct filtering estimate
  estimate.xfseq = [xi(1) xi(2)];
  estimate.Pxxfseq(:,:,1) = xi(3);
  estimate.Pxxfseq(:,:,2) = xi(4);
  estimate.pmufseq = [xi(5);1-xi(5)];

  % Do prediction
  estimate = smmkfp(u,estimate,model);

  % Do filtering
  xiNext = nan(5,size(y,2));
  for i = 1:size(y,2)
    estimate = smmkff(y(:,i),estimate,model);
    % Do merging
    [pmufseq,xfseq,Pxxfseq] = mmmerge(estimate.pmufseq,estimate.xfseq,estimate.Pxxfseq,2);
    % Construct state at next time step
    xiNext(:,i) = [xfseq(1) xfseq(2) Pxxfseq(:,:,1) Pxxfseq(:,:,2) pmufseq(1)]';
  end
end
