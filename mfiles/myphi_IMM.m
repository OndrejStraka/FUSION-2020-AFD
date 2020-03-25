function xiNext = myphi_IMM(xi,u,y,model)
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
  for i = 1:size(y,2) % FOR EACH  NODE
    estimate = smmkff_IMM(y(:,i),estimate,model);
    % Construct state at next time step
    xiNext(:,i) = [estimate.xfseq(1) estimate.xfseq(2) estimate.Pxxfseq(:,:,1) estimate.Pxxfseq(:,:,2) estimate.pmufseq(1)]';
  end
end
