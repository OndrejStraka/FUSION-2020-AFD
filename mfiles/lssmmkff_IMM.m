function estimate = lssmmkff_IMM(y,estimate,model,nodeStates,approxMETHOD)
%UNTITLED2 Large scale sitching multiple model Kalamn filter - filtering step
%   Detailed explanation goes here

% Get the number of nodes
nNode = length(y);

if nNode>2
    error('The estimation algorithm is implemented for two nodes only')
end

% Get the dimension of the state
nx = size(nodeStates,1);

% Perform filtering step for each node individualy
for iNode = 1:nNode
  if size(estimate{iNode}.xpseq,2) == 1 % first time step - only one estimate -> no IMM
    estimate{iNode} = smmkff_GPB2(y{iNode},estimate{iNode},model(iNode));
  else
    estimate{iNode} = smmkff_IMM(y{iNode},estimate{iNode},model(iNode));
  end
end

% Nodes exchange their local estimates and each computes the same set of
% means and convariance matrices using covariance intersection algorithm
% The code is written for two node only

% Get the number of filtering estimates for individual nodes
nModel = [size(estimate{1}.xfseq,2);size(estimate{2}.xfseq,2)];

%Pci = zeros(nx,nx,nModel(1),nModel(2));
%xci = zeros(nx,nModel(1),nModel(2));
% Compute covariance intersection for every model combination
for i = 1:nModel(1)
    for j = 1:nModel(2)
        % Create state
        xci(nodeStates(:,1),i,j) = estimate{1}.xfseq(:,i);
        xci(nodeStates(:,2),i,j) = estimate{2}.xfseq(:,j);
        
        Pci(nodeStates(:,1),nodeStates(:,1),i,j) = estimate{1}.Pxxfseq(:,:,i);
        Pci(nodeStates(:,2),nodeStates(:,2),i,j) = estimate{2}.Pxxfseq(:,:,j);
    end
end

% Multiply the covariance matrices to account for unknown correlation
Pci = 2*Pci;

% Gaussian sum to single gaussian for node 1
estimate{1}.xfseqfused = zeros(nx,nModel(1));
estimate{1}.Pxxfseqfused = zeros(nx,nx,nModel(1));
for i = 1:nModel(1)
    [estimate{1}.xfseqfused(:,i),estimate{1}.Pxxfseqfused(:,:,i)] = fm2mcm(estimate{2}.pmufseq',squeeze(xci(:,i,:)),squeeze(Pci(:,:,i,:)));
end

% Gaussian sum to single gaussian for node 2
estimate{2}.xfseqfused = zeros(nx,nModel(2));
estimate{2}.Pxxfseqfused = zeros(nx,nx,nModel(2));
for i = 1:nModel(2)
    [estimate{2}.xfseqfused(:,i),estimate{2}.Pxxfseqfused(:,:,i)] = fm2mcm(estimate{1}.pmufseq',squeeze(xci(:,:,i)),squeeze(Pci(:,:,:,i)));
end



% for iNode = 1:2
%     % Replace the local filtering estimates by the fused estimates
%     estimate{iNode}.xfseq = estimate{iNode}.xfseqfused(nodeStates(:,iNode),:);
%     estimate{iNode}.Pxxfseq = estimate{iNode}.Pxxfseqfused(nodeStates(:,iNode),nodeStates(:,iNode),:);
% 
% end
% 
