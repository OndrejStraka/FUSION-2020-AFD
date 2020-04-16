%function largescalesystemfd_simple(designMergeMethod,mergeMethod)
close all 
clear all
clc

localP = 'Bayes';
%localP = 'mincon';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LSS
%LSS.P = [100 5   5   1;
         %5   100 1   5;
         %5   1   100 5;
         %1   5   5   100];
LSS.P = [100 5   5   1;
         2   100 1   2;
         2   1   100 2;
         1   20  20  100];
LSS.P = LSS.P./sum(LSS.P,1);

% computation of local transition probabilities by Bayes
switch localP
  case 'Bayes'
    [tmpV,tmpD] = eig(LSS.P);
    LSS.Pi = tmpV(:,abs(diag(tmpD)-1)<eps);
    LSS.Pi = LSS.Pi/sum(LSS.Pi);   
    n1_P_Bayes = kron(eye(2),ones(1,2))*(LSS.P.*LSS.Pi')*kron(eye(2),ones(2,1));
    n1_P_Bayes = n1_P_Bayes./sum(n1_P_Bayes,1);
    n2_P_Bayes = kron(ones(1,2),eye(2))*(LSS.P.*LSS.Pi')*kron(ones(2,1),eye(2));
    n2_P_Bayes = n2_P_Bayes./sum(n2_P_Bayes,1);
  case 'mincon'
    % computation of local transition probabilities by constrained minimization
    P11 = LSS.P(1:2,1:2);
    P12 = LSS.P(1:2,3:4);
    P21 = LSS.P(3:4,1:2);
    P22 = LSS.P(3:4,3:4);
    permP = [P11(:)';P21(:)';P12(:)';P22(:)'];
    objective = @(x) norm(permP - [x(1);1-x(1);x(2);1-x(2)]*[x(3) 1-x(3) x(4) 1-x(4)],'fro');
    x0 = [0.5 0.5 0.5 0.5];
    lb = [0 0 0 0];
    ub = [1 1 1 1];
    xmin = fmincon(objective,x0,[],[],[],[],lb,ub);
    n1_P_mincon = [xmin(1) xmin(2);1-xmin(1) 1-xmin(2)];
    n2_P_mincon = [xmin(3) xmin(4);1-xmin(3) 1-xmin(4)];
  otherwise
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NODE 1
% Node 1 - Model 1
node(1).M(1).A = [0.76 0.05];
node(1).M(1).B = 0.12;
node(1).M(1).G = sqrt(0.003);
node(1).M(1).C = 0.9;
node(1).M(1).H = 0.01;
% Node 1 - Model 2
node(1).M(2).A = [0.86 0.15];
node(1).M(2).B = 0.14;
node(1).M(2).G = sqrt(0.003);
node(1).M(2).C = 1;
node(1).M(2).H = 0.01;
% Node 1 - Transition probabilities
switch localP
  case 'Bayes'
    node(1).P = n1_P_Bayes;
  case 'mincon'
    node(1).P = n1_P_mincon;
  otherwise
end

    
% Node 1 - Probabilities of initial model
node(1).Pmua = [1;0];
% Node 1 - Mean and variance of the node's initial state
node(1).xa = 0;
node(1).Pxxa = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NODE 2
% Node 2 - Model 1
node(2).M(1).A = [0.10 0.87];
node(2).M(1).B = 0.13;
node(2).M(1).G = sqrt(0.002);
node(2).M(1).C = 0.9;
node(2).M(1).H = 0.01;
% Node 2 - Model 2
node(2).M(2).A = [0.05 0.775];
node(2).M(2).B = 0.15;
node(2).M(2).G = sqrt(0.002);
node(2).M(2).C = 1;
node(2).M(2).H = 0.01;
% Node 2 - Transition probabilities
switch localP
  case 'Bayes'
    node(2).P = n2_P_Bayes;
  case 'mincon'
    node(2).P = n2_P_mincon;
  otherwise
end
% Node 2 - Probabilities of initial model
node(2).Pmua = [1;0];
% Node 2 - Mean and variance of the node's initial state
node(2).xa = 0;
node(2).Pxxa = 0.01;

LSS.Pmua = kron(node(1).Pmua,node(2).Pmua);
node(1).mu_from_LSS = kron([1;2],ones(2,1)); % [1 1 2 2]'
node(2).mu_from_LSS = kron(ones(2,1),[1;2]); % [1 2 1 2]'

% The columns define states that belongs to individual nodes
nodeStates = logical([1 0;0 1]);

%LSS.P - kron(node(1).P,node(2).P) % DBG
%LSS.Pi
%node(1).P
%node(2).P
%pause

nNode = length(node);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AFD SETTINGS
% The discount factor
%lambda = 0.98;
lambda = 0.9;

%architecture = 'decentralized';
architecture = 'distributed';

%mergeMethod = 'GPB1';
%mergeMethod = 'IMM';
mergeMethod = 'GPB2';
%mergeMethod = 'DEBUG'

activeMode = true;
design = false;
%designMergeMethod = 'GPB1';
%designMergeMethod = 'IMM';
designMergeMethod = 'GPB2';
verbose = true;
debug = false;
maxIteration = 70;

% CENTRAL NODE SETTINGS
centralDecision = true;
returnPmuToLocalnodes = true;
LSS.returnPmuToLocalnodes = returnPmuToLocalnodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION SETTINGS
% Time horizon
F = 400;
Fp1 = F + 1;
t = 0:F;
% Criterion
Ldmax = 2;
Jtail = lambda^(F+1)*Ldmax/(1-lambda);
% #Monte Carlo simulations
nmc = 1e5;
% Random number generation
rng('default');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OFF-LINE DESIGN
if design
  % The grid of states
  xrange = [-1.5:0.1:1.5];
  Pxrange = [1.9e-4:0.05e-4:2e-4];
  Prange = [0:0.01:1];

  for iSubsys = 1:2 % DEFINE VARIABLES
    % Define decentralized models for node iNode
    nodedec{iSubsys} = node(iSubsys);
    nodedec{iSubsys}.M(1).A = nodedec{iSubsys}.M(1).A(1,iSubsys);
    nodedec{iSubsys}.M(2).A = nodedec{iSubsys}.M(2).A(1,iSubsys);
    nodedec{iSubsys}.As = [nodedec{iSubsys}.M(1).A nodedec{iSubsys}.M(2).A];
    nodedec{iSubsys}.Bs = [nodedec{iSubsys}.M(1).B nodedec{iSubsys}.M(2).B];
    nodedec{iSubsys}.Gs = [nodedec{iSubsys}.M(1).G nodedec{iSubsys}.M(2).G];
    nodedec{iSubsys}.Cs = [nodedec{iSubsys}.M(1).C nodedec{iSubsys}.M(2).C];
    nodedec{iSubsys}.Css = [nodedec{iSubsys}.M(1).C nodedec{iSubsys}.M(2).C nodedec{iSubsys}.M(1).C nodedec{iSubsys}.M(2).C];
    nodedec{iSubsys}.Hs = [nodedec{iSubsys}.M(1).H nodedec{iSubsys}.M(2).H];
    nodedec{iSubsys}.Hss = [nodedec{iSubsys}.M(1).H nodedec{iSubsys}.M(2).H nodedec{iSubsys}.M(1).H nodedec{iSubsys}.M(2).H];
    
    inputs{iSubsys} = {[-1,0,1]};
    switch designMergeMethod
      case 'GPB1'
        states{iSubsys} = {xrange,Pxrange,Prange};
        mL_h{iSubsys} = @(xi,u) min(xi(3),1-xi(3));
      otherwise
        states{iSubsys} = {xrange,xrange,Pxrange,Pxrange,Prange};
        mL_h{iSubsys} = @(xi,u) min(xi(5),1-xi(5));
    end
    nS{iSubsys} = cellfun(@length,states{iSubsys});
    innerEdges{iSubsys} = generateinneredges(states{iSubsys});

      isStateAdmissible_h{iSubsys} = @(xi) true(1,size(xi,2));
      isInputAdmissible_h{iSubsys} = @(u) true(1,size(u,2));

    switch designMergeMethod
      case'GPB1'
        transition_h{iSubsys} = @(xi,u,y) myphi_GPB1_1D(xi,u,y,nodedec{iSubsys});
        disturbance_gs_h{iSubsys} = @(xi,u) mynoise_GPB1_1D(xi,u,nodedec{iSubsys});
      case'IMM'
        transition_h{iSubsys} = @(xi,u,y) myphi_IMM_1D(xi,u,y,nodedec{iSubsys});
        disturbance_gs_h{iSubsys} = @(xi,u) mynoise_1D(xi,u,nodedec{iSubsys});
      case'GPB2'
        transition_h{iSubsys} = @(xi,u,y) myphi_GPB2_1D(xi,u,y,nodedec{iSubsys});
        disturbance_gs_h{iSubsys} = @(xi,u) mynoise_1D(xi,u,nodedec{iSubsys});
      otherwise
        error('unknown merge Method')
    end

    mV_h{iSubsys} = @(innerEdges,V,xi,u) mVut_1D(innerEdges,V,xi,u,transition_h{iSubsys},disturbance_gs_h{iSubsys});
  end
  for iSubsys = 1:2 % VALUE ITERATION
    [gamma{iSubsys},Q{iSubsys},V{iSubsys},iterationConvergence{iSubsys},iterationTimes{iSubsys}] = valueiteration_discraggr2(states{iSubsys},innerEdges{iSubsys},inputs{iSubsys},mL_h{iSubsys},mV_h{iSubsys},isStateAdmissible_h{iSubsys},isInputAdmissible_h{iSubsys},lambda,[],maxIteration,[],'non-vectorized',verbose);
  end
  fileName = sprintf('decentralized-design-%s',localP);
  save(fileName,'nodedec','states','nS','innerEdges','inputs',...
    'gamma','Q','V','iterationConvergence','iterationTimes')
end

% Load decentralized controllers
fileName = sprintf('decentralized-design-%s',localP);
load(fileName,'nS','innerEdges','gamma','nodedec');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ON-LINE PART
% Initialize the predictive estimates at the first time step and each node
k = 1;
for iNode = 1:nNode
  estimate{iNode,k}.xpseq = node(iNode).xa;
  estimate{iNode,k}.Pxxpseq = node(iNode).Pxxa;
  estimate{iNode,k}.pmupseq = node(iNode).Pmua;
end

x = zeros(nNode,Fp1);
mu = zeros(nNode,Fp1);
LSS.mu = zeros(1,Fp1);

tic
for imc = 1:nmc
  fprintf('A:<strong>%s</strong>,GlobalDecision:<strong>%s</strong>,Feedback:<strong>%s</strong>\t MC simulation #:<strong>%i</strong>\n',architecture,mat2str(centralDecision),mat2str(returnPmuToLocalnodes),imc)
  k = 1;
  LSS.mu(:,k) = gendrnd(LSS.Pmua); % LSS models
  for iNode = 1:nNode % 
    x(nodeStates(:,iNode),1) = normrndm(node(iNode).xa,node(iNode).Pxxa); % initial state
    mu(iNode,k) = node(iNode).mu_from_LSS(LSS.mu(:,k)); % local nodes
    nv(iNode,1) = size(node(iNode).M(1).H,2); % measurement noise dimension
    v{iNode} = randn(nv(iNode),Fp1); % measurement noise
    nw(iNode,k) = size(node(iNode).M(1).G,2); % state noise dimension
    w{iNode} = randn(nw(iNode),F); % state noise
  end

  % prepare variables
  y = cell(nNode,Fp1);
  u = cell(nNode,Fp1); 
  d = nan(nNode,Fp1);
  
  for k = 1:Fp1
    if k>1 % Dynamics in k=2,3,...
      for iNode = 1:nNode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct information state of the i-th node 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch architecture
          case 'decentralized'
            switch mergeMethod
              case 'GPB1' % 3D -> 3D or 5D
                xi{iNode} = [estimate{iNode,k-1}.xfseq;estimate{iNode,k-1}.Pxxfseq(1,1,1); estimate{iNode,k-1}.pmufseq(1)];
                switch designMergeMethod
                  case 'GPB1' %3D->3D Bellman
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseq;estimate{iNode,k-1}.Pxxfseq(1,1,1); estimate{iNode,k-1}.pmufseq(1)];
                  otherwise %3D->5D Bellman
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseq;estimate{iNode,k-1}.xfseq;estimate{iNode,k-1}.Pxxfseq(1,1,1);estimate{iNode,k-1}.Pxxfseq(1,1,1); estimate{iNode,k-1}.pmufseq(1)];
              end
              otherwise % 5D -> 3D or 5D
                xi{iNode} = [estimate{iNode,k-1}.xfseq';estimate{iNode,k-1}.Pxxfseq(1,1,1); estimate{iNode,k-1}.Pxxfseq(1,1,2); estimate{iNode,k-1}.pmufseq(1)];
                switch designMergeMethod
                  case 'GPB1' %5D->3D Bellman
                    tmp_m = sum(estimate{iNode,k-1}.xfseq.*estimate{iNode,k-1}.pmufseq');
                    tmp_v = sum((squeeze(estimate{iNode,k-1}.Pxxfseq)' + (estimate{iNode,k-1}.xfseq-tmp_m).^2).*estimate{iNode,k-1}.pmufseq');
                    xi_u{iNode} = [tmp_m;tmp_v; estimate{iNode,k-1}.pmufseq(1)];
                  otherwise %5D->5D
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseq';estimate{iNode,k-1}.Pxxfseq(1,1,1); estimate{iNode,k-1}.Pxxfseq(1,1,2); estimate{iNode,k-1}.pmufseq(1)];
                end
            end
          case 'distributed'
            switch mergeMethod
              case 'GPB1' % 3D -> 3D or 5D
                switch designMergeMethod
                  case 'GPB1' %3D->3D Bellman
                    flagS = nodeStates(:,iNode);
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseqfused(flagS,:);estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,1); estimate{iNode,k-1}.pmufseq(1)];
                  otherwise %3D->5D Bellman
                    flagS = nodeStates(:,iNode);
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseqfused(flagS,:);estimate{iNode,k-1}.xfseqfused(flagS,:);estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,1);estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,1);estimate{iNode,k-1}.pmufseq(1)];
                end
              otherwise % 5D -> 3D or 5D
                switch designMergeMethod
                  case 'GPB1' %5D->3D Bellman
                    flagS = nodeStates(:,iNode);
                    tmp_m = sum(estimate{iNode,k-1}.xfseqfused(flagS,:).*estimate{iNode,k-1}.pmufseq');
                    tmp_v = (squeeze(estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,1)) + (estimate{iNode,k-1}.xfseqfused(flagS,1)-tmp_m).^2)*estimate{iNode,k-1}.pmufseq(1)...
                           +(squeeze(estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,2)) + (estimate{iNode,k-1}.xfseqfused(flagS,2)-tmp_m).^2)*estimate{iNode,k-1}.pmufseq(2);
                    xi_u{iNode} = [tmp_m;tmp_v; estimate{iNode,k-1}.pmufseq(1)];
                  otherwise %5D->5D
                    flagS = nodeStates(:,iNode);
                    xi_u{iNode} = [estimate{iNode,k-1}.xfseqfused(flagS,:)';estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,1); estimate{iNode,k-1}.Pxxfseqfused(flagS,flagS,2); estimate{iNode,k-1}.pmufseq(1)];
                end
            end
          otherwise
            error('unknown architecture')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute input generated by the i-th node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if activeMode
          u{iNode,k-1} =  gamma{iNode}(subidx2linidx(nS{iNode},aggregationidx(xi_u{iNode},innerEdges{iNode})));
        else % passive FD
          u{iNode,k-1} = sign(sin(0.2*(k-1)));
        end
        % State dynamics of the i-th node
        x(nodeStates(:,iNode),k) = node(iNode).M(mu(iNode,k-1)).A*x(:,k-1) +...
          node(iNode).M(mu(iNode,k-1)).B*u{iNode,k-1} + node(iNode).M(mu(iNode,k-1)).G*w{iNode}(:,k-1);
      end % input generation
      % Generate models in the next time instant
      LSS.mu(:,k) = gendrnd(LSS.P(:,LSS.mu(:,k-1))); % generate LSS models
      for iNode = 1:nNode
        mu(iNode,k) = node(iNode).mu_from_LSS(LSS.mu(:,k)); % generate local models
      end
    end % k>1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate measurements for individual subsystems
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iNode = 1:nNode
      y{iNode,k} = node(iNode).M(mu(iNode,k)).C*x(nodeStates(:,iNode),k) +...
        node(iNode).M(mu(iNode,k)).H*v{iNode}(:,k);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch architecture
      case 'decentralized'
        if k == 1
          % Estimation for k = 1 is SAME for all mergeMethods since we start from a
          % single model => no merging
          for iNode = 1:nNode
            switch mergeMethod
              %case 'GPB1'
                %estimate{iNode,k} = smmkff_GPB1(y{iNode,k},estimate{iNode,k},nodedec{iNode});
              %case 'IMM'
                %estimate{iNode,k} = smmkff_GPB2(y{iNode,k},estimate{iNode,k},nodedec{iNode});
              case 'GPB2'
                estimate{iNode,k} = smmkff_GPB2(y{iNode,k},estimate{iNode,k},nodedec{iNode});
              otherwise
                error('unknown merge method')
            end
          end
        else % Estimation for k > 1
          for iNode = 1:nNode
            switch mergeMethod
              %case 'GPB1'
                %[xiNext] = myphi_GPB1_1D(xi{iNode},u{iNode,k-1},y{iNode,k},nodedec{iNode});
                %estimate{iNode,k}.xfseq = xiNext(1)';
                %estimate{iNode,k}.Pxxfseq(1,1,1) = xiNext(2);
                %estimate{iNode,k}.pmufseq = [xiNext(3) 1-xiNext(3)]';
              %case 'IMM'
                %[xiNext] = myphi_IMM_1D(xi{iNode},u{iNode,k-1},y{iNode,k},nodedec{iNode});
                %estimate{iNode,k}.xfseq = xiNext(1:2)';
                %estimate{iNode,k}.Pxxfseq(1,1,1) = xiNext(3);
                %estimate{iNode,k}.Pxxfseq(1,1,2) = xiNext(4);
                %estimate{iNode,k}.pmufseq = [xiNext(5) 1-xiNext(5)]';
              case 'GPB2'
                [xiNext,estimate{iNode,k}.likelihood,tmp{iNode}] = myphi_GPB2_1D(xi{iNode},u{iNode,k-1},y{iNode,k},nodedec{iNode});
                if (~centralDecision) && (~returnPmuToLocalnodes)
                  estimate{iNode,k}.xfseq = xiNext(1:2)';
                  estimate{iNode,k}.Pxxfseq(1,1,1) = xiNext(3);
                  estimate{iNode,k}.Pxxfseq(1,1,2) = xiNext(4);
                  estimate{iNode,k}.pmufseq = [xiNext(5) 1-xiNext(5)]';
                end
              otherwise
                error('unknown merge method')
            end
          end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTRAL NODE BEGIN
        if centralDecision
          if k == 1
            estimateLSS{k}.pmupseq = kron(estimate{1,k}.pmupseq,estimate{2,k}.pmupseq); %(1,1) , (1,2) , (2,1) , (2,2) % combine local predictions
            estimateLSS{k}.likelihood = kron(estimate{1,k}.likelihood,estimate{2,k}.likelihood); % combine local likelihoods
            estimateLSS{k}.pmufseq = estimateLSS{k}.pmupseq .* estimateLSS{k}.likelihood; % multiply by likelihood
            estimateLSS{k}.pmufseq = estimateLSS{k}.pmufseq/sum(estimateLSS{k}.pmufseq); % normalize
            % calculate probabilties for local nodes
            estimate{1,k}.pmufseq_central = kron(eye(2),ones(1,2))*estimateLSS{k}.pmufseq; %[1 1 0 0;0 0 1 1] * [(1,1) (1,2) [2,1) (2,2)]'
            estimate{2,k}.pmufseq_central = kron(ones(1,2),eye(2))*estimateLSS{k}.pmufseq; %[1 0 1 0;0 1 0 1] * [(1,1) (1,2) [2,1) (2,2)]'
          else
            estimateLSS{k}.pmupseq_full = LSS.P(:).* kron(estimateLSS{k-1}.pmufseq,ones(4,1)); % predict 16 probabilities
            estimateLSS{k}.likelihood_full = [kron(estimate{1,k}.likelihood(1:2),estimate{2,k}.likelihood(1:2));kron(estimate{1,k}.likelihood(1:2),estimate{2,k}.likelihood(3:4));...
              kron(estimate{1,k}.likelihood(3:4),estimate{2,k}.likelihood(1:2));kron(estimate{1,k}.likelihood(3:4),estimate{2,k}.likelihood(3:4))]; % combine local likelihoods
            estimateLSS{k}.pmufseq_full = estimateLSS{k}.pmupseq_full.*estimateLSS{k}.likelihood_full; % multiply by likelihood
            estimateLSS{k}.pmufseq_full = estimateLSS{k}.pmufseq_full/sum(estimateLSS{k}.pmufseq_full); % normalize
            estimateLSS{k}.pmufseq = kron(ones(1,4),eye(4))*estimateLSS{k}.pmufseq_full; % merging
            % calculate probabilties for local nodes
            estimate{1,k}.pmufseq_central = kron(eye(2),ones(1,2))*estimateLSS{k}.pmufseq; %[1 1 0 0;0 0 1 1] * [(1,1) (1,2) [2,1) (2,2)]'
            estimate{2,k}.pmufseq_central = kron(ones(1,2),eye(2))*estimateLSS{k}.pmufseq; %[1 0 1 0;0 1 0 1] * [(1,1) (1,2) [2,1) (2,2)]'
            if returnPmuToLocalnodes
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              % RE-MERGE LOCAL FILTERING ESTIMATES USING UPDATED PMUFSEQ
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              merge = kron(ones(1,2),eye(2));
              for iNode = 1:nNode
                if iNode == 1
                  tmp{iNode}.xiFilt_weights = kron(eye(2),kron(ones(1,2),kron(eye(2),ones(1,2))))*estimateLSS{k}.pmufseq_full; % obtain local weights of node 1
                else
                  tmp{iNode}.xiFilt_weights = kron(ones(1,2),kron(eye(2),kron(ones(1,2),eye(2))))*estimateLSS{k}.pmufseq_full; % obtain local weights of node 2
                end
                tmp{iNode}.xiMerged_weights = merge*tmp{iNode}.xiFilt_weights; % merged weights
                normalizer = 1./tmp{iNode}.xiMerged_weights;
                normalizer(isnan(normalizer)) = 1;
                tmp{iNode}.xiMerged_means = merge*(tmp{iNode}.xiFilt_means.*tmp{iNode}.xiFilt_weights); % merge means
                tmp{iNode}.xiMerged_means = tmp{iNode}.xiMerged_means.*normalizer; % normalize merged means
                tmp{iNode}.xiMerged_vars = merge*((tmp{iNode}.xiFilt_vars + (tmp{iNode}.xiFilt_means-tmp{iNode}.xiMerged_means([1 2 1 2],:)).^2).*tmp{iNode}.xiFilt_weights); % merge vars
                tmp{iNode}.xiMerged_vars = tmp{iNode}.xiMerged_vars.*normalizer; % normalize merged vars
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PASSING OUTPUT
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                estimate{iNode,k}.xfseq = tmp{iNode}.xiMerged_means';
                estimate{iNode,k}.Pxxfseq(1,1,1) = tmp{iNode}.xiMerged_vars(1,:);
                estimate{iNode,k}.Pxxfseq(1,1,2) = tmp{iNode}.xiMerged_vars(2,:);
                estimate{iNode,k}.pmufseq = tmp{iNode}.xiMerged_weights;
              end
            end
          end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTRAL NODE END
        
      case 'distributed'
        % Time update (only for k>1)
        if k == 1
          estimateLSS{1,k}.pmupseq = kron(estimate{1}.pmupseq,estimate{2}.pmupseq); %(1,1) , (1,2) , (2,1) , (2,2) % combine local predictions
        else % time update in k=2,3,... only
          switch mergeMethod
            %case 'GPB1'
              %estimate(:,k) = lssmmkfp_1D_GPB1(u(:,k-1),estimate(:,k-1),node,nodeStates);
            %case 'IMM'
              %estimate(:,k) = lssmmkfp_1D_GPB2(u(:,k-1),estimate(:,k-1),node,nodeStates);
            case 'GPB2'
              [estimate(:,k),estimateLSS(:,k)] = lssmmkfp_1D_GPB2(u(:,k-1),estimate(:,k-1),estimateLSS(:,k-1),node,LSS,nodeStates);
            otherwise
              error('unknown merge method')
          end
        end
        % Measurement update
        switch mergeMethod
          %case 'GPB1'
            %estimate(:,k) = lssmmkff_GPB1(y(:,k),estimate(:,k),node,nodeStates);
          %case 'IMM'
            %estimate(:,k) = lssmmkff_IMM(y(:,k),estimate(:,k),node,nodeStates);
          case 'GPB2'
            [estimate(:,k),estimateLSS(:,k)] = lssmmkff_GPB2(y(:,k),estimate(:,k),estimateLSS(:,k),node,LSS,nodeStates);
          otherwise
            error('unknown merge method')
        end
      otherwise
        error('unknown architecture')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute decisions at individual nodes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iNode = 1:nNode
      if centralDecision && ~returnPmuToLocalnodes
        [~,d(iNode,k)] = max(estimate{iNode,k}.pmufseq_central);
      else %decentralized or hierarchical and pmufseq computed using the informaiton from the central node
        [~,d(iNode,k)] = max(estimate{iNode,k}.pmufseq);
      end
    end
  end
  toc
  J(imc) = sum(lambda.^t.*sum(d~=mu,1));
  fprintf('GLOB AVG: %f\n ',mean(J(1:imc)))
  %fprintf('CUR: %f\n',J(imc))
  TOCtime = toc/nmc;
  %plot(x')
  %pause
end

fprintf('MC estimate of criterion value: %f\n',mean(J))
%fprintf('MC estimate of criterion value: %f\n',mean(J_DBD))

Jboot = bootstrp(100,@mean,J);
mean(Jboot)
var(Jboot)

% Save results to a text file
fid = fopen('loggedResult.txt','a+');
%fid = fopen('loggedResult_time.txt','a+');
fprintf(fid,'%s\t %s\t design:%s\t on-line:%s\t approxMethod:%s\t centralDecision:%s\t returnPmu:%s\n',datetime('now'),architecture,designMergeMethod,mergeMethod,localP,mat2str(centralDecision),mat2str(returnPmuToLocalnodes));
fprintf(fid,'MC:\t\t %u\n',nmc);
fprintf(fid,'J:\t\t %f\n',mean(J));
fprintf(fid,'mean(Jb):\t %f\n',mean(Jboot));
fprintf(fid,'var(Jb):\t %f\n',var(Jboot));
fprintf(fid,'Time:\t %f\n',TOCtime); 
fprintf(fid,'======================================================================================\n');
fclose(fid)

% Save source file and workspace to files
datestamp = datestr( now, 'yyyymmdd_HHMMSS');
sourcename = [datestamp '.m'];
% Store current source file.
s = dbstack('-completenames');
copyfile( s(1).file, sourcename );
save(datestamp)

if debug
  figure
  plot(t,[u{1,:}],t,[u{2,:}])
  grid on
  xlabel('time')
  ylabel('input')

  figure
  plot(t,[y{1,:}],t,[y{2,:}])
  grid on
  xlabel('time')
  ylabel('output')

  % figure
  % plot(t,x')
  % grid on
  % hold on
  % plot(t,xest1','.')
  % plot(t,xest2','o')
  % xlabel('time')
  % ylabel('state')

  figure
  subplot(2,1,1)
  plot(t,mu(1,:))
  hold on
  grid on
  plot(t,d(1,:),'.')

  subplot(2,1,2)
  plot(t,mu(2,:))
  hold on
  grid on
  plot(t,d(2,:),'.')

  figure
  histogram(J)
end
