function [gamma,Q,V,iterationConvergence,iterationTimes] = valueiteration_discraggr2(states,innerEdges,inputs,mL_h,mV_h,isStateAdmissible_h,isInputAdmissible_h,lambda,varargin)
%function [gamma,Q,V,max_dV,iterationTimes] = valueiteration_discraggr2(states,inner_edges,inputs,mL_h,mV_h,isstateadmissible_h,lambda,epsVVI,maxiterVVI,V)



% Check the optional input arguments and set defaults
numvarargsmax = 5;
numvarargs = length(varargin);
if numvarargs > numvarargsmax
    exception = MException('MyToolbox:xxx:tooManyOptInArgs',...
        'The function ''%s'' accepts at most %i optional input arguments.',mfilename,numvarargsmax);
    throw(exception);
end

% Skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);

% Set defaults for optional input arguments
optargs = {0.01, 50, [], 'state-vectorized', false};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);

% Place optional input arguments in memorable variable names
[epsVVI,maxiterVVI,V,stateInputVectorization,verbose] = optargs{:};


%  Preallocate arrays
min_dV = nan(1,maxiterVVI);
max_dV = nan(1,maxiterVVI);
max_abs_dV = nan(1,maxiterVVI);
min_Vopt_V = nan(1,maxiterVVI);
max_Vopt_V = nan(1,maxiterVVI);
dV_symm = nan(1,maxiterVVI);
iterationTimes = nan(1,maxiterVVI);


% The number of input quantization levels at individual dimensions
nU = cellfun(@length,inputs);
nu = length(inputs);

% The total number of all discrete inputs
nInputs = prod(nU);

% Generate all inputs from quantization levels given by inputs
u = deaggregationidx(inputs,linidx2subidx(nU,1:nInputs));

% Create logical index of admissinble inputs
flagAdmissibleInputs = isInputAdmissible_h(u);

% Pick only admissible inputs
admissibleInputs = u(:,flagAdmissibleInputs);

% The number of admissible inputs
nAdmissibleInputs = size(admissibleInputs,2);

% The number of state quantization levels at individual dimensions
nS = cellfun(@length,states);

% The total number of discrete states
nStates = prod(nS);

% Generate all states from quantization levels given in states
s = deaggregationidx(states,linidx2subidx(nS,1:nStates));

% Create logical index of admissible states
flagAdmissibleStates = isStateAdmissible_h(s);

% Pick only admissible states (It is mainly because parfor loop of parallel computing toolbox  permits only indexing without gaps)
admissibleStates = s(:,flagAdmissibleStates);

% The number of admissible states
nAdmissibleStates = size(admissibleStates,2);

% Check that the lenght of array V function is nStates - this check makes sense
% only when initial values of V function are supplied as an input variable


if isempty(V)
    % If the dimension of the supplied array of initial values is
    % inconsistent, it is set to zeros for admissible states
    V = inf(nStates,1);
    V(flagAdmissibleStates) = 0;
end

sizV = size(V);
if numel(sizV)>2 || sizV(1) ~= nStates || sizV(2) ~= 1
    % If the dimension of the supplied array of initial values is
    % inconsistent, it is set to zeros for admissible states
    V = inf(nStates,1);
    V(flagAdmissibleStates) = 0;
    
    warning('MyToolbox:valueiteration:override',...
        'The provided value function has inconsistent dimensions. Reseting to zero')
end


% Initialize the value iteration algorithm
nextIter = true;
iter = 0;


while nextIter
    tic
    iter = iter + 1;
    if verbose
        fprintf('Iteration no. %i',iter)
    end
    
    % Initialize auxiliary value function Q
    Q_iter = zeros(nAdmissibleStates,nAdmissibleInputs);
    
    % Compute auxiliary value function V and \gamma only for admissible
    % states
    switch stateInputVectorization
        case 'non-vectorized'
            % State nonvectorized version
            for ui = 1:nAdmissibleInputs
                for si = 1:nAdmissibleStates
                    Q_iter(si,ui) = mL_h(admissibleStates(:,si),admissibleInputs(:,ui)) + lambda*mV_h(innerEdges,V,admissibleStates(:,si),admissibleInputs(:,ui));
                end
            end
        case 'state-vectorized'
            % State vectorized version - intended for a deterministic model
            for ui = 1:nAdmissibleInputs
                Q_iter(:,ui) = mL_h(admissibleStates,admissibleInputs(:,ui)) + lambda*mV_h(innerEdges,V,admissibleStates,admissibleInputs(:,ui));
            end
        case 'input-vectorized'
            % Input vectorized version - intended for a deterministic model
            for si = 1:nAdmissibleStates
                Q_iter(si,:) = mL_h(admissibleStates(:,si),admissibleInputs) + lambda*mV_h(innerEdges,V,admissibleStates(:,si),admissibleInputs);
            end
        otherwise
            exception = MException('MyToolbox:xxx:unknownOption',...
                'The last optional input of function ''%s'' must be one of the following: ''non-vectorized'', ''state-vectorized'', ''input-vectorized''',mfilename);
            throw(exception);
            
    end
    
    [V_iter,u_idx] = min(Q_iter,[],2);
    
    
    % Compute bounds on the value function
    min_dV(iter) = min(V_iter - V(flagAdmissibleStates));
    max_dV(iter) = max(V_iter - V(flagAdmissibleStates));
    max_abs_dV(iter) = max(abs([min_dV(iter) max_dV(iter)]));
    min_Vopt_V(iter) = min_dV(iter)/(1-lambda);
    max_Vopt_V(iter) = max_dV(iter)/(1-lambda);
    dV_symm(iter) = 2*max_abs_dV(iter)*lambda/(1-lambda);
    
    % Check the stopping criterion over admissible states
    if max_abs_dV(iter) < epsVVI || iter >= maxiterVVI
        nextIter = false;
        Q = nan(nStates,nInputs);
        Q(flagAdmissibleStates,flagAdmissibleInputs) = Q_iter;
        gamma = nan(nu,nStates);
        gamma(:,flagAdmissibleStates) = admissibleInputs(:,u_idx);
    end
    V(flagAdmissibleStates) = V_iter;
    iterationTimes(iter) = toc;
    if verbose
        fprintf(' took %f [s], improvement in Value function %f\n',iterationTimes(iter),max_abs_dV(iter));
    end
    
    
    if verbose && ((mod(iter,10) == 0) || iter == 1)
%         V_r = reshape(V,nS);
%         figure
%         stepsurf(states([1,2]),innerEdges([1,2]),V_r(:,:,1))
%         view([0,90])
%         title(sprintf('Posture: %f,Iteration %d',states{3}(1)/pi*180,iter))
%         drawnow
%         
%         figure
%         stepsurf(states([1,2]),innerEdges([1,2]),V_r(:,:,47))
%         view([0,90])
%         title(sprintf('Posture: %f,Iteration %d',states{3}(47)/pi*180,iter))
%         drawnow
        
        %pause
        %u_tmp = admissibleInputs(:,u_idx);
        %        snew = transition_h(admissibleStates,u_tmp);
        %       ds = snew-admissibleStates;
        
        %              for i = 1:1:1%length(states{3})
        %                  figure
        %                  stepsurf(states([1,2]),innerEdges([1,2]),V_r)
        %                  title(sprintf('angle %f',states{3}(i)/pi*180))
        %                  view([0,90])
        %                  drawnow
        %                  %hold on
        %                  %quiver3(admissibleStates(1,:),admissibleStates(2,:),1*ones(1,nAdmissibleStates),ds(1,:),ds(2,:),0*ones(1,nAdmissibleStates),0)
        %              end
        %subidx2linidx(nS,aggregationidx([-0.15,-2.5,320/180*pi]',innerEdges))
    end
end

% Delete all nan values from the stopping criterion
iterationConvergence = [...
    min_dV(1:iter)
    max_dV(1:iter)
    max_abs_dV(1:iter)
    min_Vopt_V(1:iter)
    max_Vopt_V(1:iter)
    dV_symm(1:iter)];
iterationTimes = iterationTimes(1:iter);
