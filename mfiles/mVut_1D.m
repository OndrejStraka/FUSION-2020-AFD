function mV_sa = mVut(innerEdges,V,xi,u,transition_h,disturbance_gs_h,kappa)
% MVUT  Mean value of the Bellman function when random variable has
% Gaussian sum pdf.
%
% mV_SA = MVUT(states,inner_edges,V,s,a,transition_h,disturbance_gs_h) returns
% the mean value of the Bellman function V for the given state S and the
% given action A. The Bellman function V is given over a discrete grid
% defined by INNEREDGES. Function handles TRANSITION_H and
% DISTURBANCE_GS_H contain functions that compute state transition and
% parameters of the Gaussian sum, respectivelly.

%keyboard

if nargin < 7
    kappa = 3;
end

% The number of inner edges at individual dimensions
nE = cellfun(@length,innerEdges);

nx = 1;

[weight,mx,Px] = disturbance_gs_h(xi,u);

% weights for all y's that will be used
W = kron(weight,[kappa 0.5 0.5]/(nx + kappa));

% Find indices of nonzero weights
idx_nz_W = find(W>0);

c = sqrt(nx + kappa);
M = [0 -c c]'*sqrt(squeeze(Px)') + mx; % matrix of sigmapoints (columns are sets of sigmapoints for each mx)
y = M(:)';
V_sax = V(subidx2linidx(nE+1,aggregationidx(transition_h(xi,u,y(:,idx_nz_W)),innerEdges)));
mV_sa = W(idx_nz_W)*V_sax;
