function [k0, flopc,time] = hutchinson_tri_iko(R,par,N)
%EFF_RANK_SOLVER estimates trace of SA(ASSA+parI)^-1SA' by using hutchinson
% estimator and AA_b_solver this value corresponds to effective rank with
% respect to given parameter par
%
% [k0, in_iter,flopc] = hutchinson_estimator_iko(R, par, N, tol, maxit)
%
% R is R-factor of SA matrix
% N, tol, maxit are optional
%   default values are:
%       N       : 2
%       tol     : 1e-2;
%       maxit   : 100;
%
%
%% initial estimate
if(nargin < 3)
    N       = 2;
end
d           = size(R,1);  
flopc       = 0;
tr_est      = 0;
tic;
%% rademacher
v                   = randi(2, d, N);
v(v==2)             = -1;
%% solver
for i=1:N
    Riv                         = R'\v(:,i);
    tr_est                      = tr_est + par*norm(Riv)^2;
    flopc                       = flopc + 2*d + d^2;
end
time = toc;
%% expectation (average)
tr_est  = d -tr_est/N;

%% roundend above result
k0              = ceil(tr_est/10)*10;
end
