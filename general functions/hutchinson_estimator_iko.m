function [k0, in_iter,flopc, time] = hutchinson_estimator_iko(SA, par, N, tol, maxit)
%EFF_RANK_SOLVER estimates trace of A(ASSA+parI)^-1A' by using hutchinson
% estimator and AA_b_solver this value corresponds to effective rank with
% respect to given parameter par
%
% [k0, in_iter,flopc] = hutchinson_estimator_iko(SA, A, par, N, tol, maxit)
%
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
    tol     = 1e-2;
    maxit   = 100;
elseif(nargin < 4)
    tol     = 1e-2;
    maxit   = 100;
elseif(nargin < 5)
    maxit   = 100;
end
[~,d]           = size(SA); 
flopc       = 0;
in_iter     = zeros(N,1);
tr_est      = 0;
tic;
%% rademacher
v                   = randi(2, d, N);
v(v==2)             = -1;
%% solver
for i=1:N
    [AA_Iiv, in_iter(i),fc]     = AA_b_solver_iko(SA, v(:,i), par, tol, maxit);
    tr_est                      = tr_est + par*sum(v(:,i).*AA_Iiv);
    flopc                       = flopc + fc(end) + 2*d;
end
time = toc;
%% expectation (average)
tr_est  = d -tr_est/N;

%% result
k0              = ceil(tr_est/10)*10;
end