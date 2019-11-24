function [k0, in_iter,flopc] = hutchinson_estimator_iko(SA, A, par, N, tol, maxit)
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
if(nargin < 4)
    N       = 2;
    tol     = 1e-2;
    maxit   = 100;
elseif(nargin < 5)
    tol     = 1e-2;
    maxit   = 100;
elseif(nargin < 6)
    maxit   = 100;
end

%% rademacher
[n,d]               = size(A);  
v                   = randi(2, size(A,1), N);
v(v==2)             = -1;
Av                  = A'*v;
flopc               = 2*n*d*N;
%% solver
in_iter     = zeros(N,1);
tr_est      = 0;
for i=1:N
    [AA_Iiv, in_iter(i),fc]     = AA_b_solver_iko(SA, Av(:,i), par, tol, maxit);
    tr_est                      = tr_est + sum(Av(:,i).*AA_Iiv);
    flopc                       = flopc + fc + 2*d;
end

%% expectation (average)
tr_est  = tr_est/N;


%% tolerance warning
if(sum(in_iter>=maxit) > 0)
    fprintf('\n\n!!!!!!!!!! ITERATION NUMBER IS NOT SUFFICENT FOR DESIRED TOLERANCE !!!!!!!!!!\n\n')
end

%% result
k0              = ceil(tr_est/10)*10;
end