function [x,xx,time] = lsrn_iko(A,b,lam,m,x1,tol,maxit,params)
%%LSRN_IKO imlementation of the algorithm in paper:
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney. 
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems." 
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%  [x,xx,time] = lsrn_iko(A,b,lam,m,x1,tol,maxit,params)
%
% params.SA = sketch matrix otherwise it is Gaussian
%

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time] = generate_SA_lsrn(A,m, true);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time] = generate_SA_lsrn(A,m, true);
    end
end
tic

%% SVD decomposition
if(lam == 0)
   [~, S, V] = svd(SA,'econ');
else
   [~, S, V] = svd([SA;sqrt(lam)*eye(d)],'econ');
end
%preconditioner
R = V*(S.^-1);

%% LSQR Solver
[x, ~, xx] = lsqr_pre_iko([A; sqrt(lam)],[b; zeros(d)],R,x1,tol,maxit);
time = toc+rp_time;


end






%% IF SA is not provided
function [SA, time, flopc] = generate_SA_lsrn(A,m)
%%GENERATE_SA generates gauusian random matrix
%
%   [SSA, time] = random_matrix_v2(A,m)
%
%   E[SA'SA] = I

[n,d]   = size(A);
tic;
S       = randn(m,n);
SA      = S*A;
time    = toc;

%% flops if you do not have lightmaster packet do not use this output
if(nargout>2)
    flopc = 0;
    flopc = flopc + flops_randnorm(m,n);
    flopc = flopc + flops_mul(m,n,d);
end
end