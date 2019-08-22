function [x,xx,time,flopc] = lsrn_cheby_iko(A,b,lam,m,tol,x1,maxit,params)
%%LSRN_IKO imlementation of the algorithm in paper:
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney. 
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems." 
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%  [x,xx,time,flopc] = lsrn_iko(A,b,lam,m,tol,x1,maxit,params)
%
% params.SA = sketch matrix otherwise it is Gaussian
% flopc is flop count

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time, f_rp] = generate_SA_lsrn(A,m, true);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time, f_rp] = generate_SA_lsrn(A,m, true);
    end
end
tic

%% SVD decomposition
if(lam == 0)
   [~, S, V] = svd(SA,'econ');
   f_svd     = 2*m*d^2 + 11*d^3;
else
   [~, S, V] = svd([SA;sqrt(lam)*eye(d)],'econ');
   f_svd     = 2*(m+d)*d^2 + 11*d^3;
end
%preconditioner
N   = V*(S^-1);
AN  = A*N;

%% LSQR Solver
% [x, ~, xx, f_lsqr] = lsqr_pre_iko([A; sqrt(lam)],[b; zeros(d)],R,tol,maxit);
[x, iter, xx, f_cheby] = chebyshev_iko(AN,b,d/m,maxit);
x   = N*x;
xx  = N*xx;

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
flopc   = f_rp + f_cheby + f_svd + iter*(d+2*d^2); 


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

%% flops count (refer to lightspeed matlab packet)
if(nargout>2)
    %randn 
    f1      = n*m*19;
    %multilication
    f2      = m*d*(2*n-1);
    %sum
    flopc   = f1 + f2; 
end
end