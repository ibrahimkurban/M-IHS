function [x,xx,time,flopc] = lsrn_v2_lsqr_iko(A,b,lam,m,tol,x1,maxit,params)
%%LSRN_v2_LSQR_IKO imlementation of the second algorithm in paper:
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney.
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems."
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%  [x,xx,time,flopc] = lsrn_v2_lsqr_iko(A,b,lam,m,tol,x1,maxit,params)
%
% params.SA = sketch matrix otherwise it is Gaussian
% flopc is flop count
%

[n,d] = size(A);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time, f_rp] = generate_SA_lsrn(A',m, true);
else
    if(~isfield(params, 'SA'))
        [SAt, rp_time, f_rp] = generate_SA_lsrn(A',m, true);
    end
end
tic

%% SVD decomposition
if(lam == 0)
    [~, S, V] = svd(SAt,'econ');
    f_svd     = 2*m*n^2 + 11*n^3;
else
    [~, S, V] = svd([SAt;sqrt(lam)*eye(n)],'econ');
    f_svd     = 2*(m+n)*n^2 + 11*n^3;
end
%preconditioner
M   = V*(S^-1);
MA  = M'*A;
Mb  = M'*b;

%% LSQR Solver
[x, iter, xx, f_lsqr] = lsqr_iko(MA,Mb,0,tol,maxit);

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
flopc   = f_rp + f_lsqr + f_svd + iter*(n+2*n^2);


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