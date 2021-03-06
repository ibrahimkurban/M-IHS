function [x,xx,time,flopc] = lsrn_cheby_iko(A,b,lam,m,x1,tol,maxit,params)
%%LSRN_IKO imlementation of the algorithm in paper:
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney.
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems."
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%  [x,xx,time,flopc] = lsrn_iko(A,b,lam,m,tol,x1,maxit,params)
%
% params.SA = sketch matrix otherwise it is Gaussian
% flopc is flop count
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%
%% Check Tikhonov
[n,d] = size(A);
time  = zeros(maxit+3,1);
if(lam ~= 0 )
    A = [A;sqrt(lam)*eye(d)];
    b = [b; zeros(d,1)];
    n = n+d;
end
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time, f_rp] = generate_SA_lsrn(A,m);
    params.k0 = d;
else
    if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_lsrn(A,m);
    else
        SA      = params.SA;
        rp_time = 0;
        f_rp    = 0;
    end
    if(~isfield(params, 'k0'))
        params.k0 = d;
    end
end
time(1) = rp_time;

%% SVD decomposition
tic;
[~, S, V] = svd(SA,'econ');
f_svd     = 2*m*d^2 + 11*d^3;

%preconditioner
N   = V*(S^-1);
AN  = A*N;
time(2) = toc;
time(3) = 0;
%% Cheybyshev Solver
[x, iter, xx, f_cheby, time(4:end)] = chebyshev_iko(AN,b,x1,d/m,tol,maxit);
tic;
x   = N*x;
xx  = N*xx;

%% complexity (flop count refer to lightspeed malab packet)
%timing
time(end) = time(end) + toc;

% flop count
flopc   = f_rp + f_cheby + f_svd + [1:iter]*(d+2*d^2);


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