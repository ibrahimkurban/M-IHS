function [x,xx,time,flopc] = lsrn_v2_lsqr_iko(A,b,lam,m,x1,tol,maxit,params)
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
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
% for under determined problems
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%
time    = zeros(maxit+3,1);
%% Check Tikhonov
[n,d] = size(A);
if(lam ~= 0 )
    A = [A/sqrt(lam), eye(n)];
    d = n+d;
end
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time, f_rp] = generate_SA_lsrn(A',m);
else
    if(~isfield(params, 'SA'))
        [SAt, rp_time, f_rp] = generate_SA_lsrn(A',m);
    end
end
time(1)  =rp_time;

%% SVD decomposition
tic;
[~, S, V] = svd(SAt,'econ');
f_svd     = 2*m*n^2 + 11*n^3;

%preconditioner
M   = V*(S^-1);
MA  = M'*A;
Mb  = M'*b;
time(2) = toc;
%% LSQR Solver
[x, iter, xx, f_lsqr, time(4:end)] = lsqr_iko(MA,Mb,0,tol,maxit);

x   = x(1:d-n)/sqrt(lam);
xx  = xx(1:d-n,:)/sqrt(lam);

%% complexity (flop count refer to lightspeed malab packet)

% flop count
flopc   = f_rp + f_lsqr + f_svd + [1:iter]*(n+2*n^2);


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