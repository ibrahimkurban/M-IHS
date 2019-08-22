function [ x, err, time, iterr] = eff_acc_rand_kaczmarz_iko(A, b, x1, maxit)
%%EFF_ACC_RAND_KACZMARZ effective randomized kacmarz algorithm at paper:
% Liu, Ji, and Stephen Wright. 
% "An accelerated randomized Kaczmarz algorithm."
% Mathematics of Computation 85.297 (2016): 153-178.
%
%
%[ x, err, time, iterr] = eff_acc_rand_kaczmarz(A, b, x1, maxit)
%
sig     = svd(A,0);
lammin  = sig(end)^2;
lam     = lammin;
%% input
[n, d] = size(A);
if(nargin < 3)
    x1      = rand(d,1);
    tol     = 1e-3;
    maxit   = length(b);
elseif(nargin < 3)
    tol     = 1e-3;
    maxit   = length(b);
elseif(nargin < 4)
    maxit   = length(b);
end

%% initialize
err      = zeros(maxit,1);
time    = zeros(maxit,1);
iterr   = zeros(maxit,1);
k       = 0;
x       = x1;
y       = x;
gamma   = 0;
c1       = n^-1;
c2       = lam/n;
tic;  %!!!!! TIC
%% MAIN LOOP
while(k < 2 || k < maxit)
    k        = k + 1;
    
    %calculate constants
    c3       = gamma^2;
    gamma    = ((c1 - c2*c3) + sqrt((c1 - c2*c3)^2 + 4*c3))*0.5;
    c3       = gamma^2;
    gamma2   = ((c1 - c2*c3) + sqrt((c1 - c2*c3)^2 + 4*c3))*0.5;
    alpha    = (n - gamma2*lam)/(gamma2*(n^2 - lam));
    
    %uni sample
    if(rem(k,n)==1), order=randperm(n); end
    i        = order(rem(k-1,n)+1);
    
    %proj
    s        = (A(i,:)*y - b(i))/sum(A(i,:).^2);
    g        = s*A(i,:)';
    gt       = (1 - alpha + alpha*gamma)*g;
    yt       = (1 -n*gamma)*alpha*x + (1 - alpha + n*alpha*gamma)*y;
    x        = y - g;
    y        = yt - gt;
    
    
    
    if(rem(k, floor(maxit/100))==0)
        err(k)  = norm(A*x-b);
        time(k)  = toc;
        iterr(k) = k;
    end
end

%% SAVE
iter    = k;
err     = nonzeros(err(1:iter));
time    = nonzeros(time(1:iter));
iterr   = nonzeros(iterr);