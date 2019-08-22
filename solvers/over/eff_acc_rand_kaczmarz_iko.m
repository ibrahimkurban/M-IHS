function [ xx, time, iter, flopc] = eff_acc_rand_kaczmarz_iko(A, b, x1, maxit)
%%KACZMARZ standard kaczmarz algorithm to solve overdetermined not
%%necessarily sysmetric Ax = b
% Liu, Ji, and Stephen Wright. 
% "An accelerated randomized Kaczmarz algorithm." 
%  Mathematics of Computation 85.297 (2016): 153-178.
%
%[ xx, time, iter, flopc] = eff_acc_rand_kaczmarz_iko(A, b, x1, maxit)
%


%% initialize
[n,d]   = size(A);
xx      = zeros(d,100);
iter    = zeros(100,1);
flopc   = zeros(100,1);
mid_iter = ceil(maxit/100);
x       = x1;
tic;  %!!!!! TIC

%% normalization constant
if(maxit > n)
    normA = sum(A.^2,2);
end
%% FIND LAMMIN
K2  = ceil(maxit/10);
K1  = max(1, K2 - 10*n);
%stage 1
for l = 1:K2+1
    i        = randi(n,1);
    x        = x + A(i,:)'*((b(i) - A(i,:)*x)/normA(i));
    %K1 iterater
    if(l == K1+1)
        err1 = norm(A*x - b);
        flops(flops + flops_mul(n,d,1) + flops_sum(b));
    end
    
    %record
    if(rem(l, floor(maxit/100))==0)
        ii           = round(l/floor(maxit/100));
        xx(:,ii)     = x;
        iter(ii)     = l;
        flopc(ii)    = (4*d+9)*mid_iter;
    end
end
%K2 iterate
err2 = norm(A*x - b);
%min lam
lam = (n*(1 - (err2/err1)^(0.5/(K2-K1))));
%% MAIN LOOP
order   = randperm(n);
k       = K2+1;
y       = x;
c1      = n^-1;
c2      = lam/n;
c3      = 0;
while(k < maxit)
    k        = k + 1;
    
    %calculate constants
    gamma    = ((c1 - c2*c3) + sqrt((c1 - c2*c3)^2 + 4*c3))*0.5;
    c3       = gamma^2;
    gamma2   = ((c1 - c2*c3) + sqrt((c1 - c2*c3)^2 + 4*c3))*0.5;
    alpha    = (n - gamma2*lam)/(gamma2*(n^2 - lam));
    
    %uni sample
    if(rem(k,n)==1), order=randperm(n); end
    i        = order(rem(k-1,n)+1);
    
    %proj
    s        = (A(i,:)*y - b(i))/normA(i);
    g        = s*A(i,:)';
    gt       = (1 - alpha + alpha*gamma)*g;
    yt       = (1 -n*gamma)*alpha*x + (1 - alpha + n*alpha*gamma)*y;
    x        = y - g;
    y        = yt - gt;
    
    %record
    if(rem(k, floor(maxit/100))==0)
        ii          = round(k/floor(maxit/100));
        xx(:,ii)     = x;
        iter(ii)    = k;
        flopc(ii)  = (60 + 9*d)*mid_iter;
    end
end
flopc = cumsum(flopc) + 6*n*d;
time = toc;