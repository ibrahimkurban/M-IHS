function [ xx, time, flipflop, iter] = eff_acc_rand_kaczmarz_flop(A, b, x1, maxit)
%%KACZMARZ standard kaczmarz algorithm to solve overdetermined not
%%necessarily sysmetric Ax = b
%
%
%[ x, xx, time, iter] = kaczmarz(A, b, x1, tol, maxit)
%
fprintf('Accelerated Radomized Kaczmarz runs...')
flops(0);

%% initialize
[n,d]   = size(A);
flipflop= zeros(100,1);
xx      = zeros(d,100);
time    = zeros(100,1);
iter    = zeros(100,1);
x       = x1;
tic;  %!!!!! TIC

%% normalization constant
if(maxit > n)
    normA = sum(A.^2,2);
    flops(flops + n*flops_mul(1,d,1));
end
%% FIND LAMMIN
K2  = ceil(maxit/10);
K1  = max(1, K2 - 10*n);
%stage 1
for l = 1:K2+1
%     if(rem(l,n)==1), order=randperm(n); end
%     i        = order(rem(l-1,n)+1);
    i        = randi(n,1);
%    i        = (rem(l-1,n)+1);
    x        = x + A(i,:).'*((b(i) - A(i,:)*x)/normA(i));
    flops(flops + flops_mul(1,d,1) + 2 + flops_sum(x));
    %K1 iterater
    if(l == K1+1)
        err1 = norm(A*x - b);
        flops(flops + flops_mul(n,d,1) + flops_sum(b));
    end
    
    %record
    if(rem(l, floor(maxit/100))==0)
        ii          = round(l/floor(maxit/100));
        xx(:,ii)     = x;
        time(ii)     = toc;
        iter(ii)    = l;
        flipflop(ii) = flops;
    end
    disp(err_x(x))
end
%K2 iterate
err2 = norm(A*x - b);
flops(flops + flops_mul(n,d,1) + flops_sum(b));
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
    flops(flops + 40);
    
    %uni sample
    if(rem(k,n)==1), order=randperm(n); end
    i        = order(rem(k-1,n)+1);
    
    %proj
    s        = (A(i,:)*y - b(i))/normA(i);
    flops(flops + flops_mul(1,d,1) + 2);
    g        = s*A(i,:)';
    flops(flops + flops_sum(x));
    gt       = (1 - alpha + alpha*gamma)*g;
    flops(flops + flops_sum(g) + 3);
    yt       = (1 -n*gamma)*alpha*x + (1 - alpha + n*alpha*gamma)*y;
    flops(flops + 3*flops_sum(g) + 10);
    x        = y - g;
    y        = yt - gt;
    flops(flops + 2*flops_sum(x));
    
    %record
    if(rem(k, floor(maxit/100))==0)
        ii          = round(k/floor(maxit/100));
        xx(:,ii)     = x;
        time(ii)     = toc;
        iter(ii)    = k;
        flipflop(ii) = flops;
    end
end
fprintf('%2.2f sec elapsed\n', time(end))