function [x, xx,par, B, V] = lsqr_iko_reg_test_v3(A,b,x,iter)
%%LSQR_IKO is my implementation of LSQR in which the QR step
% is used. the method uses lower bidiagonalization instead
% of upper bidiagonalization algorithm
%
% x = lsqr_iko(A,b,x,iter)
%

%% store states
aa          = zeros(iter+1,1);
bb          = zeros(iter+1,1);
par         = zeros(iter+1,1);
xx          = zeros(length(x), iter+1);
V           = zeros(length(x), iter);
%% bidiagonalization init.
[n,d]       = size(A);
beta1       = norm(b);
u           = b/beta1;
%% first iter
v = A'*u;
alpha   = norm(v);
v       = v/alpha;

u       = A*v - alpha*u;
beta    = norm(u);
u       = u/beta;

aa(1)   = alpha;
bb(1)   = beta;
V(:,1)  = v;
%%main loop
for i = 2:iter+1
    %% bidiagonalization
    if(i == 1)
        v = A'*u;
    else
        v       = A'*u - beta*v;
    end
    alpha   = norm(v);
    v       = v/alpha;
    
    u       = A*v - alpha*u;
    beta    = norm(u);
    u       = u/beta;
    
    aa(i)   = alpha;
    bb(i)   = beta;
    V(:,i)  = v;
    
    %% SVD of bidiag form
    B               = full(spdiags([aa(1:i), bb(1:i)],[0, -1], i+1,i));
    [Ub,sigb,Vb]    = svd(B);
    sigb            = diag(sigb);
    Utb             = Ub(1,:)'*beta1;
    %% GCV
    gamma   = @(lam)(lam./(sigb.^2+lam));
    delta0  = Utb(end)^2;
    fun     = @(lam)(norm(gamma(10^lam).*Utb(1:end-1))^2 + delta0)/(n-size(B,2) + sum(gamma(10^lam)))^2;
    [par_log, fval,~, out] = fminbnd(fun, -15, log10(sigb(1)), []);
    
    par(i)  = 10^par_log;
    %% sol
    x       = V(:,1:i)*(Vb*((sigb./(sigb.^2 + par(i))).*Utb(1:i)));
    
    xx(:,i) = x;
end
xx      = xx(:, 2:end);
par     = par(2:end);
end