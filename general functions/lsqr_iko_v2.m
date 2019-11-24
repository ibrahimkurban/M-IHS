function [x, xx, tol] = lsqr_iko_v2(A,b,x,iter)
%%LSQR_IKO_v2 is my implementation of LSQR in which the QR step
% is avoided. the method uses upper bidiagonalization instead
% of lower bidiagonalization algorithm
%
% x = lsqr_iko_v2(A,b,x,iter)
%
%

%% store states
if(nargout > 1)
    xx          = zeros(length(x), iter);
end
%% bidiagonalization init.
v           = A'*b;
theta       = norm(v);
v           = v/theta;

p           = A*v;
rho         = norm(p);
p           = p/rho;

w           = v;
phi         = theta/rho;
x           = x + (phi/rho)*w;
tol         = norm(b)^2 - phi^2;
%%main loop
for i = 2:iter
    %% bidiagonalization (lanczos)
    v       = A'*p - rho*v;
    theta   = norm(v);
    v       = v/theta;
    
    w       = v - (theta/rho)*w;
    
    p       = A*v - theta*p;
    rho     = norm(p);
    p       = p/rho;
    
    %% solution
    phi     = -(phi*(theta/rho));
    x       = x + (phi/rho)*w;
    tol     = tol - phi^2;
    %% store
    if(nargout > 1)
        xx(:,i)     = x;
    end
end
end