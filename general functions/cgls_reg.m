function [x,xx,time,iter] = cgls_reg(A,b,tol,maxit,x1, lam)
%CGLS Conjugate Gradient Least Squares
%   22 Jan 2013: Updated syntax and documentation.
%                Folkert Bleichrodt (CWI).
fprintf('CGLS pre runs...');
[~,n] = size(A);
x = x1;
% if(nargin <6)
%     M = eye(n);
%     pre = 0;
% end
xx      = zeros(n, maxit);
time    = zeros(1,maxit);
tic;
r = b - A*x;
% s = M'\(A'r);
s = A'*r - lam*x;


% Initialize
p      = s;
norms0 = norm(s);
gamma  = norms0^2;
normx  = norm(x);
xmax   = normx;
iter   = 0;
flag   = 0;



%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
while (iter < maxit) && (flag == 0)
    
    iter = iter+1;
%     q = A*(M\p);
    q = A*p;
    
    delta = norm(q)^2 + lam*norm(p)^2;
    if delta == 0, delta = eps; end
    alpha = gamma / delta;
    
    x     = x + alpha*p;
    r     = r - alpha*q;
    
    
%     s = M'\(A'*r);
    s = A'*r - lam*x;
    norms  = norm(s);
    gamma1 = gamma;
    gamma  = norms^2;
    beta   = gamma / gamma1;
    p      = s + beta*p;
    
    % Convergence
    normx = norm(x);
    xmax  = max(xmax, normx);
    flag  = (norms <= norms0 * tol) || (normx * tol >= 1);
    
    xx(:,iter) = x;
    time(iter) = toc;
end % while
% x           = M\x;
time(iter)   = toc;
xx          = xx(:,1:iter);
time        = time(1:iter);
% xx          = M\xx;
fprintf('%2.2f sec elapsed\n',toc);
end