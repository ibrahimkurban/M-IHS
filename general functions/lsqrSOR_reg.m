function [ x, xx, time, iter] = lsqr_reg(A, b, tol, maxit,x1, lam )
%-----------------------------------------------------------------------
% LSQR uses an iterative (conjugate-gradient-like) method.
% fprintf('LSQR pre runs...')
% Initialize.
atol    = tol;
btol    = tol;
[m,n]   = size(A);
iter    = 0;
istop   = 0;
xxnorm  = 0;
z       = 0;
cs2     = -1;
sn2     = 0;
Anorm   = 0;
exit    = false;
% if(nargin < 6)
%     M   = eye(n);
%     pre = 0;
% end
xx      = zeros(n, maxit);
time    = zeros(1,maxit);
tic;

% Set up the first vectors u and v for the bidiagonalization.
% These satisfy  beta*u = b,  alfa*v = A'u.

u       = b(1:m);
x       = x1;
alfa    = 0;
beta    = norm(u);
if beta > 0
    u = (1/beta)*u;
%     v = M'\(A'*u);
    v = A'*u;
    alfa = norm(v);
end

if alfa > 0
    v = (1/alfa)*v;
    w = v;
end
rhobar = alfa;
phibar = beta;
bnorm  = beta;

%------------------------------------------------------------------
%     Main iteration loop.
%------------------------------------------------------------------
while ~exit
    iter = iter + 1;
    % Perform the next step of the bidiagonalization to obtain the
    % next beta, u, alfa, v.  These satisfy the relations
    %      beta*u  =  A*v  - alfa*u,
    %      alfa*v  =  A'*u - beta*v.
%     u = A*(M\v) - alfa*u;
    u = A*v - alfa*u;
    
    beta = norm(u);
    if beta > 0
        u     = (1/beta)*u;
        Anorm = norm([Anorm alfa beta lam]);
%         v = M'\(A'*u)   - beta*v;
        v = A'*u - beta*v;
        alfa  = norm(v);
        if alfa > 0,  v = (1/alfa)*v; end
    end
    
    % Use a plane rotation to eliminate the damping parameter.
    % This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
    
    rhobar1 = norm([rhobar lam]);
    cs1     = rhobar/rhobar1;
    phibar  = cs1*phibar;
    
    % Use a plane rotation to eliminate the subdiagonal element (beta)
    % of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
    
    rho     =   norm([rhobar1 beta]);
    cs      =   rhobar1/rho;
    sn      =   beta   /rho;
    theta   =   sn*alfa;
    rhobar  = - cs*alfa;
    phi     =   cs*phibar;
    phibar  =   sn*phibar;
    tau     =   sn*phi;
    
    % Update x and w.
    
    t1      =   phi  /rho;
    t2      = - theta/rho;
    x       = x      + t1*w;
    w       = v      + t2*w;
    
    % Use a plane rotation on the right to eliminate the
    % super-diagonal element (theta) of the upper-bidiagonal matrix.
    % Then use the result to estimate  norm(x).
    
    delta   =   sn2*rho;
    gambar  = - cs2*rho;
    rhs     =   phi - delta*z;
    zbar    =   rhs/gambar;
    xnorm   =   sqrt(xxnorm + zbar^2);
    gamma   =   norm([gambar theta]);
    cs2     =   gambar/gamma;
    sn2     =   theta /gamma;
    z       =   rhs   /gamma;
    xxnorm  =   xxnorm + z^2;
    
    % Test for convergence.
    % First, estimate the condition of the matrix  Abar,
    % and the norms of  rbar  and  Abar'rbar.
    Arnorm  =   alfa*abs(tau);
    rnorm   =   abs(phibar);
    
    % Now use these norms to estimate certain other quantities,
    % some of which will be small near a solution.
    
    test1   =   rnorm /bnorm;
    test2   =   Arnorm/(Anorm*rnorm);
    t1      =   test1/(1 + Anorm*xnorm/bnorm);
    rtol    =   btol + atol*Anorm*xnorm/bnorm;
    
    % The following tests guard against extremely small values of
    % atol, btol  or  ctol.  (The user may have set any or all of
    % the parameters  atol, btol, conlim  to 0.)
    % The effect is equivalent to the normal tests using
    % atol = eps,  btol = eps,  conlim = 1/eps.
    
    if iter >= maxit,   istop = 7; end
    if 1 + test2  <= 1, istop = 5; end
    if 1 + t1     <= 1, istop = 4; end
    if  test2 <= atol,  istop = 2; end
    if  test1 <= rtol,  istop = 1; end
    
    if istop > 0; exit = true; end
    xx(:,iter) = x;
    time(iter) = toc;
end
% x           = M\x;
time(iter)  = toc;
xx          = xx(:,1:iter);
time        = time(1:iter);
% xx          = M\xx;
% fprintf('%2.2f sec elapsed\n', toc);
end
