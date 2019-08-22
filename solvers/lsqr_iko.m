function [x, xx] = lsqr_iko(A,b,lam,x,tol,maxit)
%%LSQR_IKO is my implementation of LSQR in which the QR step
% is used. the method uses lower bidiagonalization instead
% of upper bidiagonalization algorithm
%
% [x, xx] = lsqr_iko(A,b,lam,x,tol,maxit)
%

%% store states
if(nargout > 1)
    xx          = zeros(length(x), maxit);
end
%% bidiagonalization init.
beta1       = norm(b);
phibar      = beta1;
u           = b/phibar;
v           = A'*u;
alpha       = norm(v);
v           = v/alpha;
rhobar      = alpha;
w           = v;
lamres      = 0;
i           = 0;
%%main loop
while(i < 2 || (i < maxit && relres > tol))
    i = i+1;
    
    %% bidiagonalization
    u       = A*v - alpha*u;
    beta    = norm(u);
    u       = u/beta;
    
    
    v       = A'*u - beta*v;
    alpha   = norm(v);
    v       = v/alpha;
    
    %% regularization shift
    rhotil  = sqrt(rhobar^2 + lam);
    cstil   = rhobar/rhotil;
    sntil   = sqrt(lam)/rhotil;
    psi     = sntil*phibar;
    phibar  = cstil*phibar;
    
    %% Givens rotation for QR and 
    rho     = sqrt(rhotil^2 + beta^2);
    cs      = rhotil/rho;
    sn      = beta/rho;
    theta   = sn*alpha;
    rhobar  = cs*alpha; 
    phi     = cs*phibar;
    phibar  = -sn*phibar;
    %% solution
    x       = x + (phi/rho)*w;
    w       = v - (theta/rho)*w;
    
    %% error
    lamres  = lamres + psi^2;
    relres  = sqrt(lamres + phibar^2)/beta1; %sqrt(norm(Ax - b)^2 + lam*norm(x)^2);
    %% store
    if(nargout > 1)
        xx(:,i)     = x;
    end
end
end