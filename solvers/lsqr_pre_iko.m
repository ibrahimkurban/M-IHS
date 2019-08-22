function [x, i, xx, flopc] = lsqr_pre_iko(A,b,lam,R,tol,maxit)
%%LSQR_IKO is my implementation of LSQR in which the QR step
% is used. the method uses lower bidiagonalization instead
% of upper bidiagonalization algorithm
%
% [x, i, xx, flopc] = lsqr_pre_iko(A,b,lam,R,tol,maxit)
%

%% store states
[n,d]           = size(A);
if(nargout > 2)
    xx          = zeros(d, maxit);
end

%% bidiagonalization init.

beta1       = norm(b);          
phibar      = beta1;                                         
u           = b/phibar;         
v           = R'\(A'*u);        
alpha       = norm(v);          
v           = v/alpha;
rhobar      = alpha;
w           = v;
lamres      = 0;
x           = zeros(d,1);
i           = 0;
%%main loop
while(i < 2 || (i < maxit && relres > tol))
    i = i+1;
    
    %% bidiagonalization
    u       = A*(R\v) - alpha*u;
    beta    = norm(u);
    u       = u/beta;
    
    
    v       = R'\(A'*u) - beta*v;
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
    if(nargout > 2)
        xx(:,i)     = x;
    end
end
%back preconditioning
x = R\x;
if(nargout > 2)
    xx = R\xx;
end

%% flop count refer to lightspeed malab packet
if(nargout > 3)
    f_iter = 2*d^2 + 4*n*d + 5*n + 22*d + 70;
    f_init = 2*n*d + 2*d^2 + 16*d + 3*n + 7;
    flopc  = f_iter*[1:i]+f_init;
end
end