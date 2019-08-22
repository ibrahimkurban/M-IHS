function [x, k] = AA_b_solver_iko(A, b, par, tol, maxit, x1)
%%AA_b solver solves the linear system (A'A + par*eye)x = b by using GKL-2
%bidiagonalization procedure which generates upper bidiagonal matrix 
%without reorthogonalization and also without any additional cost, it finds
%    relative error : |AA_Ix - b|/|b|
%
%   [x, k] = AA_b_solver(A, b, par, tol, maxit)
%
%   tol, maxit, x1 is optional
%   default values are :
%   tol         : 1e-2
%   maxit       : 100
%   x1          : zeros
%

%% initial estimate
if(nargin < 4)
    tol     = 1e-2;
    maxit   = 100;
    x1      = zeros(size(A,2),1);
elseif(nargin < 5)
    maxit   = 100;
    x1      = zeros(size(A,2),1);
elseif(nargin < 6)
    x1      = zeros(size(A,2),1);
end
%% initial step
theta1      = norm(b);
v           = b/theta1;
p           = A*v;
rho         = norm(p);
p           = p/rho;

%% shift regularization
rhobar      = sqrt(rho^2 + par);
c           = rho/rhobar;
s           = sqrt(par)/rhobar;
w           = v/rhobar;

%% inverse bidiag
phi         = theta1/rhobar;

%% sol
x           = phi*w;% + x1;

k           = 1;
while(k == 1 || (k < maxit && tol < tol_prev_it))
    k       = k+1;
    
    v       = A'*p - rho*v;
    theta   = norm(v);
    v       = v/theta;
    
    p       = A*v - theta*p;
    rho     = norm(p);
    p       = p/rho;
    
    
    %regularization shift
    parbar2     = (par + (s*theta)^2);
    thetabar    = c*theta;
    
    
    rhobar      = sqrt(rho^2 + parbar2);
    c           = rho/rhobar;
    s           = sqrt(parbar2)/rhobar;
    
    %d
    w           = (v - thetabar*w)/rhobar;
    
    %bidiag inversion
    phi         = phi*(-thetabar/rhobar);
    
    %sol
    x           = x + phi*w;
    
    %tolerance
    tol_prev_it = abs(phi*rhobar)/theta1;
end

end