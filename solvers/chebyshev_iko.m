function [x, maxit, xx, flopc] = chebyshev_iko(A,b,x1,r,tol,maxit)
%%CHEBYSHEV_IKO chebyshev implementation for LSRN algortihm refer to
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney. 
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems." 
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%   [x, maxit, xx, flopc] = chebyshev_iko(A,b,x1,r,maxit)
%   r is the d/m
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%
[n,d]   = size(A);
%% momentum weights
Ksup    = 1/(1-sqrt(r))/sqrt(min(n,d)/r);% + 0.28;
Kinf    = 1/(1+sqrt(r))/sqrt(min(n,d)/r);% - 0.030;

dc = (Ksup^2 + Kinf^2)/2;
cc = (Ksup^2 - Kinf^2)/2;

%% init
xx      = zeros(d, maxit);
r       = b - A*x1;
bnrm    = norm(b);
x       = x1;
v       = zeros(d,1);
i       = 0;
while(i < 2 || (norm(r)/bnrm >= tol && i < maxit))
    i       = i+1;
    %beta
    switch i-1
        case 0
            beta = 0;
        case 1
            beta = 0.5*(cc/dc)^2;
        otherwise
            beta = (alha*cc*0.5)^2;
    end
    
    %alpha
    switch i-1
        case 0
            alha = 1/dc;
        case 1
            alha = dc - cc^2/(2*dc);
        otherwise
            alha = 1/(dc - alha*cc^2/4);
    end
    
    %update
    v = beta*v + A'*r;
    x = x + alha*v;
    r = r - alha*(A*v);
    
    %save
    xx(:,i) = x;
end

%% flop count refer to lightspeed malab packet
if(nargout > 3)
    f_iter = 20 + 4*d + 4*n*d + 2*n;
    flopc  = f_iter*[1:maxit];
end

end