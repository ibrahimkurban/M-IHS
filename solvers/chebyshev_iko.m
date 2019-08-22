function [x, maxit, xx, flopc] = chebyshev_iko(A,b,r,maxit)
%%CHEBYSHEV_IKO chebyshev implementation for LSRN algortihm refer to
% Meng, Xiangrui, Michael A. Saunders, and Michael W. Mahoney. 
% "LSRN: A parallel iterative solver for strongly over-or underdetermined systems." 
% SIAM Journal on Scientific Computing 36.2 (2014): C95-C118.
%
%   [x, maxit, xx, flopc] = chebyshev_iko(A,b,r,maxit)
%   r is the d/m
%% momentum weights
Ksup    = 1/(1-sqrt(r))^2;% + 0.28;
Kinf    = 1/(1+sqrt(r))^2;% - 0.030;

dc = (Ksup + Kinf)/2;
cc = (Ksup - Kinf)/2;

%% init
[n,d]   = size(A);
xx      = zeros(d, maxit);
r       = b;
x       = zeros(d,1);
v       = zeros(d,1);

for i=1:maxit
    %beta
    switch i
        case 0
            beta = 0;
        case 1
            beta = 0.5*(cc/dc)^2;
        otherwise
            beta = (alpha*cc*0.5)^2;
    end
    
    %alpha
    switch i
        case 0
            alpha = 1/dc;
        case 1
            alpha = dc - cc^2/(2*dc);
        otherwise
            alpha = 1/(d - alpha*c^2/4);
    end
    
    %update
    v = beta*v + A'*r;
    x = x + alpha*v;
    r = r - alpha*(A*v);
    
    %save
    xx(:,i) = x;
end

%% flop count refer to lightspeed malab packet
if(nargout > 3)
    f_iter = 20 + 4*d + 4*n*d + 2*n;
    flopc  = f_iter*[1:maxit];
end
