function [R,V,D] = solver_gkl_v2_iko_c(A, v, k, cgs)
%%SOLVER_GKL_V2_IKO_C produces upper bidiagonal matrix by using mex files
%%established by Larsen. For
%
%   [R,V,D] = solver_gkl_v2_iko_c(A, v, k)
%
%
% References:
%  Aake Bjorck, "Numerical Methods for Least Squares Problems",
%  SIAM, Philadelphia, 1996, pp. 68-69.
%
%  J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart,
%  ``Reorthogonalization and Stable Algorithms Updating the
%  Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
%  136, pp. 772-795.
%
% [R,V,D] = solver_gkl_v2_iko_c(A, v, k, cgs)
%
% cgs 0: modified gram schmidt
% cgs 1: classic gram schmidth default


[n,d]   =size(A);
rr      = zeros(k,1);
tt      = zeros(k,1);
V       = zeros(d,k);

%% firs iteration
v       = v/norm(v);
p       = A*v;
rr(1)   = norm(p);
p       = p/rr(1);

D       = zeros(d, k);
d       = v/rr(1);
D(:,1)  = d;
V(:,1)  = v;

%% Check input arguments.
ind     = 1:k;
alpha   = 1/sqrt(2);
ii      = k;
if(nargin <= 2)
    k = d;
    cgs = true;
elseif(nargin <= 3)
    cgs= true;
end

%% main iteration
for i = 2:k
    v           = A'*p - rr(i-1)*v;
    [v,tt(i)]   = reorth(V,v,norm(v),ind(1:i),alpha,cgs);
    v           = v/tt(i);
    
    p           = A*v - tt(i)*p;
    rr(i)       = norm(p);
    p           = p/rr(i);
    
    %forward substituion
    d = (v - tt(i)*d)/rr(i);
    
    V(:,i) = v;
    D(:,i) = d;
    
    %check accuracy
    % this is due to unstable implemantation of GKL procedure
    % possible reason is limited double precision
    % this is a temporary solution
    % need some work here
    % there exist some packets that solves this issue in C
    if(i > 2 && rr(i)*tt(i) > 10*rr(i-1)*tt(i-1))
        ii = i-2;
        break;
    end
end
R = (spdiags([rr, tt],[0, 1], ii,ii));
if(ii < k)
    V = V(:,1:ii);
    D = D(:,1:ii);
end
end




