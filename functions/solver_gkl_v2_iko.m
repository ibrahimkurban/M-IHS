function [R,V,D] = solver_gkl_v2_iko(A, v, k)
%%SOLVER_GKL_V2
%
%   [R,V,D] = solver_gkl_v2(A, v, k)
%
[n,d]   =size(A);
rr      = zeros(k,1);
tt      = zeros(k,1);
V       = zeros(d,k);

v       = v/norm(v);
p       = A*v;
rr(1)   = norm(p);
p       = p/rr(1);



D       = zeros(d, k);
d       = v/rr(1);
D(:,1)  = d;


V(:,1)  = v;
VV      = v;
ii      = k;
for i = 2:k
    v       = A'*p - rr(i-1)*v;
    v       = v -VV*(VV.'*v);
    tt(i)   = norm(v);
    v       = v/tt(i);
    
    p       = A*v - tt(i)*p;
    %             p  = p - PP*(PP'*p);
    rr(i)   = norm(p);
    p       = p/rr(i);
    
    d = (v - tt(i)*d)/rr(i);
    
    V(:,i) = v;
    D(:,i) = d;
    VV     = V(:,1:i);
    
    %check accuracy
    % this is due to unstable implemantation of GKL procedure
    % possible reason is limited double precision
     % this is a temporary solution, may possible be wrong
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




