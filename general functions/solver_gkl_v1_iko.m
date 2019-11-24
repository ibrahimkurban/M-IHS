function [V,B] =solver_gkl_v1_iko(A, u, k)
%%SOLVER_GKL_V1 produces lower bidiagonal matrix and applies
%%reorthogonalization on the right base vector
%
%
% [V,B] =solver_gkl_v1_iko(A, u, k)
%

[n,d]   =size(A);
aa      = zeros(k,1);
bb      = zeros(k,1);
V       = zeros(d,k);
% U       = zeros(n,k+1);

%first iteraiton
u       = u/norm(u);
v       = A'*u;
aa(1)   = norm(v);
v       = v/aa(1);

% U(:,1)       = u;
% UU      = u;

u       = A*v - aa(1)*u;
% u       = u - UU*(UU'*u);
bb(1)   = norm(u);
u       = u/bb(1);


V(:,1)  = v;
VV      = v;
% U(:,2)  = u;
% UU      = U(:,1:2);
for i = 2:k
    v       = A'*u - bb(i-1)*v;
    [v,aa(i)] = reorth(VV,v,norm(v),[],1/sqrt(2),0);
    v       = v/aa(i);
    if(aa(i) == 0)
        i = i-1;
        break
    end
    %     v       = v - VV*(v'*VV)';
    %     aa(i)   = norm(v);
    %     v       = v/aa(i);
    
    u       = A*v - aa(i)*u;
    %     [u,bb(i)] = reorth(UU,u,norm(u),[],1/sqrt(2),0);
    %     u       = u/bb(i);
    %u       = u - UU*(UU'*u);
    bb(i)   = norm(u);
    u       = u/bb(i);
    %         if(bb(i) == 0)
    %         i = i-1;
    %         break
    %     end
    V(:,i)  = v;
    VV      = V(:,1:i);
    %     U(:,i+1)  = u;
    %     UU      = U(:,1:i);
end
B = (spdiags([aa, bb],[0, -1], i+1,i));

end




