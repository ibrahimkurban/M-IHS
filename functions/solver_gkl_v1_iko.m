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

%first iteraiton
u       = u/norm(u);
v       = A'*u;
aa(1)   = norm(v);
v       = v/aa(1);

u       = A*v - aa(1)*u;
%u       = u - UU*(UU'*u);
bb(1)   = norm(u);
u       = u/bb(1);


V(:,1)  = v;
VV      = v;
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
    %u       = u - UU*(UU'*u);
    bb(i)   = norm(u);
    u       = u/bb(i);
    
    V(:,i)  = v;
    VV      = V(:,1:i);
end
B = (spdiags([aa, bb],[0, -1], i+1,i));
% figure; subplot(1,2,1);imagesc(abs(U(:,1:i)'*U(:,1:i)), [1e-15 1e-9]); title('U^TU'); axis square;
% subplot(1,2,2); imagesc(abs(Vt(:,1:i)'*Vt(:,1:i)), [1e-15 1e-9]); title('V^TV'); axis square;
end




