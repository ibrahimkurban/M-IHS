function [x0e, k_star, obj_par] = LS_truncated(U,sig,V, b, x1, x0)
%%LS_EFFECTIVE finds effective oracle information by using oracle
%%information. Simply, it extracts the oracle part that oracle energy is
%%higher than noise energy. While doing that TSVD anaylsis is used
%
%   Author  : Ibrahim Kurban Ozaslan
%   Date    : 14.07.2018
%   v1      : finds k* by TSVD
%
%   [x0e, par_0e, obj_par] = TSVD_effective(A, b, x1, x0 )
%
%   x1 : prior
%   x0: oracle info
%   x0e: effective oracle info

%% take SVD of matrix and variables
sigi    = sig.^-1;

Ub      = U'*b;
Vx1     = V'*x1;
g       = Ub - sig.*Vx1;

%% function to minimize
x0_til  = x0 - x1;
Vsg     = cumsum(V.*reshape(sigi.*g, 1, []), 2);
obj     = @(k)(norm(x0_til - Vsg(:,k)));

%% minimize
K       = min(length(x1), length(b));
obj1    = obj(1);
obj_par = obj1;
for k = 1:K
    objj = obj(k);
    if(objj <= obj_par)
        obj_par = objj;
        k_star = k;
    elseif(objj > obj1)
        break
    end
end

% obj_par = zeros(K,1);
% for k = 1:K
%     obj_par(k) = obj(k);
% end
% [~, k_star] = min(obj_par);


%% solution
x0e = V(:,1:k_star)*(V(:,1:k_star)'*x0);

end