function [par, dev_est] = LS_sgcv_lower_iko(U, sig,b)
%%LS_SGCV calcualtes GCV for given U,sig b
%
%
%  [par, dev_est] = LS_gcv(U, sig,b)
%
%   x1 : prior
%

%% take SVD of matrix and variables
sig2    = sig.^2;
Ub2     = (U'*b).^2;
delta2  = norm(b)^2 - sum(Ub2);
if((delta2) < 1e-14)
    delta2 = 0;
end
%% function to minimize
[m,d]   = size(U);

md      = (m-d);
beta    = @(lam)( lam./(sig2 + lam));
obj     = @(lam)((sum((beta(lam).^2).*Ub2) + delta2)...
    /(md + sum(beta(lam)))^2);

%% minimization of curve
options                 = optimset('TolX', 1e-4, 'Display','off');
[par_log,gcv_val]       = fminbnd(@(lam)obj(10^lam), -14, log10(sig(1)), options);
par                     = 10^par_log;
dev_est = sqrt(gcv_val)*sqrt(md + sum(beta(par)));

end

