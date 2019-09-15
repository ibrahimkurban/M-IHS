function [par, dev_est] = LS_sgcv_corrected_iko(U,sig,b,n)
%%LS_SGCV calcualtes GCV for given U,sig b
%
%
%  [par, dev_est] = LS_sgcv_corrected(U,sig,b,n)
%
%   n =  size(A,1);
%
%% take SVD of matrix and variables
[m,d]   = size(U);
sig2    = sig.^2;
Ub2      = (U'*b).^2;
delta2  = norm(b)^2 - sum(Ub2);
if((delta2) < 1e-14)
    delta2 = 0;
end

%% function to minimize
beta    = @(lam)( lam./(sig2 + lam));
obj     = @(lam)(sum((beta(lam).^2).*Ub2) + delta2)...
    /( (m-d + sum(beta(lam)))*(n-d + sum(beta(lam))) );

%% minimization of curve
options                 = optimset('TolX', 1e-4, 'Display','off');
[par_log,gcv_val]       = fminbnd(@(lam)obj(10^lam), -14, log10(sig(1)), options);
par                     = 10^par_log;
dev_est                 = sqrt(gcv_val)*sqrt(( (m-d + sum(beta(par)))*(n-d + sum(beta(par))) ));

end

