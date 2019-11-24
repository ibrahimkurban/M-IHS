function [x_ls_gcv, par, dev_est,time, obj, lambdas] = LS_gcv_iko(U, sig,V,b, x1)
%%LS_GCV calcualtes RegLS solutionwith paramter found via GCV of G.Golub.
%
%   Author  : Ibrahim Kurban Ozaslan
%   Date    : 14.07.2018
%   v2      : finds lambda by fminbnd
%
%   [x_ls_gcv, par_ls_gcv, obj_ls_gcv_par, lambdas_ls_gcv] = LS_gcv(A,b, x1)
%
%   x1 : prior
%

tic;
%% take SVD of matrix and variables
sig2    = sig.^2;

Ub      = U'*b;
Vx1     = V'*x1;
g       = Ub - sig.*Vx1;
g2      = g.^2;
delta2  = norm(b)^2 - norm(Ub)^2;
if((delta2) < 1e-14)
    delta2 = 0;
end
%% function to minimize
d       = length(x1);
n       = length(b);

MN      = (n-d).*(n>d);
beta    = @(lambda)( lambda./(sig2 + lambda));
obj     = @(lambda)((sum((beta(lambda).^2).*g2) + delta2)...
                                            /(MN + sum(beta(lambda)))^2);

%% Parameter Grid iteration
if(nargout <= 3)
    %minimization
    options             = optimset('TolX', 1e-4, 'Display','off');
    [par_log,gcv_val]          = fminbnd(@(lam)obj(10^lam), -14, log10(sig(1)), options);
    par = 10^par_log;
else
    %% Search Grid num of points at 1. and 2. iterate
    grid1 = 200;
    grid2 = 100;
    lambdas1    = logspace(-16, log10(max(sig)), grid1);
    [par, gcv_val, obj, lambdas] = grid_search(obj, lambdas1, grid2);
end
dev_est = sqrt(gcv_val)*sqrt(MN + sum(beta(par)));
%% solution
x_ls_gcv     = x1 + V*( (sig./(sig2 + par)).*g );
time = toc;
end

