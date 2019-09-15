function [x_ls_oracle, par_ls_oracle, err_ls_oracle, err_lam, lambdas] = LS_ridge(U, sig, V, b, x1,obj )
%%LS_METRIC calcualtes RegLS solution and Tikhonov regularization parameter
%%by using given metric.
%
%   Author  : Ibrahim Kurban Ozaslan
%   Date    : 06.12.2018
%
%   [x_ls_oracle, par_ls_oracle, err_ls_oracle] = LS_oracle(A, b, x1, x0 )
%
%   x1 : prior
%   obj : selection metric: function object whose input is x estimate
%


%% take SVD of matrix and variables
sig2    = sig.^2;

g       = U'*b - sig.*(V'*x1);
%% function to minimize
fun     = @(lambda)(obj(x1+ V*((sig./(sig2 + lambda)).*g)));

%% Parameter Grid iteration
if(nargout <= 3)
    %minimization
    options             = optimset('TolX', 1e-14, 'Display','off');
    par_ls_oracle       = fminbnd(fun, eps, 1e8, options);
else
    %% Search Grid num of points at 1. and 2. iterate
    grid1 = 100;
    grid2 = 100;
    lambdas1    = logspace(-16, log10(max(sig)), grid1);
    [par_ls_oracle, err_ls_oracle, err_lam, lambdas] = grid_search(fun, lambdas1, grid2);
end

%% oracle solution
x_ls_oracle     = x1 + V*( (sig./(sig2 + par_ls_oracle)).*g );

end