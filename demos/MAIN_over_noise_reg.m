clear all
close all
clc

%%
n       = 2^16;
d       = 1000;
m       = 2*d;


TOL     = 0;
MAXIT   = 20;
% TOL2    = 1e-10;
% MAXIT2  = 35;


PROFILE = 5;%5,3
MAX     = 1;
KAPPA   = 8;%7,2
LEVEL   = 0.01;
MC      = 1;

%% GENERATE DATA
[A, b, x0, x1, dev, U, sigA, V] = generate_data(n,d,...
    'correlated',               true,...
    'matrix correlation',       0.9,...
    'matrix deviation',         5,...
    'hansen singular values',   PROFILE,...
    'hansen max',               MAX,...
    'kappa',                    KAPPA,...
    'noise level',              LEVEL,...
    'hansen input',             true,...
    'input dev',                1,...
    'prior dev',                0,...
    'prior mean',               0,...
    'structure',                false,...
    'plot',                     false);

%% ORACLE Regularized SOL
[xxOR, parOR, errOR, err_rid] = oracle_methods(U, sigA, V,b,x0,x1);
err_x0  = @(xx)(errOR{2}(xx,1));
k0      = parOR(1);
par0    = parOR(2);
fprintf('effective rank : %d\n', k0);
fprintf('Oracle reg. par: %1.2e\n', par0);

lam = par0;
%%
N_met   = 7;
err_mc  = zeros(MAXIT,N_met,MC);
flop_mc = zeros(MAXIT,N_met,MC);
names   = cell(N_met,1);
for mc = 1:MC
    fprintf('mc: %d\n', mc);
    err   = zeros(MAXIT,N_met);
    flop  = zeros(MAXIT,N_met);
    k     = 1;
    %% blendenpik
    [~,xx,~, flopc] = blendenpik_iko(A,b,lam,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'Blendenpik';
    k            = k+1;
     
    %% lsrn + lsqr
    [~,xx,~, flopc] = lsrn_lsqr_iko(A,b,lam,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'LSRN-lsqr';
    k            = k+1;
    
    %% lsrn + chebyshev
    [~,xx,~, flopc] = lsrn_cheby_iko(A,b,lam,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'LSRN-cheby.';
    k            = k+1;

    %% Precond. Chebyshev
    [x,xx,~, flopc] = chebyshev_pre_iko(A,b,lam,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'Chebyshev-precond.';
    k            = k+1;
    
    %% Acc IHS exact
    [~,xx,~, flopc] = acc_ihs_exact_iko(A,b,lam,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'Acc. IHS';
    k            = k+1;
    
    %% M-IHS exact
    params1      = struct('k0', k0);
    [~,xx,~, flopc] = mihs_exact_iko(A,b,lam,m,x1,TOL,MAXIT, params1);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'M-IHS-exact';
    k            = k+1;
    
    %% M IHS inexact
    params2      = struct('k0', k0, 'subtol', 1e-1, 'submaxit', k0);
    [~,xx,~, flopc] = mihs_inexact_iko(A,b,lam,m,x1,TOL,MAXIT, params2);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    names{k}     = 'M-IHS-inexact';
    k            = k+1;
    
    %% PD M IHS exact
%     params3      = struct('k0', k0);
%     [~,xx,~, flopc] = pd_mihs_over_exact_iko(A,b,lam,[3*k0, 3*k0],x1,[TOL TOL2],[MAXIT MAXIT2], params3);
%     err(:,k)     = err_x0(xx);
%     flop(:,k)    = flopc;
%     names{k}     = 'PD-MIHS-exact';
%     k            = k+1;
    
    %% PD M IHS inexact
%     params4      = struct('k0', k0, 'subtol', 1e-1, 'submaxit', k0);
%     [~,xx,~, flopc] = pd_mihs_over_inexact_iko(A,b,lam,[3*k0, 3*k0],x1,[TOL TOL2],[MAXIT MAXIT2], params4);
%     err(:,k)     = err_x0(xx);
%     flop(:,k)    = flopc;
%     names{k}     = 'PD-MIHS-inexact';
%     k            = k+1;
    
    %% mc
    err_mc(:,:,mc)  = err;
    flop_mc(:,:,mc) = flop;
end
err_m   = mean(err_mc,3);
flop_m  = mean(flop_mc,3);
err_s   = std(err_mc,0,3);
flop_s  = std(flop_mc,0,3);   
figure;
plot(log10(flop_m), log10(err_m));
legend(names{:})
