clear all
close all
clc

%%
n       = 2^16;
d       = 4000;
m       = 4000;


TOL     = 0;
MAXIT   = 100;


PROFILE = 5;%5,3
MAX     = 1;
KAPPA   = 8;%7,2
MC      = 32;

%% GENERATE DATA
[A, b, x0, x1, dev, U, sigA, V] = generate_data(n,d,...
    'correlated',               true,...
    'matrix correlation',       0.9,...
    'matrix deviation',         5,...
    'hansen singular values',   PROFILE,...
    'hansen max',               MAX,...
    'kappa',                    KAPPA,...
    'noise level',              0,...
    'hansen input',             false,...
    'input dev',                1,...
    'prior dev',                0,...
    'prior mean',               0,...
    'structure',                false,...
    'plot',                     false);

err_x0  = @(xx)(sqrt(sum((x0 - xx).^2, 1))/norm(x0));
%%
N_met   = 6;
err_mc  = zeros(MAXIT,N_met,MC);
flop_mc = zeros(MAXIT,N_met,MC);
for mc = 1:MC
    fprintf('mc: %d\n', mc);
    err   = zeros(MAXIT,N_met);
    flop  = zeros(MAXIT,N_met);
    k     = 1; 
    %% blendenpik
    [~,xx,~, flopc] = blendenpik_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
    %% lsrn + lsqr
    [~,xx,~, flopc] = lsrn_lsqr_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
    %% lsrn + chebyshev
    [~,xx,~, flopc] = lsrn_cheby_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
    %% Precond. Chebyshev
    [x,xx,~, flopc] = chebyshev_pre_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
    %% Acc IHS exact
    [~,xx,~, flopc] = acc_ihs_exact_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
    %% M-IHS exact
    [~,xx,~, flopc] = mihs_exact_iko(A,b,0,m,x1,TOL,MAXIT);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    k            = k+1;
    
  
    %% mc
    err_mc(:,:,mc)  = err;
    flop_mc(:,:,mc) = flop;
end
err_m   = mean(err_mc,3);
flop_m  = mean(flop_mc,3);
err_s   = std(err_mc,0,3);
flop_s  = std(flop_mc,0,3);

%%
figure;
plot(log10(flop_m), log10(err_m));
legend('Blendenpik', 'LSRN-lsqr', 'LSRN-cheby.', 'Chebyshev-precond.', 'Acc. IHS', 'M-IHS')
xlabel('flop count');
ylabel('log$_{10}\ \|x^i - x^*\|_{2}\big/\|x^*\|_{2}$')

