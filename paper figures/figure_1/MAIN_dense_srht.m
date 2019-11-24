clear all
close all
clc

%%
n       = 2^15;
d       = 1000;
 

TOL     = 0;
MAXIT   = 20;

PROFILE = 5;%5,3
MAX     = 1e2;
KAPPA   = 8;%7,2
LEVEL   = 0.01;
MC      = 32;

% %% GENERATE DATA
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
N_met   = 2;
err_mc  = zeros(MAXIT,N_met,MC);
flop_mc = zeros(MAXIT,N_met,MC);
rr_mc   = zeros(MAXIT,N_met,MC);
in_iter_mc=zeros(MAXIT,MC);
names   = cell(N_met,1);

m = 5*k0;
%% Convergence rate vs Dimension for different probs
W = V*(sqrt((sigA.^2 + par0)).*V');
rate_x = @(xx)(sqrt(sum((W*(xx - xxOR{2})).^2)))/norm(W*xxOR{2});

%% MC
for mc = 1:MC
    fprintf('mc: %d\n', mc);
    err   = zeros(MAXIT,N_met);
    flop  = zeros(MAXIT,N_met);
    rr    = zeros(MAXIT,N_met);
    k     = 1;   
    
    %% M-IHS exact
    params1      = struct('k0', k0);
    [~,xx,~, flopc] = mihs_exact_iko(A,b,lam,m,x1,TOL,MAXIT, params1);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    rr(:,k)      = rate_x(xx);
    names{k}     = 'M-IHS-exact';
    k            = k+1;
    
    %% M IHS inexact
    params2      = struct('k0', k0, 'subtol', 1e-1, 'submaxit', k0);
    [~,xx,~, flopc, in_iter_mc(:,mc)] = mihs_inexact_iko(A,b,lam,m,x1,TOL,MAXIT, params2);
    err(:,k)     = err_x0(xx);
    flop(:,k)    = flopc;
    rr(:,k)      = rate_x(xx);
    names{k}     = 'M-IHS-inexact';
    k            = k+1;
    
    %% mc
    err_mc(:,:,mc)  = err;
    flop_mc(:,:,mc) = flop;
    rr_mc(:,:,mc)   = rr;
end
err_m   = mean(err_mc,3);
flop_m  = mean(flop_mc,3);
err_s   = std(err_mc,0,3);
flop_s  = std(flop_mc,0,3);  

figure; hold on;
pp1 = plot(log10(permute(rr_mc(:,1,:), [1 3 2])), 'b', 'linewidth', 0.5);
pp2 = plot(log10(permute(rr_mc(:,2,:), [1 3 2])), 'r', 'linewidth', 0.2);
p3 = plot(log10(sqrt(k0/m).^[1:MAXIT]), 'c^-', 'linewidth', 2, 'MarkerSize', 9);
grid minor
legend([p3, pp1(1), pp2(1)], 'theoretical', 'Exact M-IHS', 'Inexact M-IHS')
xlabel('iteation number, $i$'); ylabel('log$_{10}\ \|x^i - x^*\|_{D^{-1}}\big/\|x^*\|_{D^{-1}}$')

