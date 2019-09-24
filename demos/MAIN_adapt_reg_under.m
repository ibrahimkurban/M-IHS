clear all
close all
clc

%% SPEC
%size
n       = 2^9;
d       = 2^13;

LEVEL   = 0.0001;
KAPPA   = 6;
MAX     = 1;
PROFILE = 5;
HYBRID  = 150;


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
    'plot',                     true);
fprintf('noise deviation: %2.2e\n', dev)


%% ORACLE Regularized SOL
[xxOR, parOR, errOR, err_rid] = oracle_methods(U, sigA, V,b,x0,x1);
err_x   = @(xx)(errOR{1}(xx,1));
err_xreg= @(xx)(errOR{2}(xx,1));
k0      = parOR(1);
par     = parOR(2);
fprintf('effective rank : %d\n', k0);
fprintf('Oracle reg. par: %1.2e\n', par);
%% Full LS Soution
[x_gcv, par_gcv, dev_est1] = LS_gcv_iko(U, sigA,V,b, x1);
% [~, par_lc, ~, x_lc] = L_curve_analysis(U, sigA, V, b, x1, logspace(-15,1, 100));
% [reg_corner,rho,eta,reg_param] = l_curve(U,sigA,b,'Tikh',V,eye(d));
fprintf('Reg LS GCV reg. par    : %1.2e\n', par_gcv);
% fprintf('Reg LS Lcurve reg. par : %1.2e\n', par_lc);

%% Hybrid
options         = HyBRset_iko('RegPar', 'WGCV', 'Iter', HYBRID, 'Reorth', 'on', 'x_true', xxOR{1});
fprintf('\nHyBR runs...'); tic;
[x_hy, output] = HyBR_iko(A, b, [], options);
fprintf('%2.1e sec elapsed\n\n', toc);

options         = HyBRset_iko('RegPar', 'WGCV', 'Iter', HYBRID, 'Omega', 'corrected', 'Reorth', 'on', 'x_true', xxOR{1});
fprintf('\nHyBR runs...'); tic;
[x_hy2, output2] = HyBR_iko(A, b, [], options);
fprintf('%2.1e sec elapsed\n\n', toc);
%% RP
m       = 3*n;
% SAtt = SAt;
SAt     = generate_SA_iko(A',m);

%% methods
xtol        = 1e-2;
maxit       = 25;
ptol        = 0;
params.SAt  = SAt;
L           = k0+200;
params.L    = L;


% [xx3, ~, pari3]                 = IS_under_svd_reg(A, b, x1, m,xtol, maxit, SAt);
% [xx4, ~, pari4, in_iter4]       = IS_under_gkl_v2_reg(A, b, x1, m,xtol,maxit, L,SAt);
[~, xx1, pari1, in_iter1]       = reg_dual_mihs_svd_iko(A,b,m,x1,xtol,ptol,maxit,params);
[~, xx2, pari2, in_iter2]       = reg_dual_mihs_gkl_iko(A,b,m,x1,xtol,ptol,maxit,params);

%% plot
figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('oracle error');
plot(log10(err_rid)*ones(HYBRID,1), 'k:', 'linewidth', 2);
plot(log10(err_x(x_gcv))*ones(HYBRID,1), 'r:', 'linewidth', 2);
plot(log10(output.Enrm(1:end)), '-.', 'linewidth', 2);
plot(output.iterations, log10(output.Enrm(output.iterations)), 'o', 'markersize', 10)
plot(log10(output2.Enrm(1:end)), '-.', 'linewidth', 2);
plot(output2.iterations, log10(output2.Enrm(output2.iterations)), 'o', 'markersize', 10)
plot(log10(err_x(xx1)), 'linewidth', 2);
plot(log10(err_x(xx2)), 'linewidth', 2);
% plot(log10(err_x(xx3)), 'x-', 'linewidth', 2);
% plot(log10(err_x(xx4)), 'o-', 'linewidth', 2);
legend('Oracle Reg', 'LS GCV',  'Hybrid wgcv', 'GCV stop', 'Hybrid wgcv, corr.', 'Dual MIHS (svd)', 'Dual MIHS (gkl)')

figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('\lambda');
plot(log10(parOR(2))*ones(5,1), 'k:', 'linewidth', 2);
plot(log10(par_gcv)*ones(5,1), 'r:', 'linewidth', 2);
plot(log10(output.Alpha)*2, '-.', 'linewidth', 2)
plot(output.iterations, 2*log10(output.Alpha(output.iterations)), 'o', 'markersize', 10)
plot(log10(output2.Alpha)*2, '-.', 'linewidth', 2)
plot(output2.iterations, 2*log10(output2.Alpha(output2.iterations)), 'o', 'markersize', 10)
plot(log10(pari1), 'linewidth', 2);
plot(log10(pari2), 'linewidth', 2);
% plot(log10(pari3), 'x-', 'linewidth', 2);
% plot(log10(pari4), 'o-', 'linewidth', 2);
legend('Oracle Reg', 'LS GCV',  'Hybrid wgcv', 'GCV stop', 'Hybrid wgcv, corr.', 'gcv STOP corr.', 'Dual MIHS (svd)', 'Dual MIHS (gkl)')


