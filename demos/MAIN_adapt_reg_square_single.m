clear all
close all
clc
%%% THIS MAIN FILES TAKES SVD
%% SPEC
n       = 64;

LEVEL   = 0.01;
PROFILE = 2;
EXAMPLE = 4;

hy_iter = 200;
SVD_A   = true;
PLOT_A  = true;
%% GENERATE DATA
[A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 2, 4, LEVEL, SVD_A, PLOT_A, 'type', 2);

% [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
%     = generate_data_IRtool_iko(n, 1, 7, LEVEL, SVD_A, PLOT_A, 'type', 1, 'level', 2);

[n,d]       = size(A);

%% Full LS Soution
[x_gcv, par_gcv, dev_est]  = LS_gcv_iko(U, sigA,V,b, x1);
fprintf('Reg LS GCV reg. par    : %1.2e\n', par_gcv);

%% Hybrid
options          = HyBRset_iko('RegPar', 'WGCV', 'Iter', hy_iter, 'Omega', 'adapt', 'Reorth', 'on', 'x_true', xxOR{1});
fprintf('\nHyBR runs...'); tic;
[x_hy, output] = HyBR_iko(A, b, [], options);
fprintf('%2.1e sec elapsed\n\n', toc);

%% METHODS
XTOL    = [min(LEVEL,5e-3) 0];
MAXIT   = [15 25];
PTOL    = [0 0];
m1      = 2*k0;
m2      = 5*k0;
L       = min(k0+300, m1);

[~, xx1, pari1, pars1]      = reg_pd_mihs_gkl_lower_iko(  A,b,[m1, m2],x1,XTOL,PTOL, MAXIT, struct('L', L));
[~, xx2, pari2, pars2]      = reg_pd_mihs_svd_lower_iko(  A,b,[m1, m2],x1,XTOL,PTOL, MAXIT);
% [xx3, ~, pari3, in_iter]    = IS_square_svd_reg_fmin_log( A,b, x1, m1, m2, XTOL(1), PTOL(1), MAXIT(1), MAXIT(2));
% [xx4, ~, pari4]             = IS_square_gkl_v2_reg(       A,b, x1, m1, m2, XTOL(1), PTOL(1), MAXIT(1), MAXIT(2), L);

%% plot
figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('oracle error');
plot(log10(err_x(xxOR{2}))*ones(hy_iter,1), 'k:', 'linewidth', 2);
plot(log10(err_x(x_gcv))*ones(hy_iter,1), 'r:', 'linewidth', 2);
plot(log10(output.Enrm(1:end)), '-.', 'linewidth', 2);
plot(output.iterations, log10(output.Enrm(output.iterations)), 'o', 'markersize', 10)
plot(log10(err_x(xx1)), 'linewidth', 2);
plot(log10(err_x(xx2)), 'linewidth', 2);
% plot(log10(err_x(xx3)), 'x-', 'linewidth', 2);
% plot(log10(err_x(xx4)), 'o-', 'linewidth', 2);
legend('Oracle Reg', 'LS GCV', 'HYBRID', 'HYBRID stop by GCV', 'PD gkl', 'PD svd')

figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('\lambda');
plot(log10(par0)*ones(hy_iter,1), 'k:', 'linewidth', 2);
plot(log10(par_gcv)*ones(hy_iter,1), 'r:', 'linewidth', 2);
plot(log10(output.Alpha.^2), '-.', 'linewidth', 2)
plot(log10(pari1(:, end)), 'linewidth', 2);
plot(log10(pari2(:,end)), 'linewidth', 2);
% plot(log10(pari3(:,end)), 'x-', 'linewidth', 2);
% plot(log10(pari4(:,end)), 'o-', 'linewidth', 2);
legend('Oracle Reg', 'LS GCV', 'HYBRID','PD-gkl', 'PD svd')

