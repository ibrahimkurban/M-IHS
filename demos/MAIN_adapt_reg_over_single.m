clear all
close all
clc
%%% THIS MAIN FILES TAKES SVD
%% SPEC
%size
n       = 64;
% d       = 2^9;

LEVEL   = 0.01;
PROFILE = 3;
EXAMPLE = 3;

HYBRID          = 100;

%% GENERATE DATA

[A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
           = generate_data_IRtool_iko(n, 3, EXAMPLE, LEVEL, 1, 1, 'angles', 0:1:179, 'p', 2*n, 'd', 5*n);

[n,d] = size(A);
%% Full LS Soution
[x_gcv, par_gcv, dev_est]  = LS_gcv_iko(U, sigA,V,b, x1);
[x_gcv2, par_gcv2, dev_est2]  = LS_gcv_iko(speye(d), sigA,V,U'*b, x1);
% [~, par_lc, ~, x_lc]                            = L_curve_analysis(U, sigA, V, b, x1, logspace(-15,1, 100));

fprintf('Reg LS GCV reg. par            : %1.2e\n', par_gcv);
fprintf('Reg LS GCV reg. par(eff) : %1.2e\n', par_gcv2);
% fprintf('Reg LS Lcurve reg. par         : %1.2e\n', par_lc);


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
fprintf('RP runs...\n')
m           = 3*d;
% [SA, rpt1]  = generate_SA_count([A, b, b - A*x0],m);
[SA, rpt1]  = generate_SA_iko([A, b, b - A*x0],m);
SA          = full(SA);
Sw          = SA(:,end);
Sb          = SA(:,end-1);
SA          = SA(:,1:end-2);


fprintf('Sketch Size 1            : %d\n', m)
fprintf('\n%2.2e sec elapsed\n\n', rpt1);


%% Classical Sketching Lower Bound
[Us, sigs, Vs]  = dsvd(SA);
sc              = sqrt(n/m);
[x_sgcv, par_sgcvp, dev_sest]  = LS_gcv_iko(Us, sigs,Vs,Sb, x1);
par_sgcv        = par_sgcvp/sc^2;
dev_est3        = dev_sest/sc;

%% METHODS
XTOL    = 1e-2;
MAXIT   = 25+ ceil(log(max(size(A))));
PTOL    = 0;
L       = d;%min([k0+300, d, n, m]); %L<d creates problem on lower bound not resolved

[~, xx1, pari1, dev_est4]  = reg_mihs_svd_lower_iko(A,b,m,x1,XTOL,PTOL,MAXIT);
[~, xx2, pari2, dev_est5]  = reg_mihs_gkl_lower_iko(A,b,m,x1,XTOL,PTOL,MAXIT,struct('L', L));


%% plot
figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('oracle error');
plot(log10(err_x(xxOR{2}))*ones(HYBRID,1), 'k:', 'linewidth', 2);
plot(log10(err_x(x_gcv))*ones(HYBRID,1), 'r:', 'linewidth', 2);
plot(log10(err_x(x_gcv2))*ones(HYBRID,1), 'g:', 'linewidth', 3);
plot(log10(err_x(x_sgcv))*ones(HYBRID,1), 'c:', 'linewidth', 2);
plot(log10(output.Enrm(1:end)), '-.', 'linewidth', 2);
plot(output.iterations, log10(output.Enrm(output.iterations)), 'o', 'markersize', 10)
plot(log10(output2.Enrm(1:end)), '-.', 'linewidth', 2);
plot(output2.iterations, log10(output2.Enrm(output2.iterations)), 'o', 'markersize', 10)
plot(log10(err_x(xx1)), 'linewidth', 2);
plot(log10(err_x(xx2)), 'linewidth', 2);

legend('Oracle Reg', 'LS GCV', 'LS GCV(eff)', 'CS lower bound', 'HYBRID', 'HYBRID stop by GCV', 'HYBRID corrected', 'HYBRID corrected stop by GCV', 'M-IHS (svd)', 'M-IHS (gkl)')

figure;
hold on; grid on;
xlabel('iterate, k'); ylabel('\lambda');
plot(log10(par0)*ones(HYBRID,1), 'k:', 'linewidth', 2);
plot(log10(par_gcv)*ones(HYBRID,1), 'r:', 'linewidth', 2);
plot(log10(par_gcv2)*ones(HYBRID,1), 'g:', 'linewidth', 3);
plot(log10(par_sgcv)*ones(HYBRID,1), 'c-.', 'linewidth', 3);
plot(log10(output.Alpha)*2, '-.', 'linewidth', 2)
plot(output.iterations, 2*log10(output.Alpha(output.iterations)), 'o', 'markersize', 10)
plot(log10(output2.Alpha)*2, '-.', 'linewidth', 2)
plot(output2.iterations, 2*log10(output2.Alpha(output2.iterations)), 'o', 'markersize', 10)
plot(log10(pari1(:, end)), 'linewidth', 2);
plot(log10(pari2(:,end)), 'linewidth', 2);
legend('Oracle Reg', 'LS GCV', 'LS GCV(eff)', 'CS lower bound', 'HYBRID', 'HYBRID stop by GCV', 'HYBRID corrected', 'HYBRID corrected stop by GCV', 'M-IHS (svd)', 'M-IHS (gkl)')

