clear all
close all
clc
%%% THIS MAIN FILES TAKES SVD
%% SPEC
%size
n       = 30;
LEVEL   = 0.001;%repelem([0.001, 0.005, 0.01 0.04 0.08 0.1 0.12 0.16 0.2], 5);  
N_noi = numel(LEVEL);
EXAMPLE = 6;
PROBLEM = 4;

Ns      = 5;

hy_iter = 300;
MAXIT   = 25;
PTOL    = 0;
TOL     = 1e-2;
%% GENERATE DATA
switch PROBLEM
    case 1
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 2, EXAMPLE, LEVEL, true, true, 'type', 1, 'level', 10, 'p', 3*n);
    case 2
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 2, EXAMPLE, LEVEL, true, true, 'type', 2, 'level', 10, 'p', 10*n);
    case 3
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 2, EXAMPLE, LEVEL, true, true, 'type', 2, 'level', 1, 'p', 10*n);
    case 4
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 3, EXAMPLE, LEVEL, true, true, 'type', 1);
    case 5
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 3, EXAMPLE, LEVEL, true, true, 'angles', 0:10:179, 'p', 30*n);
    case 6
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
           = generate_data_IRtool_iko(n, 3, EXAMPLE, LEVEL, true, true, 'angles', 0:10:179, 'p', 40*n, 'd', 15*n);
    case 7
        [A, b, x1, x0, err_x0, dev, prob_info, err_x, k0, par0, xxOR, U, sigA, V] ...
            = generate_data_IRtool_iko(n, 3, EXAMPLE, LEVEL, true, true, 'angles', 0:10:179, 'p', 30*n);
end
[n,d] = size(A);
m     = 3*d;

%% pre allcoation for noise mc
N_met           = 10; %num of methods in the for bloop below
iter_noise      = zeros(N_met, N_noi);
err_noise       = zeros(N_met, N_noi);
par_noise       = zeros(N_met, N_noi);
dev_noise       = zeros(N_met, N_noi);
err_iter_noise  = cell(N_noi,1);
par_iter_noise  = cell(N_noi,1);
Utb             = U'*b;
for i = 1:N_noi
    fprintf('Noise mc: %d/%d\n', i, N_noi)
    err     = zeros(10,1);
    pare    = zeros(10,1);
    deve    = zeros(10,1);
    iter    = ones(10,1);
    err_iter= cell(10,1);
    par_iter= cell(10,1);
    k       = 1;
    
    %% 1 LS GCV
    [x, part, devt]  = LS_gcv_iko(U, sigA,V,b(:,i), x1);
    err(k)           = err_x(x,i);
    pare(k)          = part;
    deve(k)          = devt;
    err_iter{k}      = err(k)*ones(hy_iter,1);
    par_iter{k}      = pare(k)*ones(hy_iter,1);
    k                = k+1;
    
    %% 2 LS GCV eff
    [x, part, devt]  = LS_gcv_iko(speye(d), sigA,V, Utb(:,i), x1);
    err(k)           = err_x(x,i);
    pare(k)          = part;
    deve(k)          = devt;
    err_iter{k}      = err(k)*ones(hy_iter,1);
    par_iter{k}      = pare(k)*ones(hy_iter,1);
    k                = k+1;
    
    %% 3 HYBRID
    options          = HyBRset_iko('RegPar', 'WGCV', 'Iter', hy_iter, 'Omega', 'adapt', 'Reorth', 'on', 'x_true', xxOR{1}(:,i));
    [~, out]         = HyBR_iko(A, b(:,i), [], options);
    
    %STOP
    iter(k)          = out.iterations;
    err(k)           = out.Enrm(out.iterations);
    pare(k)          = out.alpha;
    err_iter{k}      = out.Enrm(1:out.iterations);
    par_iter{k}      = out.Alpha(1:out.iterations);
    k                = k+1;
    
    %NON-STOP
    iter(k)          = hy_iter;
    err(k)           = out.Enrm(end);
    pare(k)          = out.alpha(end);
    err_iter{k}      = out.Enrm;
    par_iter{k}      = out.alpha;
    k                = k+1;
    
    %% 5 HYBRID NON-STOP NON ADAPTIVE
    options2         = HyBRset_iko('RegPar', 'WGCV', 'Iter', hy_iter, 'Omega', 'corrected', 'Reorth', 'on', 'x_true', xxOR{1}(:,i));
    [~, out]         = HyBR_iko(A, b(:,i), [], options2);
    
    %STOP
    iter(k)          = out.iterations;
    err(k)           = out.Enrm(out.iterations);
    pare(k)          = out.alpha;
    err_iter{k}      = out.Enrm(1:out.iterations);
    par_iter{k}      = out.Alpha(1:out.iterations);
    k                = k+1;
    
    %NON-STOP
    iter(k)          = hy_iter;
    err(k)           = out.Enrm(end);
    pare(k)          = out.alpha(end);
    err_iter{k}      = out.Enrm;
    par_iter{k}      = out.alpha;
    k                = k+1;
    
    %% 7 LS SGCV + MIHS
    devs             = zeros(Ns,1);
    pars             = zeros(Ns,1);
    errs             = zeros(MAXIT,Ns);
    for j  =1:Ns
        SA                   = generate_SA_iko([A b(:,i)],m, false);
        Sb                   = SA(:,end);
        SA                   = SA(:,1:end-1);
        [Us,sigs,Vs]         = dsvd(SA);
        [~, pars(j), ~, devs(j)]  = LS_gcv_iko(Us, sigs,Vs,Sb, x1);
        
        %lower bound correction
        pars(j)  = pars(j)*m/n;
        devs(j)  = devs(j)*sqrt(m/n);
        
        %solution
        params   = struct('SA', SA);
        [~,xx]   = mihs_inexact_iko(A,b(:,i),pars(j),m,x1,TOL,MAXIT,params);
        errs(:,j)= [err_x(xx,i), nan*ones(1,MAXIT-size(xx,2))];
    end
    errs             = mean(errs,2); errs(isnan(errs))= [];
    pars             = mean(pars);
    devs             = mean(devs);
    iter(k)          = nnz(errs);
    err(k)           = errs(end);
    pare(k)          = pars;
    deve(k)          = devs;
    err_iter{k}      = errs;
    par_iter{k}      = pars*ones(hy_iter,1);
    k                = k+1;
    
    %% 8 Reg MIHS SVD + SGCV
    pars = zeros(MAXIT, Ns);
    errs = zeros(MAXIT, Ns);
    devs = zeros(Ns,1);
    for j = 1:Ns
        [~, xx, part, devt]  = reg_mihs_svd_lower_iko(A,b(:,i),m,x1,TOL,PTOL,MAXIT);
        pars(:,j)            = part;
        errs(:,j)            = err_x(xx,i);
        devs(j)              = devt;
    end
    errs             = mean(errs,2); errs(isnan(errs))= [];
    pars             = mean(pars,2); pars(isnan(pars))= [];
    devs             = mean(devs);
    iter(k)          = length(errs);
    err(k)           = errs(end);
    pare(k)          = pars(end);
    deve(k)          = devs;
    err_iter{k}      = errs;
    par_iter{k}      = pars;
    k                = k+1;
    
    %% 9 Reg MIHS SVD + SGCV - Corrected
    pars = zeros(MAXIT, Ns);
    errs = zeros(MAXIT, Ns);
    devs = zeros(Ns,1);
    for j = 1:Ns
        [~, xx, part, devt]  = reg_mihs_svd_corrected_iko(A,b(:,i),m,x1,TOL,PTOL,MAXIT);
        pars(:,j)            = part;
        errs(:,j)            = err_x(xx,i);
        devs(j)              = devt;
    end
    errs             = mean(errs,2); errs(isnan(errs))= [];
    pars             = mean(pars,2); pars(isnan(pars))= [];
    devs             = mean(devs);
    iter(k)          = length(errs);
    err(k)           = errs(end);
    pare(k)          = pars(end);
    deve(k)          = devs;
    err_iter{k}      = errs;
    par_iter{k}      = pars;
    k                = k+1;
    
    
    %% 10 Reg M-IHS GKL + SGCV
    pars = zeros(MAXIT, Ns);
    errs = zeros(MAXIT, Ns);
    devs = zeros(Ns,1);
    params = struct('L', min([k0(i)+100, n, d, m]));
    for j = 1:Ns
        [~, xx, part, devt]  = reg_mihs_gkl_lower_iko(A,b(:,i),m,x1,TOL,PTOL,MAXIT, params);
        pars(:,j)            = part;
        errs(:,j)            = err_x(xx,i);
        devs(j)              = devt;
    end
    errs             = mean(errs,2); errs(isnan(errs))= [];
    pars             = mean(pars,2); pars(isnan(pars))= [];
    devs             = mean(devs);
    iter(k)          = length(errs);
    err(k)           = errs(end);
    pare(k)          = pars(end);
    deve(k)          = devs;
    err_iter{k}      = errs;
    par_iter{k}      = pars;
    k                = k+1;
    
    %% Save data
    iter_noise(:,i)  = iter;
    err_noise(:,i)   = err;
    par_noise(:,i)  = pare;
    dev_noise(:,i)  = deve;
    err_iter_noise{i}= err_iter;
    par_iter_noise{i}= par_iter;
end


%% Statistics



%% Plot





