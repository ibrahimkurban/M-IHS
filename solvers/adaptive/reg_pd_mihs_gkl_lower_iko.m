function [x, xx, pars1, pars2, dev_est, in_iter, time] = reg_pd_mihs_gkl_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%%REG_PD_MIHS_GKL_LOWER adaptive regularizaiton of m-ihs algorithm
%-  gkl is used (upper bidiagonal factorization
%- lower parameter is set by gcv of (WASt,Wb)
%
% [x, xx, pars, dev_est, in_iter, time] = reg_pd_mihs_svd_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%
%   params.SAt   = sketch matrix
%   params.WASt  = sketch matrix
%   params.Wb    = sketch mea.
%   params.L     = size of gkl procedure
%   
%
% I. Kurban Özaslan
% Bilkent EE
% September 2019



%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time1]     = generate_SA_mihs(A',m(1), false);
    [WASt, rp_time2]    = generate_SA_mihs([SAt', b],m(2), false);
    Wb                  = WASt(:,end);
    WASt(:,end)         = [];
    params.L            = min(m);
else
    if(~isfield(params, 'SAt'))
        [SAt, rp_time1]     = generate_SA_mihs(A',m(1), false);
        [WASt, rp_time2]    = generate_SA_mihs([SAt', b],m(2), false);
        Wb                  = WASt(:,end);
        WASt(:,end)         = [];
    else
        SAt     = params.SAt;
        WASt    = params.WASt;
        Wb      = params.Wb;
        rp_time1= 0;
        rp_time2= 0;   
    end
    if(~isfield(params, 'L'))
        params.L = min(m);
    else
        params.L = min(params.L, min(m));
    end
end

%% some spec
[n,d]   = size(A);
pars2    = zeros(maxit(1), maxit(2));
pars1    = zeros(maxit(1), 1);
xx      = zeros(d,maxit(1));
in_iter = zeros(maxit(1), maxit(2));
tic;


%% SGCV
L               = params.L;
SAWWb           = SAt*b;
[R,V,D]         = solver_gkl_v2_iko_c(WASt, SAWWb, L);
L               = size(R,1);
RR              = R*R';
RR_I            = @(lam)(RR + lam*speye(L));
params.L        = L;

%% largest sing value estimate (greshgorin disc can be used instead)
sig1_log = log10(max(svd(full(R(1:10,1:10)))));

%% lower bound
% theta1          = norm(SAWWb);
% [par_low, dev]  = LS_sgcv_lower_gkl_iko(R,Wb, theta1);
% 
% %%scale
% par_low         = par_low*m(2)/n;
% par_low_log     = log10(par_low);
% dev_est         = dev*sqrt(m(2)/n);
dev_est = nan;

%% Trace estimator
lambdas     = linspace(-10, 4, 100);
tr_log_fun  = trace_estimator_fun_iko(lambdas, RR);

%% initialize
x       = x1;
nu      = zeros(n,1);
nup     = nu;
z       = zeros(m(1),1);
zp      = z;
% y       = zeros(L,1);
% yp      = y;


%% MAIN ITERATIONS
i       = 0;
par_p1  = 1e-10;
par_p2  = 1e5;
EXIT    = false;
par_low_log = -16;
while(~EXIT)
    i       = i +1;
    %gradient
    b_i     = b - A*x;
    SAtnu   = SAt*nu;
    
    % take second dual and apply rp
    j       = 0;
    EXIT2   = false;
    par_log = sig1_log;
    while(~EXIT2)
        j       = j+1;
        
        %dual vectors
        Dg_ij   = D'*(SAt*(b_i - SAt'*z));
        Dz_ij   = D'*(z + SAtnu);
        
        %gcv vector and function
        f_ij    = Dg_ij + RR*Dz_ij;
        fun     = @(lam)( norm(RR_I(lam)\f_ij) )/( tr_log_fun(log10(lam)) );
        
        %minimization
        options             = optimset('TolX', 1e-2, 'MaxFunEvals', 15, 'Display','off');%, 'PlotFcns', @optimplotx
        [par_log,~, ~, out] = fminbnd(@(lam)fun(10^lam), max(par_log-2, par_low_log), par_log+1, options);  
        in_iter(i,j)        = out.funcCount;
        

        %parameter in decima
        par     = 10^par_log;
    
        %check lower bound
%         if(par > par_p2*10)
%             par         = 0.5*(par_p2 + par_low);
%             par_log     = log10(par);
%         end
        RR_par   = RR_I(par);
        
        %solution of gcv
        dy          = R'*(RR_par\(Dg_ij - par*Dz_ij));
    
        %momentum
        eff_rank    = L - par*sum(tr_log_fun(par_log));
        r           = eff_rank/m(2);
        alpha       = (1-r)^2;
        beta        = r;
        
        %solution with momentum
        zn          = z + alpha*(V*dy) + beta*(z - zp);
        
        %update
        zp          = z;
        z           = zn;
        pars2(i,j)   = par;
        
        %krylov solution
%         z           = V*y;
        
        %termination criteria
        XFLAG = norm(z - zp)/norm(zp) <= xtol(2);
        KFLAG = j >= maxit(2); 
        PFLAG = abs(par - par_p2)/par_p2 < ptol(2);

        EXIT2 = XFLAG || KFLAG || PFLAG;
        par_p2 = par;
    end
    %dual measurement corrected by selected lambda
    b_i         = b_i - par*nu;

    %recover dual problem
    dnu         = (b_i - SAt'*z)/par;
    
    %momentum
    r           = eff_rank/m(1);
    alpha       = (1-r)^2;
    beta        = r;
        
    %solution
    nun         = nu + alpha*dnu + beta*(nu - nup);
    
    %update
    nup         = nu;
    nu          = nun; 
    
    x           = A'*nu;
    xx(:,i)     = x;
    pars1(i)    = par;
    %exit flag
    XFLAG = norm(nu - nup)/norm(nup) <= xtol(1);
    KFLAG = i >= maxit(1); 
    PFLAG = abs(par - par_p1)/par_p1 < ptol(1);
   
    EXIT  = XFLAG || KFLAG || PFLAG;
    par_p1= par;
    
    par_low_log = max(par_log, par_low_log);
end
% xx            = xx(:,1:i);
% pars          = pars(1:i);
% in_iter       = in_iter(1:i);
xx(:,i+1:end)   = nan;
pars1(i+1:end)  = nan;
pars2(i+1:end,:)   = nan;
in_iter(i+1:end,:)= nan;
time    = toc + rp_time1 + rp_time2;

end

%% IF SA is not provided
function [SA, time, flopc] = generate_SA_mihs(A,m,wrep)
%%GENERATE_SA generates ROS sketch matrix
%
%   [SA, time, flopc] = generate_SA_mihs(A,SSIZE,wrep)
%
%   E[SA'SA] = I

if(issparse(A))
    A = full(A);
end
if(~exist('wrep', 'var'))
    wrep = false;
end
[n,d]   = size(A);
nt      = ceil(n/1000)*1000; %BLENDENPIK sugg. for FFTW
tic;
radem   = (randi(2, n, 1) * 2 - 3);                     % rademacher
DA      = A .* radem;                                   % one half+1 and rest -1
HDA     = dct(DA,nt);                                      % DCT transform
idx     = randsample(nt, m, wrep);                 % sampling pattern
SA      = HDA(idx, :)*(sqrt(nt)/sqrt(m));              % subsampling
time    = toc;

%% flop count refer to lightspeed malab packet
if(nargout > 2)
    f_DA    = 18*n + n*d;
    f_HDA   = ceil(nt*d*log2(7*m));
    f_SA    = 18*nt + m*d; 
    
    flopc   = f_DA + f_HDA + f_SA;
end
end
