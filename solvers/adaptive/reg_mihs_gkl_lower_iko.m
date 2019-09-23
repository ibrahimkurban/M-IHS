function [x, xx, pars, dev_est, in_iter, time] = reg_mihs_gkl_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%%REG_MIHS_GKL_LOWER adaptive regularizaiton of m-ihs algorithm
%- gkl is used
%- lower parameter is set by gcv of (SA,Sb)
%
% [x, xx, pars, dev_est, in_iter, time] = reg_mihs_gkl_lower(A,b,m,x1,xtol,ptol,maxit,params)
%
%   params.SA   = sketch matrix
%   params.Sb   = sketch mea.
%   params.L    = gkl size
%
% I. Kurban Özaslan
% Bilkent EE
% September 2019
%


%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time]   = generate_SA_mihs([A b],m, false);
    Sb              = SA(:,end);
    SA              = SA(:,1:end-1);
    params.L        = min(size(A));
else
    if(~isfield(params, 'SA'))
        [SA, rp_time] = generate_SA_mihs([A b],m, false);
        Sb      = SA(:,end);
        SA      = SA(:,1:end-1);
    else
        SA      = params.SA;
        Sb      = params.Sb;
        rp_time = 0;
    end
    if(~isfield(params, 'L'))
        params.L= min(size(A));
    end
end

%% some spec
[n,d]   = size(A);
pars    = zeros(maxit,1);
xx      = zeros(d,maxit);
in_iter = zeros(maxit,1);
L       = params.L;

%% SGCV
ASSb    = SA'*Sb;
theta1  = norm(ASSb);
[R,V,D]         = solver_gkl_v2_iko_c(SA, ASSb, L);
L               = size(R,1);
RR              = R*R';
RR_I            = @(lam)(RR + lam*speye(L));
%% lower bound
[par_low, dev]  = LS_sgcv_lower_gkl_iko(R,Sb, theta1);

%%scale
par_low_log     = log10(par_low*m/n);
dev_est         = dev*sqrt(m/n);

%% Start MAIN ITERATION
x       = x1;
% xp      = x1*0;
y       = V'*x;
yp      = zeros(L,1);

%% MAIN ITERATIONS
i       = 0;
par_p   = par_low;
EXIT    = false;
par_log = par_low_log+2;
while(~EXIT)
    i       = i +1;
    %gradient
    grad    = A'*(b - A*x);

    %lanczos
    f_gcv   = D'*grad + R*y;
    fun     = @(lam)(norm(((RR_I(lam))\f_gcv)))/(solver_trace_bidiag_inv_iko(RR_I(lam)));
    
    %minimization
    options                 = optimset('TolX', 1e-3, 'MaxFunEvals', 15, 'Display','off');
    [par_log,~, ~, out]     = fminbnd(@(lam)fun(10^lam), max(par_log-2, par_low_log), par_log+1, options);
    in_iter(i)              = out.funcCount;
    %        [pari, ~, gcv, lambdas1] = grid_search_log(fun, lambdas);

    %parameter in decima
    par     = 10^par_log;
    
    %solution by gcv
    RR_Ip   = RR_I(par);
    g_temp  = D'*(grad - par*x);
    dy      = R'*(RR_Ip\g_temp);
    
    %momentum
    eff_rank= L - par*solver_trace_bidiag_inv_iko(RR_Ip);
    r       = eff_rank/m;
    alpha   = (1-r)^2;
    beta    = r;
    
    %solution to iterate
    yn      = y + alpha*dy + beta*(y - yp);
%     xn      = x + alpha*dx + beta*(x - xp);
    
    %save
    yp      = y;
    y       = yn;
    x       = V*y;
    xx(:,i) = x;
    pars(i) = par;
    
    %exit flag
    XFLAG = norm(y - yp)/norm(yp) <= xtol;
    KFLAG = i >= maxit; 
    PFLAG = abs(par - par_p)/par_p <= ptol;
   
    EXIT  = XFLAG || KFLAG || PFLAG;
    par_p = par;
end
% xx      = xx(:,1:i);
% pars    = pars(1:i);
% in_iter = in_iter(1:i);
xx(:,i+1:end)   = nan;
pars(i+1:end)   = nan;
in_iter(i+1:end)= nan;
time    = toc + rp_time;

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

