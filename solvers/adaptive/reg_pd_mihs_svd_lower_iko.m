function [x, xx, pars1, pars2, dev_est, in_iter, time] = reg_pd_mihs_svd_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%%REG_PD_MIHS_SVD_LOWER adaptive regularizaiton of m-ihs algorithm
%- svd is used
%- lower parameter is set by gcv of (WASt,Wb)
%
% [x, xx, pars, dev_est, in_iter, time] = reg_pd_mihs_svd_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%
%   params.SAt   = sketch matrix
%   params.WASt  = sketch matrix
%   params.Wb   = sketch mea.
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
end

%% some spec
[n,d]   = size(A);
pars2   = zeros(maxit(1), maxit(2));
pars1   = zeros(maxit(1),1);
xx      = zeros(d,maxit(1));
in_iter = zeros(maxit(1), maxit(2));
tic;

%% SGCV
[U, sig, V]     = dsvd(full(WASt));
sigi            = sig.^-1;
sig2            = sig.^2;
gamma           = @(lam)(1./(sig2+lam));

%% lower bound
% [par_low, dev]  = LS_sgcv_lower_iko(U, sig,Wb);
% 
% %%scale
% par_low         = par_low*m(2)/n;
% par_low_log     = log10(par_low);
% dev_est         = dev*sqrt(m(2)/n);
dev_est = nan;
%% initialize
x       = x1;
nu      = zeros(n,1);
nup     = nu;
z       = zeros(m(1),1);
zp      = z;

%% MAIN ITERATIONS
i       = 0;
par_p1  = sig(end);
par_p2  = sig(1);
EXIT    = false;
par_low_log_t   = -16;
while(~EXIT)
    i       = i +1;
    %gradient
    b_i     = b - A*x;
    SAtnu   = SAt*nu;
    
    % take second dual and apply rp
    j       = 0;
    EXIT2   = false;
    par_log = log10(sig(1));
    while(~EXIT2)
        j       = j+1;
        
        %dual vectors
        Vg_ij   = V'*(SAt*(b_i - SAt'*z));
        Vz_ij   = V'*(z + SAtnu);
        
        %gcv vector and function
        f_ij    = sigi.*Vg_ij + sig.*Vz_ij;
        fun     = @(lam)( norm( gamma(lam).*f_ij )/sum( gamma(lam) ));
        
        %minimization
        options             = optimset('TolX', 1e-2, 'MaxFunEvals', 15, 'Display','off');%, 'PlotFcns', @optimplotx
        [par_log,~, ~, out] = fminbnd(@(lam)fun(10^lam), max(par_log-2, par_low_log_t), par_log+1, options);  
        in_iter(i,j)        = out.funcCount;


        %parameter in decima
        par     = 10^par_log;
    
        %check lower bound
%         if(par > par_p2*10)
%             par         = 0.5*(par_p2 + par_low);
%             par_log     = log10(par);
%         end
        gamma_par   = gamma(par);
        
        %solution of gcv
        dz          = V*(gamma_par.*(Vg_ij - par*Vz_ij));
    
        %momentum
        eff_rank    = m(1) - par*sum(gamma_par);
        r           = eff_rank/m(2);
        alpha       = (1-r)^2;
        beta        = r;
        
        %solution with momentum
        zn          = z + alpha*dz + beta*(z - zp);
        
        %update
        zp          = z;
        z           = zn;
        pars2(i,j)   = par;
        
        %termination criteria
        XFLAG = norm(z - zp)/norm(zp) <= xtol(2);
        KFLAG = j >= maxit(2); 
        PFLAG = abs(par - par_p2)/par_p2 < ptol(2);

        EXIT2 = XFLAG || KFLAG || PFLAG;
        par_p2 = par;
    end
    %dual measurement corrected by selected lambda
    b_i       = b_i - par*nu;

    %recover dual problem
    dnu     = (b_i - SAt'*z)/par;
    
    %momentum
    r       = eff_rank/m(1);
    alpha   = (1-r)^2;
    beta    = r;
        
    %solution
    nun         = nu + alpha*dnu + beta*(nu - nup);
    
    %update
    nup         = nu;
    nu          = nun; 
    
    x       = A'*nu;
    xx(:,i) = x;
    pars1(i)= par;
    %exit flag
    XFLAG = norm(nu - nup)/norm(nup) <= xtol(1);
    KFLAG = i >= maxit(1); 
    PFLAG = abs(par - par_p1)/par_p1 < ptol(1);
   
    EXIT  = XFLAG || KFLAG || PFLAG;
    par_p1= par;
    par_low_log_t = max(par_low_log_t, par_log);
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
