function [x, xx, pars, dev_est, in_iter, time] = reg_mihs_svd_lower_iko(A,b,m,x1,xtol,ptol,maxit,params)
%%REG_MIHS_SVD_LOWER adaptive regularizaiton of m-ihs algorithm
%- svd is used
%- lower parameter is set by gcv of (SA,Sb)
%
% [x, xx, pars, dev_est, in_iter, time] = reg_mihs_svd_lower(A,b,m,x1,xtol,ptol,maxit,params)
%
%   params.SA   = sketch matrix
%   params.Sb   = sketch mea.
%   
%
% I. Kurban Özaslan
% Bilkent EE
% September 2019



%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time]   = generate_SA_mihs([A b],m, false);
    Sb              = SA(:,end);
    SA              = SA(:,1:end-1);
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
end

%% some spec
[n,d]   = size(A);
pars    = zeros(maxit,1);
xx      = zeros(d,maxit);
in_iter = zeros(maxit,1);


%% SGCV
[U, sig, V]     = dsvd(full(SA));
[par_low, dev]  = LS_sgcv_lower_iko(U, sig,Sb);

%% some functions
sigi            = sig.^-1;
sig2            = sig.^2;
gamma           = @(lam)(1./(sig2+lam));

%% scale
par_low         = par_low*m/n;
par_low_log     = log10(par_low);
dev_est         = dev*sqrt(m/n);

%% Start MAIN ITERATION
x       = x1;
xp      = x1*0;

%% MAIN ITERATIONS
i       = 0;
par_p   = sig(1);
EXIT    = false;
par_log = par_low_log+2;
while(~EXIT)
    i       = i +1;
    %gradient
    grad    = A'*(b - A*x);
    
    %gcv vectors
    Vg      = V'*grad;
    Vx      = V'*x;
    f_gcv   = sigi.*Vg + sig.*Vx; 

    %minimization
    options             = optimset('TolX', 1e-3, 'MaxFunEvals', 15, 'Display','off');
    [par_log,~, ~, out] = fminbnd(@(lam)gcv_f(10^lam,sig2,f_gcv), ...
        max(par_log-2, par_low_log), par_log+1, options);    
    in_iter(i)          = out.funcCount;
    
    %parameter in decima
    par     = 10^par_log;
    
    %check lower bound
%     if(par > par_p*10)
%         par         = 0.5*(par_p + par_low);
%         par_log     = log10(par);
%     end
    g_par   = gamma(par);
    
    %solution by gcv
    dx      = V*(g_par.*(Vg - par*Vx));
    
    %momentum
    eff_rank= d - par*sum(g_par);
    r       = eff_rank/m;
    alpha   = (1-r)^2;
    beta    = r;
    
    %solution to iterate
    xn      = x + alpha*dx + beta*(x - xp);
    
    %save
    xp      = x;
    x       = xn;
    xx(:,i) = x;
    pars(i) = par;
    
    %exit flag
    XFLAG = norm(x - xp)/norm(xp) <= xtol;
    KFLAG = i >= maxit; 
    PFLAG = abs(par - par_p)/par_p < ptol;
   
    EXIT  = XFLAG || KFLAG || PFLAG;
    par_p = par;
end
% xx            = xx(:,1:i);
% pars          = pars(1:i);
% in_iter       = in_iter(1:i);
xx(:,i+1:end)   = nan;
pars(i+1:end)   = nan;
in_iter(i+1:end)= nan;
time    = toc + rp_time;

end

%%

function gcv = gcv_f(lam, sig2, f)

beta   = lam./(sig2+lam);
gcv    = norm(beta.*f)/sum(beta);

      
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
