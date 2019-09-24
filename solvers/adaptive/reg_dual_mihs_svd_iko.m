function [x, xx, pars, in_iter, time] = reg_dual_mihs_svd_iko(A,b,m,x1,xtol,ptol,maxit,params)
%%REG_MIHS_SVD_CORRECTED adaptive regularizaiton of dual m-ihs algorithm
%- svd is used
%- NOT READY: lower parameter is set by gcv of (SA,Sb)
%
% [x, xx, pars, in_iter, time] = reg_dual_mihs_svd_iko(A,b,m,x1,xtol,ptol,maxit,params)
%
%   params.SAt   = sketch matrix
%   NO NEED for now : params.Sb   = sketch mea.
%
% I. Kurban Özaslan
% Bilkent EE
% September 2019
%


%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time]   = generate_SA_mihs(A', m, false);
else
    if(~isfield(params, 'SAt'))
        [SAt, rp_time] = generate_SA_mihs(A',m, false);
    else
        SAt      = params.SAt;
        rp_time = 0;
    end
end

%% some spec
[n,d]   = size(A);
pars    = zeros(maxit,1);
xx      = zeros(d,maxit);
in_iter = zeros(maxit,1);


%% SVD
[~, sig, V] = dsvd(SAt);
sig2       = sig.^2;
gamma       = @(lam)(1./(sig2 + lam));

%% initialization
x       = x1;
% xp      = x1*0;
nu      = zeros(n,1);
nup     = zeros(size(nu));


%% MAIN ITERATIONS
i       = 0;
EXIT    = false;
par_log = log10(sig(end));
par_p   = sig(end);
while(~EXIT)
    i       = i +1;
    %gradient
    grad    = b - A*x;
    
    %gcv vectors
    Vg      = V'*grad;
    Vnu     = V'*nu;
    f_i     = Vg + sig2.*Vnu;
    %minimization
    options             = optimset('TolX', 1e-3, 'MaxFunEvals', 15, 'Display','off');
    [par_log,~, ~, out] = fminbnd(@(lam)gcv(10^lam, sig2, f_i), ...
        par_log-1, par_log+5, options);    
    in_iter(i)          = out.funcCount;

    %parameter in decima
    par     = 10^par_log;
    g_par   = gamma(par);
    
    %solution by gcv
    dnu      = V*(g_par.*(Vg - par*Vnu));
    
    %momentum
    eff_rank= n - par*sum(g_par);
    r       = eff_rank/m;
    alpha   = (1-r)^2;
    beta    = r;
    
    %solution to iterate
    nun      = nu + alpha*dnu + beta*(nu - nup);
    
    %save
    nup      = nu;
    nu       = nun;
    x        = A'*nu;
    xx(:,i)  = x;
    pars(i)  = par;
    
    %exit flag
    XFLAG = norm(nu - nup)/norm(nup) <= xtol;
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

%%

function gcv_val = gcv(lam, sig2, f)

beta   = lam./(sig2+lam);
gcv_val    = norm(beta.*f)/sum(beta);

      
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