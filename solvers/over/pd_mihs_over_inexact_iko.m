function [x,xx, time,flopc,in_iter,params] = pd_mihs_over_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%%PD_MIHS_v1_EXACT_IKO solver for square-over cases
%
% [x,xx, time,flopc,in_iter,params] = pd_mihs_over_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SA   1. sketch amtrix
%   params.WASt  2. sketch matrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%
%   m,tol, maxit are all two length array, e.g. m = [m1,m2]
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
%   params.submaxit sets number max iteration for subsolver (25)
%   params.subtol   sets tolerance for subsolver (relative residual)(1e-2)
%
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%

time = zeros(maxit+3, 1);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time1,f_rp1]    = generate_SA_mihs(A,m(1), false);
    [WASt, rp_time2,f_rp2]  = generate_SA_mihs(SA',m(2), false);
    noparam                 = true;
else
    noparam                 = false;
    if(~isfield(params, 'SAt'))
        [SA, rp_time1,f_rp1] = generate_SA_mihs(A,m(1), false);
        [WASt, rp_time2,f_rp2] = generate_SA_mihs(SA',m(2), false);
    else
        SA          = params.SA;
        WASt        = params.WASt;
        rp_time1    = 0;
        rp_time2    = 0;
        f_rp1       = 0;
        f_rp2       = 0;
    end
end
time(1) = rp_time1+rp_time2;
%% data
[n,d]   = size(A);
xx      = zeros(d, maxit(1));
flopc   = zeros(1, maxit(1));
in_iter = zeros(1, maxit(1));

%% effective rank
if(noparam || ~isfield(params, 'k0'))
    [k0, ~, f_tr, tr_time] = hutchinson_estimator_iko(WASt, lam, 2, 1e-1, 50);
    params.k0   = k0;
else
    k0          = params.k0;
    f_tr        = 0;
    tr_time     = 0;
end
time(2) = 0;
time(3) = tr_time;
%% inexact tolerance
if(noparam || ~isfield(params, 'subtol'))
    params.subtol = 1e-2;
end
if(noparam|| ~isfield(params, 'submaxit'))
    params.submaxit = max(25, k0);
end

%% momentum
r       = k0./m;
alpha   = (1-r).^2;
beta    = r;

%% initialization
x       = x1;
xp      = x;
lami    = lam^-1;
nu      = zeros(m(1),1);
nup     = nu;
dnu     = nu;
i       = 0;
while(i < 2 || (norm(x - xp)/norm(xp) >= tol(1) && i < maxit(1)))
    i      = i+1; tic;
    bi     = A'*(b - A*x) - lam*x;
    
    %inexact subsolver
    j      = 0;
    while(j < 2 || (norm(nu - nup)/norm(nup) >= tol(2) && j < maxit(2)))
        j    = j+1;
        %sub prob gradient
        g    = SA*(bi - SA'*nu) - lam*nu;
        
        %sub solve problem
        [dnu, in_iter(i),f_dx]   = AA_b_solver_iko(WASt, g, lam, params.subtol, params.submaxit,dnu);
        
        %momentum
        nun  = nu + alpha(2)*dnu + beta(2)*(nu - nup);
        
        %update
        nup     = nu;
        nu      = nun;
        flopc(i)= flopc(i) + f_dx(end) + 4*m(1)*d + 7*m(1) + d; 
    end
    in_iter(i) = j;
    dx      = lami*(bi - (SA'*nu));
    xn      = x + alpha(1)*dx +beta(1)*(x - xp);
    xp      = x;
    x       = xn;
    xx(:,i) = x;
    time(i+3) = toc;
end
xx      = xx(:,1:i);
in_iter = in_iter(1:i);

%% complexity (flop count refer to lightspeed malab packet)
% flop count
f_iter  = 4*n*d + 2*m(1)*d + 9*d;
flopc   = [1:i]*f_iter + cumsum(flopc) + f_rp1 + f_rp2 + f_tr;
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
