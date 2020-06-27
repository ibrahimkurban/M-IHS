function [x,xx,time,flopc,params] = dual_mihs_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%%DUAL_MIHS solver for underdetermined systems
%
% [x,xx,time,flopc,params] = dual_mihs_inexact_iko(A,b,lam,m,tol,x1,maxit,params)
%
%   params.SAt is sketch matrix
%   params.submaxit sets number max iteration for subsolver (25)
%   params.subtol   sets tolerance for subsolver (relative residual)(1e-2)
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University
%   MSc in EEE Dept.
%   November 2019
%
time    = zeros(maxit+3,1);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time,f_rp]   = generate_SA_mihs(A',m, false);
    noparam = true;
else
    noparam = false;
    if(~isfield(params, 'SAt'))
        [SAt, rp_time,f_rp] = generate_SA_mihs(A',m, false);
    else
        rp_time = 0;
        f_rp    = 0;
    end
end
time(1) = rp_time;
%% inexact tolerance
if(noparam||~isfield(params, 'subtol'))
    params.subtol = 1e-2;
end
if(noparam||~isfield(params, 'submaxit'))
    params.submaxit = 25;
end

%% data
[n,d]   = size(A);
xx      = zeros(d, maxit);
flopc   = zeros(maxit,1);
in_iter = zeros(maxit,1);

%% effective rank
if(noparam || ~isfield(params, 'k0'))
    [k0, ~, f_tr, tr_time] = hutchinson_estimator_iko(SAt, lam,2, params.subtol, params.submaxit);
    params.k0       = k0;
else
    k0      = params.k0;
    f_tr    = 0;
    tr_time = 0;
end
time(3) = tr_time;
%% momentum weights
r       = k0/m;

%%robust ones
% Ksup    = 1/(1-sqrt(r))^2 + 0.28;
% Kinf    = 1/(1+sqrt(r))^2 - 0.030;
% %
% alpha   = 4/( sqrt(Ksup) + sqrt(Kinf) )^2;
% beta    = (  ( sqrt(Ksup) - sqrt(Kinf) )/( sqrt(Ksup) + sqrt(Kinf) )  )^2;
%%theoric ones
alpha   = (1-r)^2;
beta    = r;

%% initialize
nu      = zeros(n,1);
nup     = zeros(size(nu));
dnu     = nup;
k = 0;
while(k < 2 || (norm(nu - nup)/norm(nup) >= tol && k < maxit))
    k           = k+1; tic;
    %grad
    grad        = b - A*(A'*nu) - lam*nu;
    
    %solve small system
    [dnu, in_iter(k), f_dx] = AA_b_solver_iko(SAt,grad, lam, params.subtol, params.submaxit, dnu);
    nun                     = nu + alpha*dnu + beta*(nu - nup);
    
    %update
    nup         = nu;
    nu          = nun;
    
    %store and count flop
    xx(:,k)     = A'*nu;
    flopc(k)    = f_dx(end) + 4*n*d + 8*n;
    time(k+3)   = toc;
end
xx  = xx(:,1:k);
x   = xx(:,end);
%% complexity (flop count refer to lightspeed malab packet)

% flop count
f_iter = 2*(n^2) + 4*n*d + 22*n;
flopc   = cumsum(flopc(:)') + [1:k]*f_iter + f_rp + f_tr;

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