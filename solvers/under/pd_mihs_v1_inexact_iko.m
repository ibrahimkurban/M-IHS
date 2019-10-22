function [x,xx,time,flopc,in_iter] = pd_mihs_v1_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%%PD_MIHS_v1_EXACT_IKO solver for square-under cases
%
%  [x,xx,time,flopc] = pd_mihs_v2_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SAt   1. sketch amtrix
%   params.WASt  2. sketch matrix
%
%   m,tol, maxit are all two length array, e.g. m = [m1,m2]
%
%   params.submaxit sets number max iteration for subsolver (25)
%   params.subtol   sets tolerance for subsolver (relative residual)(1e-2)

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time1,f_rp1]    = generate_SA_mihs(A',m(1), false);
    [WASt, rp_time2,f_rp2]   = generate_SA_mihs(SAt',m(2), false);
    noparam                  = true;
else
    noparam                  = false;
    if(~isfield(params, 'SAt'))
        [SAt, rp_time1,f_rp1] = generate_SA_mihs(A',m(1), false);
        [WASt, rp_time2,f_rp2] = generate_SA_mihs(SAt',m(2), false);
    else
        SAt     = params.SAt;
        WASt    = params.WASt;
        rp_time1= 0;
        rp_time2= 0;
        f_rp1   = 0;
        f_rp2   = 0;
    end
end

%% effective rank
if(noparam || ~isfield(params, 'k0'))
    [k0, ~, f_tr] = hutchinson_estimator_iko(WASt, WASt, lam, 2, 1e-1, 100);
    params.k0   = k0;
else
    k0          = params.k0;
    f_tr        = 0;
end

%% inexact tolerance
if(noparam || ~isfield(params, 'subtol'))
    params.subtol = 1e-2;
end
if(noparam|| ~isfield(params, 'submaxit'))
    params.submaxit = max(25, k0);
end

%% data
[n,d]   = size(A);
xx      = zeros(d, maxit(1));
flopc   = zeros(1, maxit(1));
in_iter = zeros(1, maxit(1));
tic;

%% momentum
r       = k0./m;
alpha   = (1-r).^2;
beta    = r;

%% initialization
x       = x1;
nu      = zeros(n,1);
nup     = nu;
z       = zeros(m(1),1);
zp      = z;
dz      = z;
i       = 0;
while(i < 2 || (norm(nu - nup)/norm(nup) >= tol(1) && i < maxit(1)))
    i      = i+1;
    bi     = b - A*x - lam*nu;
    
    %take dual and apply RP
    j      = 0;
    while(j < 2 || (norm(z - zp)/norm(zp) >= tol(2) && j < maxit(2)))
        j   = j+1;
        %subgrad
        g   = SAt*(bi - SAt'*z) - lam*z;
        
        %inexact subsolver
        [dz, in_iter(i),f_dx]   = AA_b_solver_iko(WASt, g, lam, params.subtol, params.submaxit,dz);
        zn                      = z + alpha(2)*dz + beta(2)*(z - zp);
        
        %update
        zp          = z;
        z           = zn;
        flopc(i)    = flopc(i) + f_dx(end) + 4*m(1)*n + 6*m(1)+n;
    end
    in_iter(i)  = j;
    %convert
    dnu         = (bi - SAt'*z)/lam;
    %momentum
    nun         = nu + alpha(1)*dnu +beta(1)*(nu - nup);
    %update
    nup         = nu;
    nu          = nun;
    
    x           = A'*nu;
    xx(:,i)     = x;
end
xx      = xx(:,1:i);
x       = xx(:,end);
in_iter = in_iter(1:i);

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time1 + rp_time2;

% flop count
f_iter1 = 4*n*d + 2*m(1)*n + 9*n;
flopc   = [1:i]*f_iter1 + cumsum(flopc) + f_rp1 + f_rp2 + f_tr;
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