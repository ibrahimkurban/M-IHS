function [x,xx,time,flopc,in_iter] = pd_mihs_v2_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%PD_MIHS_v2_EXACT_IKO solver for square-under cases
%
%  [xx, time] = pd_mihs_v2_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SAt   1. sketch amtrix
%   params.WASt  2. sketch matrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%
%   m,tol, maxit are all two length array, e.g. m = [m1,m2]
%

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time1,f_rp1]    = generate_SA_mihs(A',m(1), false);
    [WASt, rp_time2,f_rp2]   = generate_SA_mihs(SAt',m(2), false);
else
    if(~isfield(params, 'SAt'))
        [SAt, rp_time1,f_rp1] = generate_SA_mihs(A',m(1), false);
    end
    if(~isfield(params, 'WASt'))
        [WASt, rp_time2,f_rp2] = generate_SA_mihs(SAt',m(2), false);
    end
end
tic;

%% QR decomposition
[n,d]   = size(A);
if(lam == 0)
    [~, R] = qr(WASt,0);
    f_qr   = ceil(2*m(2)*m(1)^2-2/3*m(1)^3);
else
    [~, R] = qr([WASt;sqrt(lam)*speye(m(1))],0);
    f_qr   = ceil(2*(m(2)+m(1))*m(1)^2-2/3*m(1)^3);
end

%% momentum
r       = params.k0./[m1, m2];
alpha   = (1-r).^2;
beta    = r;

%% initialization
xx      = zeros(d,maxit(1));
in_iter = zeros(1,maxit(1));
nu      = zeros(n,1);
nup     = nu;
z       = zeros(m(1),1);
zp      = z;
i       = 0;
while(i < 2 || (norm(nu - nup)/norm(nup) >= tol(1) && i < maxit(1)))
    i      = i+1;
    bi     = b - A*(A'*nu) - lam*nu;
    
    %take dual and apply RP
    j      = 0;
    while(j < 2 || (norm(z - zp)/norm(zp) >= tol(2) && j < maxit(2)))
        j   = j+1;
        g   = SAt*(bi - SAt'*z) - lam*z;
        dz  = R\(R'\g);
        zn  = z + alpha(2)*dz + beta(2)*(z - zp);
        
        zp  = z;
        z   = zn;
    end
    in_iter(i)  = j;
    %convert
    dnu         = bi - SAt'*z;
    %momentum
    nun         = nu + alpha(1)*dnu +beta(1)*(nu - nup);
    %update
    nup         = nu;
    nu          = nun;
    xx(:,i)     = A'*nu;
end
xx      = xx(:,1:i);
x       = xx(:,end);
in_iter = in_iter(1:i);
%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time1 + rp_time2;

% flop count
f_iter1 = 4*n*d + 2*m(1)*n + 9*n;
f_iter2 = 4*m(1)*n + 2*m(1)^2 + 21*m(1) + n;
flopc   = [1:k]*f_iter1 + cumsum(in_iter)*f_iter2 + f_rp1 + f_rp2 + f_qr;
end




%% IF SA is not provided
function [SA, time, flopc] = generate_SA_mihs(A,SSIZE,wrep)
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
idx     = randsample(nt, SSIZE*N, wrep);                 % sampling pattern
SA      = HDA(idx, :)*(sqrt(nt)/sqrt(SSIZE));              % subsampling
time    = toc;

%% flop count refer to lightspeed malab packet
if(nargout > 2)
    f_DA    = 18*n + n*d;
    f_HDA   = ceil(nt*d*log2(7*m));
    f_SA    = 18*nt + m*d;
    
    flopc   = f_DA + f_HDA + f_SA;
end
end