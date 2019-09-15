function [x,xx, time,flopc,in_iter] = pd_mihs_v1_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%PD_MIHS_v1_EXACT_IKO solver for square-over cases
%
% [x,xx, time,flopc,in_iter] = pd_mihs_v1_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SA   1. sketch amtrix
%   params.WASt  2. sketch matrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%
%   m,tol, maxit are all two length array, e.g. m = [m1,m2]
%

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time1,f_rp1]    = generate_SA_mihs(A,m(1), false);
    [WASt, rp_time2,f_rp2]   = generate_SA_mihs(SA',m(2), false);
else
    if(~isfield(params, 'SAt'))
        [SA, rp_time1,f_rp1] = generate_SA_mihs(A,m(1), false);
        [WASt, rp_time2,f_rp2] = generate_SA_mihs(SA',m(2), false);
    else
        SA          = params.SA;
        WASt        = params.WASt;
        rp_time1    = 0;
        rp_time2    = 2;
        f_rp1       = 0;
        f_rp2       = 0;
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
r       = params.k0./m;
alpha   = (1-r).^2;
beta    = r;

%% initialization
xx      = zeros(d,maxit(1));
in_iter = zeros(1,maxit(1));
x       = x1;
xp      = x;
lami    = lam^-1;
nu      = zeros(m(1),1);
nup     = nu;
i       = 0;
while(i < 2 || (norm(x - xp)/norm(xp) >= tol(1) && i < maxit(1)))
    i      = i+1;
    bt     = A'*(b - A*x) - lam*x;
    
    %inexact subsolver
    j      = 0;
    while(j < 2 || (norm(nu - nup)/norm(nup) >= tol(2) && j < maxit(2)))
        j    = j+1;
        %sub prob gradient
        g    = SA*(bt - SA'*nu) - lam*nu;
        %sub solve problem
        dnu  = R\(R'\g);
        %momentum
        nun  = nu + alpha(2)*dnu + beta(2)*(nu - nup);
        
        %update
        nup  = nu;
        nu   = nun;
    end
    in_iter(i) = j;
    dx      = lami*(bt - (SA'*nu));
    xn      = x + alpha(1)*dx +beta(1)*(x - xp);
    xp      = x;
    x       = xn;
    xx(:,i) = x;

end
xx      = xx(:,1:i);
in_iter = in_iter(1:i);

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time1 + rp_time2;

% flop count
f_iter1 = 4*n*d + 2*m(1)*d + 9*d;
f_iter2 = 4*m(1)*d + 2*m(1)^2 + 21*m(1) + n;
flopc   = [1:i]*f_iter1 + cumsum(in_iter)*f_iter2 + f_rp1 + f_rp2 + f_qr;
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