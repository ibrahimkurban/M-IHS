function [x,xx,time,flopc] = acc_ipds_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%ACC_IPDS_EXACT_IKO solver for square-under cases
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

%% initizlize
xx      = zeros(d,maxit(1));
in_iter = zeros(1, maxit(1));
aD      = zeros(n,1);
rD      = (-lam)*b;
zt      = zeros(m(1),1);
ztp     = zt;
rP      = SAt*rD;
uP      = R\(R'\rP);
pP      = -uP;
vP      = SAt*(SAt'*pP) + lam*pP;
j       = 0;
while(j < 2 || (norm(zt - ztp)/norm(ztp) >= tol(2) && j < maxit(2)))
    j       = j+1;
    ztp     = zt;
    alphaP  = sum(rP.*uP)/sum(pP.*vP);
    zt      = zt + alphaP*pP;
    norm2rP = sum(rP.^2);
    rP      = rP + alphaP*vP;
    betaP   = sum(rP.*uP)/norm2rP;
    uP      = R\(R'\rP);
    pP      = -uP + betaP*pP;
    vP      = SAt*(SAt'*pP) + lam*pP;
end
init_iter = j;
uD  = (rD - SAt'*zt);
pD  = -uD;
vD  = A*(A'*pD) + lam*pD;
k   = 0;
while(k < 2 || (norm(alphaD*pD)/norm(aDp) >= tol(1) && k < maxit(1)))
    k       = k+1;
    alphaD  = sum(rD.*uD)/sum(pD.*vD);
    aDp     = aD;
    aD      = aD + alphaD*pD;
    x       = A'*aD/lam;
    norm2rD = sum(rD.^2);
    rD      = rD + alphaD*vD;
    %%
    zt      = zeros(m(1),1);
    ztp     = zt;
    rP      = -SAt*rD;
    uP      = R\(R'\rP);
    pP      = -uP;
    vP      = SAt*(SAt'*pP) + lam*pP;
    j       = 0;
    while(j < 2 || (norm(zt - ztp)/norm(ztp) >= tol(2) && j < maxit(2)))
        j       = j+1;
        ztp     = zt;
        alphaP  = sum(rP.*uP)/sum(pP.*vP);
        zt      = zt + alphaP*pP;
        norm2rP = sum(rP.^2);
        rP      = rP + alphaP*vP;
        betaP   = sum(rP.*uP)/norm2rP;
        uP      = R\(R'\rP);
        pP      = -uP+betaP*pP;
        vP      = SAt*(SAt'*pP) + lam*pP;
    end
    in_iter(k) = j;
    betaD   = sum(rD.*uD)/norm2rD;
    uD      = rD - SAt'*zt;
    pD      = -uD + betaD*pD;
    vD      = A*(A'*pD) + lam*pD;
    
    xx(:,k) = x;
end
xx      = xx(:,1:k);
x       = xx(:,end);
in_iter = in_iter(1:k);
%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time1 + rp_time2;

% flop count
f_init1 = 8*m(1)*n + 4*n*d + 2*m(1)^2 + 17*m(1) + 6*n;
f_init2 = init_iter*(30*m(1) + 4*m(1)*n + 2*m(1)^2 + 8);
f_iter1 = 34*m(1) + 2*m(1)^2 + 8*m(1)*n + 4*n*d+2*m(1)^2+16;
f_iter2 = (28 + 2*m(1) + 4*n)*m(1)+8;
flopc   = [1:k]*f_iter1 + cumsum(in_iter)*f_iter2 + f_rp1 + f_rp2 + f_qr + f_init1 + f_init2;
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


