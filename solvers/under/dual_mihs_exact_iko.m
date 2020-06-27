function [x,xx,time,flopc] = dual_mihs_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%DUAL_MIHS solver for underdetermined systems
%
% [x,xx,time,flopc] = dual_mihs_exact_iko(A,b,lam,m,tol,x1,maxit,params)
%
%   params.SAt sketch amtrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
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
    if(~isfield(params, 'SAt'))
        [SAt, rp_time,f_rp] = generate_SA_mihs(A',m, false);
    else
        rp_time = 0;
        f_rp    = 0;
    end
    noparam = false;
end
time(1) = rp_time;

%% QR decomposition
[n,d]   = size(A); tic;
if(lam == 0)
   [~, R] = qr(SAt,0);
   f_qr   = ceil(2*m*n^2-2/3*n^3);
else
   [~, R] = qr([SAt;sqrt(lam)*eye(n)],0);
   f_qr   = ceil(2*(m+n)*n^2-2/3*n^3);
end
time(2) = toc;
%% effective rank
if(noparam || ~isfield(params, 'k0'))
    [k0, f_tr, tr_time] = hutchinson_tri_iko(R,lam);
    params.k0 = k0;
else
    f_tr    = 0;
    tr_time = 0;
end
time(3) = tr_time;
%% momentum weights
r       = params.k0/m;

%%robust ones
% Ksup    = 1/(1-sqrt(r))^2 + 0.28;
% Kinf    = 1/(1+sqrt(r))^2 - 0.030;
% % 
% alpha   = 4/( sqrt(Ksup) + sqrt(Kinf) )^2;
% beta    = (  ( sqrt(Ksup) - sqrt(Kinf) )/( sqrt(Ksup) + sqrt(Kinf) )  )^2;
%%theoric ones
alpha   = (1-r)^2;
beta    = r;

%% data
x       = x1*0;
xx      = zeros(d, maxit);
nu      = zeros(n,1);
nup     = zeros(size(nu));
k = 0;
while(k < 2 || (norm(nu - nup)/norm(nup) >= tol && k < maxit))
    k           = k+1; tic;
    %grad
    grad        = b - A*x - lam*nu;
    
    %solve small system
    dnu         = R\(R'\grad);
    nun         = nu + alpha*dnu + beta*(nu - nup);
    
    %update
    nup         = nu;
    nu          = nun;
    
    x           = A'*nu;
    xx(:,k)     = x;
    time(k+3)   = toc;
end     
xx  = xx(:,1:k);
x   = xx(:,end);
%% complexity (flop count refer to lightspeed malab packet)

% flop count
f_iter = 2*(n^2) + 4*n*d + 22*n;
flopc   = [1:k]*f_iter + f_rp + f_qr+f_tr;  

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
