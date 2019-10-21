function [x,xx,time,flopc] = chebyshev_pre_iko(A,b,lam,m,x1,tol,maxit,params)
%%CHEBYSHEV_IKO chebyshev semiconvergence method with preconditioner
% (A'S'SA + lam I)^{-1}
%
% [x,xx,time,flopc] = chebyshev_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SA sketch amtrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%

[n,d]   = size(A);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time,f_rp]   = generate_SA_mihs(A,m, false);
    params.k0            = d;
else
    if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_mihs(A,m, false);
    else
        SA      = params.SA;
        rp_time = 0;
        f_rp    = 0;
    end
    if(~isfield(params, 'k0'))
        params.k0 = d;
    end
end
tic;
%% QR decomposition
if(lam == 0)
    [~, R] = qr(SA,0);
    f_qr   = ceil(2*m*d^2-2/3*d^3);
else
    [~, R] = qr([SA;sqrt(lam)*speye(d)],0);
    f_qr   = ceil(2*(m+d)*d^2-2/3*d^3);
end

%% momentum weights
r       = params.k0/m;

Ksup    = 1/(1-sqrt(r))^2;% + 0.28;
Kinf    = 1/(1+sqrt(r))^2;% - 0.030;

%%cosntants
dc = (Ksup + Kinf)/2;
cc = (Ksup - Kinf)/2;

%% iteration
xx      = zeros(d, maxit);
x       = x1;
xp      = x1*0;
i       = 0;
while(i < 2 || (norm(xp - x)/norm(xp) >= tol && i < maxit))
    i = i+1;
    %gradient
    grad    = A'*(b-A*x) - lam*x;
    
    %preconditioning
    z       = (R\(R'\grad));
    
    %parameters
    if (i ==1)
        p       = z;
        alpha   = 2/dc;
    else
        beta    = (cc*alpha/2)^2;
        alpha   = 1/(dc- beta);
        p       = z + beta*p;
    end
    %update
    xp= x;
    x = x + alpha*p;
    
    xx(:,i) = x;
end

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
f_iter  = 2*(d^2) + 4*n*d + n + 18*d;
flopc   = [1:maxit]*f_iter + f_rp + f_qr;

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

