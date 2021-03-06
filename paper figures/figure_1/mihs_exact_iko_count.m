function [x,xx,time,flopc] = mihs_exact_iko_count(A,b,lam,m,x1,tol,maxit,params)
%%MIHS_EXACT_IKO solver for over-determined systems
%
% [x, xx, time, flopc] = mihs_exact_iko(A,b,lam,m,tol,maxit,params)
%
%   params.SA sketch amtrix
%   params.k0 effective rank (reqirement if not known run eff_rank_solver)
%

[n,d]   = size(A);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time]   = generate_SA_count(A,m);
    params.k0            = d;
    f_rp = nnz(A);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time] = generate_SA_count(A,m);
        f_rp = nnz(A);
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
xx      = zeros(d, maxit);

%% iteration
xp      = x1*0;
x       = x1;
k       = 0;
while(k < 2 || (norm(x - xp)/norm(xp) >= tol && k < maxit))
    k       = k+1;
    grad    = A'*(b-A*x) - lam*x;
    xn      = x + alpha*(R\(R'\grad)) + beta*(x - xp);
        
    %update
    xp      = x;
    x       = xn;
    xx(:,k) = x; 
end
xx      = xx(:,1:k);


%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
f_iter = 2*(d^2) + 4*n*d + n + 21*d;
flopc   = [1:k]*f_iter + f_rp + f_qr;  

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