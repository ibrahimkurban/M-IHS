function [x,xx,time,flopc,in_iter,params] = mihs_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%%
%
% [x,xx,time,flopc,in_iter,params] = mihs_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%   params.SA sketch amtrix
%
%           params.submaxit sets number max iteration for subsolver (k0)
%           params.subtol   sets tolerance for subsolver (relative
%           residual)(1e-2)

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time,f_rp]   = generate_SA_mihs(A,m, false);
    noparam         = true;
else
    if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_mihs(A,m, false);
    else
        SA      = params.SA;
        rp_time = 0;
        f_rp    = 0;
    end
    noparam     = false;
end
%% effective rank
if(noparam || ~isfield(params, 'k0'))
    [k0, ~, f_tr] = hutchinson_estimator_iko(SA, SA, lam);
    params.k0 = k0;
else
    k0      = params.k0;
    f_tr    = 0;
end

%% inexact tolerance
if(noparam ||~isfield(params, 'subtol'))
    params.subtol = 1e-2;
end
if(noparam || ~isfield(params, 'submaxit'))
    params.submaxit = max(k0, 25);
end


%% data
[n,d]   = size(A);
xx      = zeros(d, maxit);
flopc   = zeros(maxit,1);
in_iter = zeros(maxit,1);


tic;
%% momentum weights
r       = k0/m;
alpha   = (1-r)^2;
beta    = r;

%% iteration
xp      = x1*0;
dx      = xp;
x       = x1;
k       = 0;
while(k < 2 || (norm(x - xp)/norm(xp) >= tol && k < maxit))
    k       = k+1;
    grad    = A'*(b-A*x) - lam*x;
    
    %solve subproblem
    [dx, in_iter(k), f_dx]      = AA_b_solver_iko(SA,grad, lam, params.subtol, params.submaxit, dx);
    xn                          = x + alpha*dx + beta*(x - xp);
    
    %update
    xp      = x;
    x       = xn;
    
    %store and count flop
    xx(:,k) = x;
    flopc(k)= f_dx(end) + 4*n*d + n + 7*d;
end
xx      = xx(:,1:k);
flopc   = flopc(1:k);
in_iter = in_iter(1:k);

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
flopc   = cumsum(flopc) + f_rp + f_tr;

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