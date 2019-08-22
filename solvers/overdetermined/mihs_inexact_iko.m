function [xx, time] = mihs_inexact_iko(A,b,lam,m,tol,maxit,params)
%%
%
% [xx, time] = mihs_inexact(A,b,x1,lam,m,tol,maxit,k0, cgtol, cgmaxit)
%
%   params.SA sketch amtrix
%   params.exact set inexact or exact solver
%
%   if params.exact = false then 
%           params.submaxit sets number max iteration for subsolver
%           params.subtol   sets tolerance for subsolver (relative residual)

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time]   = generate_SA_mihs(A,m, false);
    noparam         = true;
else
    if(~isfield(params, 'SA'))
        [SA, rp_time] = generate_SA_mihs(A,m, false);
    end
    noparam     = false;
end

%% inexact or not
if(noparam || ~isfield(params, 'exact'))
    params.exact = true;
end

%% if inexact tolerance 
if(~params.exact)
    if(~isfield(params, 'subtol'))
        params.subtol = 1e-2;
    end
    if(~isfield(params, 'submaxit'))
        params.maxit = 25;
    end
end
   
tic
%% data
[n,d]   = size(A);
xx      = zeros(d, maxit);
time    = zeros(maxit,1);
tic;

%% effective rank
if(noparam || ~isfield(params, 'k0'))
    k0 = eff_rank_solver(SA, SA, lam);
end

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
    
    switch params.exact
        case false
            %inexact version
            dx      = AA_b_solver(SA,grad, lam, cgtol, cgmaxit, dx);
        case true
            %exact version
            dx      = (SA'*SA + lam*speye(d))\grad;
    end
    xn      = x + alpha*dx + beta*(x - xp);
    
    %update
    xp      = x;
    x       = xn;
    xx(:,k) = x;
    time(k) = toc;
    
end
xx      = xx(:,1:k);
time    = time(1:k);
fprintf('%3.0e elapsed\n', toc);

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

%% flops
if(nargout>2)
    flopc = 0;
    flopc = flopc + flops_randnorm(n);
    flopc = flopc + flops_fft(DA)/4;
    flopc = flopc + flops_sum(SA)/N;
end
end