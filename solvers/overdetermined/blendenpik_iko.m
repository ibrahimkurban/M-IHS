function [x,xx,time] = blendenpik_iko(A,b,lam,m,x1,tol,maxit,params)
%%BLENDENPIK_IKO imlementation of the algorithm in paper:
% Avron, Haim, Petar Maymounkov, and Sivan Toledo. 
% "Blendenpik: Supercharging LAPACK's least-squares solver." 
% SIAM Journal on Scientific Computing 32.3 (2010): 1217-1236.
%
% [x,xx,time] = blendenpik_iko(A,b,lam,m,x1,tol,maxit,params)
%
% params.SA = sketch matrix otherwise it is ROS with DCT
%

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time] = generate_SA_blendenpik(A,m, true);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time] = generate_SA_blendenpik(A,m, true);
    end
end
tic

%% QR decomposition
if(lam == 0)
   [~, R] = qr(SA,0);
else
   [~, R] = qr([SA;sqrt(lam)*eye(d)],0);
end

%% LSQR Solver
[x, ~, xx] = lsqr_pre_iko([A; sqrt(lam)],[b; zeros(d)],R,x1,tol,maxit);
time = toc+rp_time;


end






%% IF SA is not provided
function [SA, time, flopc] = generate_SA_blendenpik(A,SSIZE,wrep)
%%GENERATE_SA generates ROS sketch matrix
%
%   [SA, time, flopc] = generate_SA_blendenpik(A,SSIZE,wrep)
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