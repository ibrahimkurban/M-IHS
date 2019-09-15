function [x,xx,time, flopc] = blendenpik_iko(A,b,lam,m,x1,tol,maxit,params)
%%BLENDENPIK_IKO imlementation of the algorithm in paper:
% Avron, Haim, Petar Maymounkov, and Sivan Toledo. 
% "Blendenpik: Supercharging LAPACK's least-squares solver." 
% SIAM Journal on Scientific Computing 32.3 (2010): 1217-1236.
%
% [x,xx,time, flopc] = blendenpik_iko(A,b,lam,m,x1,tol,maxit,params)
%
% params.SA = sketch matrix otherwise it is ROS with DCT
%

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time, f_rp] = generate_SA_blendenpik(A,m, true);
else
        if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_mihs(A,m, false);
    else
        SA      = params.SA;
        rp_time = 0;
        f_rp    = 0;
    end
end
tic

%% QR decomposition
[~,d]   = size(A);
if(lam == 0)
   [~, R] = qr(SA,0);
   f_qr   = ceil(2*m*d^2-2/3*d^3);
else
   [~, R] = qr([SA;sqrt(lam)*eye(d)],0);
   f_qr   = ceil(2*(m+d)*d^2-2/3*d^3);
end

%% LSQR Solver
if(lam == 0)
    [x, ~, xx, f_lsqr] = lsqr_pre_iko(A,b,0,R,tol,maxit);
else
    [x, ~, xx, f_lsqr] = lsqr_pre_iko([A; sqrt(lam)*speye(d)],[b; zeros(d,1)],0,R,tol,maxit);
end
%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
flopc   = f_rp + f_lsqr + f_qr; 

end






%% IF SA is not provided
function [SA, time, flopc] = generate_SA_blendenpik(A,m,wrep)
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