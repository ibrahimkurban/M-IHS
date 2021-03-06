function [x,xx,time, flopc] = blendenpik_v2_iko(A,b,lam,m,x1,tol,maxit,params)
%%BLENDENPIK_IKO imlementation of the algorithm in paper:
% Avron, Haim, Petar Maymounkov, and Sivan Toledo. 
% "Blendenpik: Supercharging LAPACK's least-squares solver." 
% SIAM Journal on Scientific Computing 32.3 (2010): 1217-1236.
%
% [x,xx,time, flopc] = blendenpik_iko(A,b,lam,m,x1,tol,maxit,params)
%
% params.SA = sketch matrix otherwise it is ROS with DCT
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
%
% This is dual version,  i.e. udnerdetermined blendenpik
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%

time    = zeros(maxit+3,1);
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time, f_rp] = generate_SA_blendenpik(A',m, true);
else
    if(~isfield(params, 'SA'))
        [SAt, rp_time, f_rp] = generate_SA_blendenpik(A',m, true);
     else
        rp_time = 0;
        f_rp    = 0;
    end
end
time(1) = rp_time;

%% QR decomposition
[n,d]   = size(A);
tic;
if(lam == 0)
   [~, R] = qr(SAt,0);
   f_qr   = ceil(2*m*n^2-2/3*n^3);
else
   [~, R] = qr([SAt/sqrt(lam); eye(n)],0);
   f_qr   = ceil(2*(m+n)*n^2-2/3*n^3);
end
time(2) = toc;
%% pre condition
if(lam== 0)
    RA = R'\A;
else
    RA = R'\[A/sqrt(lam), eye(n)];
end
Rb  = R'\b;
time(2) = toc;
%% LSQR solver
[x, iter, xx, f_lsqr, time(4:end)] = lsqr_iko(RA,Rb,0,tol,maxit);

x   = x(1:d)/sqrt(lam);
xx  = xx(1:d,:)/sqrt(lam);
%% complexity (flop count refer to lightspeed malab packet)

% flop count
flopc   = f_rp + f_lsqr + f_qr + [1:iter]*(2*n^2 + 14*n);  %last term for R\

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