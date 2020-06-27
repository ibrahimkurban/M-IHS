function [x,xx,time,flopc] = acc_ihs_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%ACC_IHS_IKO implementation of the paper in
%
% Wang, Jialei, et al. 
% "Sketching meets random projection in the dual: 
% A provable recovery algorithm for big and high-dimensional data."
% Electronic Journal of Statistics 11.2 (2017): 4896-4944.
%
% [x,xx,time,flopc] = acc_ihs_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%
%
%   time(1) = rp time
%   time(2) = SA decomposition time
%   time(3) = trace estiamtion time
%   time(i) = time of ith iter
%   use cumsum
%
% there are some mistakes in the algorithm presented in the paper
% this corrected version, I have tried to use same notation
%
% in this version we use QR decompositon to solve small systems
% params.SA is sketch matrix
%
%   Ibrahim Kurban Ozaslan
%   Bilkent University 
%   MSc in EEE Dept. 
%   November 2019
%
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time,f_rp]   = generate_SA_mihs(A,m, false);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_mihs(A,m, false);
    end
end
%% data
[n,d]   = size(A);
xx      = zeros(d, maxit);
time    = zeros(maxit+3,1);
time(1) = rp_time;
%% QR decomposition
tic;
if(lam == 0)
   [~, R] = qr(SA,0);
   f_qr   = ceil(2*m*d^2-2/3*d^3);
else
   [~, R] = qr([SA;sqrt(lam)*eye(d)],0);
   f_qr   = ceil(2*(m+d)*d^2-2/3*d^3);
end
time(2) = toc;
time(3) = 0; %sd estimation
%%
tic;
w       = x1;
wp      = x1*0;
r       = - (A'*b);
u       = R\(R'\r);
p       = -u;
v       = A'*(A*p)+lam*p;
k       = 0;
while(k < 2 || (norm(w - wp)/norm(wp) >= tol && k < maxit))
    k       = k+1;
    wp      = w;
    a       = (r'*u)/(p'*v);
    w       = w + a*p;
    rp      = r;
    r       = r + a*v;
    beta    = (r'*u)/(rp'*rp);
    u       = R\(R'\r);
    p       = -u + beta*p;
    v       = A'*(A*p) + lam*p;
    
    xx(:,k) = w;
    time(k+3) = toc; tic;
end
xx      = xx(:,1:k);
x       = xx(:,end);
%% complexity (flop count refer to lightspeed malab packet)
% flop count
f_init = 2*d^2 + 6*n*d + 17*d;
f_iter = 2*(d^2) + 4*n*d + 16 + 31*d;
flopc  = [1:k]*f_iter + f_init + f_rp + f_qr;   

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