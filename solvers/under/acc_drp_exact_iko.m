function [x,xx,time,flopc] = acc_drp_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%%ACC_IHS_IKO implementation of the paper in
%
% Wang, Jialei, et al. 
% "Sketching meets random projection in the dual: 
% A provable recovery algorithm for big and high-dimensional data."
% Electronic Journal of Statistics 11.2 (2017): 4896-4944.
%
% [xx, time, flopc] = acc_ihs_exact_iko(A,b,lam,m,x1,tol,maxit,params)
%
% there are some mistakes in the algorithm presented in the paper
% this corrected version, I have tried to use same notation
%
% in this version we use QR decompositon to oslve small systems
% params.SAt = sketch matrix

%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SAt, rp_time,f_rp]   = generate_SA_mihs(A',m, false);
else
    if(~isfield(params, 'SA'))
        [SAt, rp_time,f_rp] = generate_SA_mihs(A',m, false);
    end
end
tic;

%% QR decomposition
[n,d]   = size(A);
if(lam == 0)
   [~, R] = qr(SAt,0);
   f_qr   = ceil(2*m*n^2-2/3*n^3);
else
   [~, R] = qr([SAt;sqrt(lam)*speye(n)],0);
   f_qr   = ceil(2*(m+n)*n^2-2/3*n^3);
end

%%
xx      = zeros(d,maxit);
alpha   = zeros(n,1);
r       = (-lam)*b;
z       = SAt*(SR\(SR.'\r));
u       = r - SAt'*z;
p       = -u;
v       = A*(A'*p) + lam*p;

k = 0;
while(k < 2 || (norm(dffxp)/norm(xp) >= tol && k < maxit))
    k       = k+1;
    a       = (r'*u)/(p'*v);
    dffxp   = a*p; %for tolerance
    xp      = alpha;
    alpha   = alpha + a*p;
    rp      = r;
    r       = r + a*v;
    z       = SAt*(SR\(SR'\r));
    beta    = (r'*u)/(rp'*rp);
    u       = r - SAt'*z;
    p       = -u+beta*p;
    v       = A*(A'*p) + lam*p;
    
    xx(:,k)  = (A'*alpha)/lam;
end
xx      = xx(:,1:k);
x       = xx(:,end);


%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
f_init  = 2*n^2+ 4*n*d + 4*m*n + 12*n + 9*d + 16 ;
f_iter  = 2*(n^2) + 4*n*d + 4*m*n + 5*n;
flopc   = [1:k]*f_iter + f_rp + f_qr + f_init;    

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

%% flop count refer to lightspeed malab packet
if(nargout > 2)
    f_DA    = 18*n + n*d;
    f_HDA   = ceil(nt*d*log2(7*m));
    f_SA    = 18*nt + m*d; 
    
    flopc   = f_DA + f_HDA + f_SA;
end
end