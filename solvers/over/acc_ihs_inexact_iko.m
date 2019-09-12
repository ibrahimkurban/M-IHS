function [xx, time, flopc] = acc_ihs_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%%ACC_IHS_IKO implementation of the paper in
%
% Wang, Jialei, et al. 
% "Sketching meets random projection in the dual: 
% A provable recovery algorithm for big and high-dimensional data."
% Electronic Journal of Statistics 11.2 (2017): 4896-4944.
%
% [xx, time, flopc] = acc_ihs_inexact_iko(A,b,lam,m,x1,tol,maxit,params)
%
% there are some mistakes in the algorithm presented in the paper
% this corrected version, I have tried to use same notation
%
% in this version we use AAb inexact solver proposed by IKO which is not suggested in the paper
%
%           params.submaxit sets number max iteration for subsolver (25)
%           params.subtol   sets tolerance for subsolver (relative
%           residual)(1e-2)

%% inexact tolerance 
    if(~isfield(params, 'subtol'))
        params.subtol = 1e-2;
    end
    if(~isfield(params, 'submaxit'))
        params.maxit = 25;
    end
    
%% generate sketch matrix or not
if(~exist('params', 'var'))
    [SA, rp_time,f_rp]   = generate_SA_mihs(A,m, false);
else
    if(~isfield(params, 'SA'))
        [SA, rp_time,f_rp] = generate_SA_mihs(A,m, false);
    end
end
tic;

%% data
[n,d]   = size(A);
xx      = zeros(d, maxit);
f_iterr = zeros(1,maxit);
%%
w       = x1;
wp      = x1*0;
r       = - (A'*b);
[u,~,f_idx]= AA_b_solver(SA,r, lam, cgtol, cgmaxit, dx);
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
    [u,~,f_dx]= AA_b_solver(SA,r, lam, cgtol, cgmaxit, dx);
    p       = -u + beta*p;
    v       = A'*(A*p) + lam*p;
    
    xx(:,k) = w;
    f_iterr(k) = f_dx(end);
end
xx      = xx(:,1:k);

%% complexity (flop count refer to lightspeed malab packet)
%timing
time    = toc+rp_time;

% flop count
f_init = 2*d^2 + 6*n*d + 17*d;
f_iter = 2*(d^2) + 4*n*d + 16 + 31*d;
flopc  = [1:k]*f_iter + f_init + f_rp + f_qr + f_iterr(1:k) + f_idx(end);   


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
