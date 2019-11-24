function [SSA, time, flopc] = generate_SA(A,SSIZE,N, wrep)
%%GENERATE_SA generates random projected matrix SA for N times by using only one
%%rademacher randomization, uses ROS sketch
%
%   [SSA, time] = random_matrix_v2(A,m,N)
%
%   SA = SSA(:,:,i), i = 1,2,3,...N
%   E[SA'SA] = I

if(issparse(A))
    A = full(A);
end
if(~exist('wrep', 'var'))
    wrep = false;
end
[n,d]   = size(A);
SSA     = zeros(SSIZE,d,N); 
nt      = n;ceil(n/1000)*1000; %BLENDENPIK sugg. for FFTW
tic;
radem   = (randi(2, n, 1) * 2 - 3);                     % rademacher
DA      = A .* radem;                                   % one half+1 and rest -1
HDA     = dct(DA,nt);                                      % DCT transform
idx     = randsample(nt, SSIZE*N, wrep);                 % sampling pattern
SSAt    = HDA(idx, :)*(sqrt(nt)/sqrt(SSIZE));              % subsampling
time    = toc;
for k = 1:N
    SSA(:,:,k) = SSAt(1+(k-1)*SSIZE:k*SSIZE, :);
end

%% flops
if(nargout>2)
    flopc = 0;
    flopc = flopc + flops_randnorm(n);
    flopc = flopc + flops_fft(DA)/4;
    flopc = flopc + flops_sum(SSAt)/N;
end
end