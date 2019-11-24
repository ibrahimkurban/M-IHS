function [SSA, time, flopc] = generate_SA_op(A,SSIZE,N, type, sparseS, wrep)
%%GENERATE_SA generates random projected matrix SA for N times by using only one
%%rademacher randomization, uses ROS sketch
%
%   [SSA, time, flopc] = generate_SA_op(A,SSIZE,N, type, sparseS, wrep)
%
%   SA = SSA(:,:,i), i = 1,2,3,...N
%   E[SA'SA] = I
%
%   type    : 'dct' or 'hadamard'
%   sparseS : boolean makes output sparse or not
%   wrep    : boolean (false, default)without replacement or (true)with replacement
%
if(~exist('wrep', 'var'))
    wrep = false;
end
[n,d]       = size(A);
SSA         = zeros(SSIZE,d,N);
nt          = 2^(ceil(log2(n)));
idx         = randsample(nt, SSIZE*N, wrep);                 % sampling pattern
idx_c       = 1:nt;
idx_c(idx)  = [];
idx_p       = randperm(numel(idx), numel(idx));
P           = opPermutation(idx_p);
switch type
    case 'hadamard'
        H = P*opExcise(opHadamard(nt)/sqrt(nt), idx_c, 'rows');
    case 'dct'
        H = P*opExcise(opDCT(nt), idx_c, 'rows');
end

radem   = (randi(2, n, 1) * 2 - 3);                     % rademacher
DA      = A .* radem;                                   % one half+1 and rest -1                                    % DCT transform

if(n == nt)
    SSAt    = (H*DA)*(sqrt(nt)/sqrt(SSIZE));              % subsampling
else
    SSAt    = (H*[DA; zeros(nt-n, d)])*(sqrt(nt)/sqrt(SSIZE));              % subsampling
end
if(sparseS)
    SSAt = sparse(SSAt);
end
time = 0;
if(N > 1)
    for k = 1:N
        SSA(:,:,k) = SSAt(1+(k-1)*SSIZE:k*SSIZE, :);
    end
else
    SSA = SSAt;
end

%% flops
if(exist('flops_randnorm', 'file') > 0 && nargout>2)
    flopc = 0;
    flopc = flopc + flops_randnorm(n);
    flopc = flopc + flops_fft(DA)/4;
    flopc = flopc + flops_sum(SSAt)/N;
end
end