function [SA, time, flopc] = generate_SA_gaussian(A,m)
%%GENERATE_SA generates gauusian random matrix
%
%   [SA, time, flopc] = generate_SA_gaussian(A,m)
%
%   E[SA'SA] = I

[n,d]   = size(A);
tic;
S       = randn(m,n);
SA      = S*A/sqrt(m);
time    = toc;

%% flops if you do not have lightmaster packet do not use this output
if(nargout>2)
    flopc = 0;
    flopc = flopc + flops_randnorm(m,n);
    flopc = flopc + flops_mul(m,n,d);
end
end