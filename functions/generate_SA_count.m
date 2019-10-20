function [SA, time] = generate_SA_count(A, m)
%%GENERATE_SA generates Count sketch matrix
%
%   [SA, time, flopc] = generate_SA_mihs(A,SSIZE,wrep)
%
%   E[SA'SA] = I
tic;
n       = size(A,2);
sgn     = randi(2, n, 1) * 2 - 3;                    % one half are +1 and the rest are -1
A       = sgn.*A;
ind_r   = randperm(m);
ll      = [ind_r, randsample(m, n - m, true)'];     % sample n items from [s] with replacement
S       = sparse(ll, 1:n,1);
SA      = S*A;

time = toc;
end