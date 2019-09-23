function [xxOR, parOR, errOR, err_rid] = oracle_methods(U, sig, V,b,x0,x1)
% input data is structure
%     data.A = A;
%     data.b = b;
%     data.x0 = x0;
%     data.x1 = x1;
%     data.dev = NOISE_DEV;
%     data.MC = numel(unique(NOISE_DEV));
fprintf('Oracle Methods run....'); tic;
N       = size(b,2);
d       = length(x1);
xxOR    = {2,1};
parOR   = zeros(2,N);
err_rid = zeros(1,N);
for i=1:2, xxOR{i} = zeros(d,N); end
% for i=1:2, parOR{i} = zeros(N,1); end

%% EFFECTIVE SOLUTION
% [U, sig, V]   = dsvd(A);

for i = 1:N
    [xxOR{1}(:,i), parOR(1, i)]   = LS_truncated_iko(U,sig,V, b(:,i), x1, x0 );
end
errOR{1} = @(xx, dev)(sqrt(sum((xxOR{1}(:,dev) - xx).^2, 1))/norm(xxOR{1}(:,dev)));
for i =1:N
    [xxOR{2}(:,i), parOR(2, i)]   = LS_ridge_iko(U,sig,V, b(:,i), x1, @(x)errOR{1}(x, i));
    err_rid(i) = errOR{1}(xxOR{2}(:,i), i);
end
errOR{2} = @(xx, dev)(sqrt(sum((xxOR{2}(:,dev) - xx).^2, 1))/norm(xxOR{2}(:,dev)));
fprintf('done %2.2f sec elapsed\n', toc);
end