function [U, s, V] = dsvd(A)
if(issparse(A))
    A = full(A);
end
[U, S, V] = svd(A, 'econ');
s = diag(S);
end