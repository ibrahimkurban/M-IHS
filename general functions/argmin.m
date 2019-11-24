function [max_arg] = argmin(A, dim)
if(nargin > 1)
[~, max_arg] = min(A, [], dim);
else
    [~, max_arg] = min(A);
end