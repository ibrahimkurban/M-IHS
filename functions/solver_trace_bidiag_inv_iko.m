function trace = solver_trace_bidiag_inv_iko(RR)
%%TRACE_LANCZOS finds trace((R*R'+lam*eye)^-1) in 7k operations
%
% trace = trace_lanczoz(RR)
%
%   RR = R*R'*lam*eye(k), can be sparse matrix
%
tt          = diag(RR, 1).^2;
rr          = diag(RR);
n           = length(rr);
trace       = 1/rr(n);
trace_i     = trace;
rho2        = rr(n);
for i = n-1:-1:1
    theta2  =  tt(i)/rho2;
    rho2    = rr(i) - theta2;
    trace_i = (1+theta2*trace_i)/rho2;
    trace   = trace + trace_i;
end
end