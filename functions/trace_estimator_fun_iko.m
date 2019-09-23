function tr_est_fun = trace_estimator_fun_iko(lambdas, RR)

%% check whether lambdasi logarithmic or not
if(diff(diff(lambdas(end-2:end))) < 1e-10)
    islog = true;
else
    islog = false;
end

L       = length(lambdas);
tr      = zeros(L,1);
I_L     = speye(size(RR,1));

switch islog
    %% logarithmic
    case true
        for i = 1:L
            lam     = lambdas(i);
            RR_I    = RR + 10^lam*I_L;
            tr(i)   = solver_trace_bidiag_inv(RR_I);
        end
        pp          = interp1((lambdas), tr, 'pchip', 'pp');
        tr_est_fun  = @(lam)ppval(pp, lam);
        
        
    otherwise
        %% not
        for i = 1:L
            lam     = lambdas(i);
            RR_I    = RR + lam*I_L;
            tr(i)   = solver_trace_bidiag_inv(RR_I);
        end
        pp          = interp1(log10(lambdas), tr, 'pchip', 'pp');
        tr_est_fun  = @(lam)ppval(pp, log10(lam));
end
end


%%
function trace = solver_trace_bidiag_inv(RR)
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