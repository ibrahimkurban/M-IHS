function [par, dev_est] = LS_sgcv_lower_gkl_iko(R,b, theta1)
%%LS_SGCV calcualtes GCV for given U,sig b
%
%
%  [par, dev_est] = LS_sgcv_lower_gkl(R,b)
%
%

%% constuct bidiagonal system Ry =f
L       = size(R,1);
RR      = R*R';
n       = length(b);

bnrm2   = norm(b)^2;
f       = inv_ubidiag_e1(theta1, diag(R,1), diag(R));
delta2  = bnrm2 - norm(f)^2;
if((delta2) < 1e-14)
    delta2 = 0;
end

%% minimization of curve
options                 = optimset('TolX', 1e-4, 'Display','off');
[par_log,gcv_val]       = fminbnd(@(lam)gcv_gkl(10^lam, RR, f, delta2,n-L), -10, 3, options);
par                     = 10^par_log;
dev_est                 = sqrt(gcv_val)*sqrt(n-L+ par*solver_trace_bidiag_inv_iko(RR+par*speye(L)));

end


function r = inv_ubidiag_e1(t1, tt,rr)

hh = [t1; tt]./rr;
hh(2:end) = -hh(2:end);
r = cumprod(hh);
end

function gcv = gcv_gkl(lam, RR, f, delta2,nd)

RR_I    = RR + lam*speye(size(RR));
tr      = solver_trace_bidiag_inv_iko(RR_I);
gcv     = (norm(lam*(RR_I\f))^2 + delta2)/(nd + lam*tr)^2;

end