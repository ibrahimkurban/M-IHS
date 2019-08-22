function [x, iter, flopc] = cg_templates(A,b,x,M,tol,maxit,n2)

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
%
% cg.m solves the symmetric positive definite linear system Ax=b
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix, [] if not used
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%
% if A'A =b is tried to solve then n2 is the  size(A,1) otherwise ignore it
%
flag    = 0;                                 % initialization
pre     = ~isempty(M);
iter    = 0;
bnrm2   = norm( b );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end


r = b - A*x;
error = norm( r ) / bnrm2;
if ( error < tol ) return, end

for iter = 1:maxit                       % begin iteration
    if(pre)
        z  = M \ r;
    else
        z = r;
    end
    rho = (r'*z);
    
    if ( iter > 1 )                      % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
    else
        p = z;
    end
    
    q = A*p;
    alpha = rho / (p'*q );
    x = x + alpha * p;                    % update approximation vector
    
    r = r - alpha*q;                      % compute residual
    error = norm( r ) / bnrm2;            % check convergence
    if ( error <= tol ), break, end
    
    rho_1 = rho;
    
end

if ( error > tol ) flag = 1; end         % no convergence

%% flop count refer to lightspeed malab packet
if(nargout > 2)
    n = size(b);
    if(pre)
        f_iter = 3*n^2 + 17*n + 16;
    else
        f_iter = 2*n^2 + 10*n + 16;
        if(exist('n2', 'var'))
            f_iter = 4*n*n2 + + 12*n + 16;
        end
    end
    f_init    = 2*n^2 + 8*n + 14;
    if(exist('n2', 'var'))
        f_init = 4*n*n2 + + 10*n + 16;
    end
    flopc   = [1:i]*f_iter + f_init;
end

end
% END cg.m