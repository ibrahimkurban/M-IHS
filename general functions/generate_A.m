function [A, SINGULAR_VAL,x0] = generate_A(n, d, varargin)
%%DATA_GENERATE generates data with given specs
%
%    A = mxn
%
%    [A, SINGULAR_VAL] = generate_A(m, n, varargin)
%
%   VARARGIN
%
%       'singular values'(randn singular val)      :   double(n,1)
%       'correlated'(false)                        :   boolean
%       'matrix correlation'(.6)                   :   double, 0< <1
%       'matrix deviation'(3)                      :   double deviation ofA
%                                                       sugg.(3*n)
%       'hansen singular values'(0)DEF2            : double scaler in 1,..10
%       'hansen max'(1e6)                          : double value of sigma_1
%       'kappa'(6)                                 : double in log10 scale
%       'hasnen input'(false)                       : return x0 if exist
%
%
%       DEF2
%       1)deriv2,       w/   input
%       2)gravity       w/   input
%       3)heat          w/   input
%       4)show          w/   input
%       5)philips       w/   input
%       6)i_lap.        w/   input
%       7)foxgood       w/  input
%       8)baart         w   input
%       9)blurred image w/o input
%       10)parallax     w/o  input
fprintf('Data Generation runs...')
%% DEFAULT values
KAPPA           = 6;
CORRELATED      = false;
CORRELATION     = .6;
A_DEV           = 3;
HANSEN          = 0;
HANSEN_MAX      = 1e6;
HANSEN_IN       = false;
%% user input
for i = 1:2:length(varargin)
    str = varargin{i};
    switch str
        case 'singular values'
            SINGULAR_VAL = varargin{i+1};
        case 'hansen singular values'
            HANSEN = varargin{i+1};
        case 'hansen input'
            HANSEN_IN = varargin{i+1};
        case 'hansen max'
            HANSEN_MAX = varargin{i+1};
        case 'kappa'
            KAPPA = varargin{i+1};
        case 'correlated'
            CORRELATED = varargin{i+1};
        case 'matrix correlation'
            CORRELATION = varargin{i+1};
        case 'matrix deviation'
            A_DEV = varargin{i+1};
        otherwise
            error(['generate data: un-identified input!!!:  ' varargin{i}])
    end
end

%% HANSEN DATA
if(HANSEN > 0)
    HAN_DATA = load('hansen_data.mat', 'HANSEN_CUT');
    HAN_DATA = HAN_DATA.HANSEN_CUT;
    sing_han = HAN_DATA{HANSEN,1};
    
    %%%SAMPLE HANSEN DATA
    n_han = size(sing_han,1);
    ind_han = linspace(1, n_han, min(n,d));
    SINGULAR_VAL = interp1(1:n_han, sing_han, ind_han, 'linear');
    SINGULAR_VAL = SINGULAR_VAL*HANSEN_MAX/max(SINGULAR_VAL);
    SIG_LOG = log10(SINGULAR_VAL);
    SIG_REC = KAPPA*SIG_LOG/abs(SIG_LOG(end));
    SINGULAR_VAL = 10.^SIG_REC;
end
switch CORRELATED
    case 0 % un correlated
        A = A_DEV*randn(n,d);
    case 1 % correlated
        n1       = ones(n,d);
        R        = ihs_corr(d, CORRELATION, A_DEV);
        T        = chol(R);
        A        = randn(n,d)*T+n1;
end
[U,~,V]  = svd(A, 'econ');
A        = U * diag(SINGULAR_VAL) * V';
%% INPUT
if(HANSEN_IN)
    han_x = size(HAN_DATA{HANSEN,2},1);
    if(han_x > 1)
        x_han = HAN_DATA{HANSEN,2};
        x0 = interp1(1:n_han, x_han, linspace(1, n_han, d), 'linear');
        x0 = reshape(x0, d, 1);
    else
        x0 = rand(d,1);
    end
end


fprintf('done\n')
end

function R = ihs_corr(n, a,dev)
%%IHS_CORR generate correlation matrix defined in IHS paper
%
%   R_jk = dev*a^|j-k|
%
%   R = ihs_corr(n, a)
%
%   R = nxn
%
P = zeros(n);
ind = 0:n-1;
for i = 1:n
    P(i:end,i) = ind(1:end-i+1);
end
P = P + P.';
R = dev*(a*ones(n)).^P;
end