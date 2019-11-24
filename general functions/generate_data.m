function [data, b, x0, x1, DEV, U, sigA, V ] = generate_data(n, d, varargin)
%%DATA_GENERATE generates data with given specs
%
%    A = mxn
%
%    [NAME, SAVE_DIR, A, y, x0, x00, A_corr, y_corr, err_x] = generate_data(m, n, varargin)
%
%   VARARGIN
%
%       'singular values'(randn singular val)      :   double(n,1)
%       'correlated'(false)                        :   boolean
%       'matrix correlation'(.6)                   :   double, 0< <1
%       'matrix deviation'(3)                      :   double deviation ofA
%                                                       sugg.(3*n)
%       'hansen singular values'(0)DEF2            : double scaler in 1,..10
%       'hansen in'(false)                         : input is created if hansen input exists
%       'hansen max'(1e6)                          : double value of sigma_1
%       'kappa'(6)                                 : double in log10 scale
%       'input dev'(1)                                : double
%       'prior dev'(0)
%       'prior mean'(0)
%       'noise level'(0)
%       'structure'(false)
%       'plot'(false)                               : plot or not
%
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
fprintf('Data Generation runs...'); tic
%% DEFAULT values
KAPPA           = 6;
CORRELATED      = false;
CORRELATION     = .6;
A_DEV           = 3;
HANSEN_IN       = false;
HANSEN          = 0;
HANSEN_MAX      = 1e6;
INPUT_DEV       = 1;
PRIOR_MEAN      = 0;
PRIOR_DEV       = 0;
NOISE_LEVEL     = 0;
STRUCT          = false;
PLOT            = false;
%% user input
for i = 1:2:length(varargin)
    str = varargin{i};
    switch str
        case 'singular values'
            sigA = varargin{i+1};
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
        case 'input dev'
            INPUT_DEV = varargin{i+1};
        case 'prior mean'
            PRIOR_MEAN = varargin{i+1};
        case 'prior dev'
            PRIOR_DEV = varargin{i+1};
        case 'noise level'
            NOISE_LEVEL = varargin{i+1};
        case 'structure'
            STRUCT = varargin{i+1};
        case 'plot'
            PLOT = varargin{i+1};
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
    sigA = interp1(1:n_han, sing_han, ind_han, 'linear');
    sigA = sigA/max(sigA);
    SIG_LOG = log10(sigA);
    SIG_REC = KAPPA*SIG_LOG/abs(SIG_LOG(end));
    sigA = (10.^SIG_REC)*HANSEN_MAX;
end
switch CORRELATED
    case 0 % un correlated
        A = A_DEV*randn(n,d);
    case 1 % correlated
        n1       = ones(n,d);
        R        = ihs_corr(min(n,d), CORRELATION, A_DEV);
        T        = chol(R);
        if(n>=d)
            A        = randn(n,d)*T+n1;
        else
            A        = T*randn(n,d)+n1;
        end
end
[U,~,V]  = svd(A, 'econ');
A        = U * (sigA'.* V');

%% INPUT
x1 = PRIOR_DEV*randn(d,1) + PRIOR_MEAN;
x1 = V*(V'*x1);
if(HANSEN_IN)
    han_x = size(HAN_DATA{HANSEN,2},1);
    if(han_x > 1)
        x_han = HAN_DATA{HANSEN,2};
        x0 = interp1(1:n_han, x_han, linspace(1, n_han, d), 'linear');
        x0 = reshape(x0, d, 1);
    else
        x0 = x1 + INPUT_DEV*randn(d,1);
    end
else
    x0 = x1 + INPUT_DEV*randn(d,1);
%      x0 = x1 + INPUT_DEV*sum(V,2);
end
x0 = V*(V'*x0);
%% MEASUREMENT
b0      = A*x0;
b0_nrm  = norm(b0);
%%NOISE
w = zeros(n, numel(NOISE_LEVEL));
DEV = zeros(numel(NOISE_LEVEL),1);
for i = 1:numel(NOISE_LEVEL)
    r           = randn(n,1);
    DEV(i)      = (NOISE_LEVEL(i)*b0_nrm)/norm(r);
    w(:,i)      = DEV(i)*r;
end

%% NOISY MEASUREMENT
b = b0 + w;
sigA = sigA(:);

%% OUTPUT
if(STRUCT)
    data.size       = [n, d];
    data.A          = A;
    data.b          = b;
    data.x0         = x0;
    data.x1         = x1;
    data.dev        = DEV;
    data.MC         = numel(unique(NOISE_LEVEL));
    data.metric     = @(xx)(sqrt(sum((x0 - xx).^2, 1))/norm(x0));
    data.nlevel     = NOISE_LEVEL;
    data.matdev     = A_DEV;
    data.matcor     = CORRELATION;
else
    data = A;
end

%% PLOT
if(PLOT)
    D = length(DEV);
    if D ==1
        D = [];
    end
    if(STRUCT)
        data.fig = figure('Name', 'Theoric Plot');
    else
        fig = figure('Name', 'Theoric Plot');
    end
    subplot(1,3,1);
    %coherence
    imagesc(normA(A), [0 1]); colorbar;
    title('coherence of A: $<\hat{a_i}, \hat{a_j}>$', 'Interpreter', 'latex');
    xlabel('column index i'); ylabel('column index j')
    axis square
    set(gca,'fontsize',16)
    
    %singular values
    subplot(1,3,2);
    semilogy(1:min(d,n), sigA, 'linewidth', 3); hold on;
    semilogy(1:min(d,n), abs(U'*b(:,[D, 1])), 'linewidth', .5); hold on;
    semilogy(1:min(d,n), DEV([D, 1]).*ones(min(d,n),1), 'linewidth', .5); hold on;
    title('Energy of Measurement'); xlabel('index i');
    if(numel(DEV) == 1)
        legend('\sigma_i', '|u_i^Ty|', '\sigma_\omega');
    else
        legend('\sigma_i', 'L:|u_i^Ty|', 'S:|u_i^Ty|', 'L:\sigma_\omega', 'S:\sigma_\omega', 'Location', 'northeastoutside');
    end
    set(gca,'fontsize',16)
    
    
    subplot(1,3,3);
    semilogy(1:min(d,n), sigA, 'linewidth', 3); hold on;
    semilogy(1:min(d,n), abs((V'*x0)), 'linewidth', 0.5)
    semilogy(1:min(d,n), abs((U'*w(:, [D, 1]))./sigA), 'linewidth', 0.5);
    title('Eneryg of Iput and Noise'); xlabel('index i');
    if(numel(DEV) == 1)
        legend('\sigma_i', '|v_i^Tx_0|', '|u_i^T\omega/\sigma_i|');
    else
        legend('\sigma_i', '|v_i^Tx_0|', 'L: |u_i^T\omega/\sigma_i|', 'S: |u_i^T\omega/\sigma_i|', 'Location', 'northeastoutside');
        set(gca,'fontsize',16)
        
    end
end
fprintf('done %2.2f sec elapsed\n', toc);
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

function [ o ] = normA( A )
%UNT?TLED2 Summary of this function goes here
%   Detailed explanation goes here
[n,d]   = size(A);
if(n>d)
    A = A ./ repmat(sqrt(sum(A.^2, 1)), n,1);
    o=A.'*A;
else
    A = A ./ repmat(sqrt(sum(A.^2, 2)), 1,d);
    o=A*A.';
end
end