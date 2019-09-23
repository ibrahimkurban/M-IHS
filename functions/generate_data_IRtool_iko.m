function [A, b, x1, x0, err_x0, dev, prob_info, err_x,k0, par0, xxOR, U, sig, V] = generate_data_IRtool_iko(n, problem, example, level, svd_analysis, plotting, varargin)
%%GENERATE_DATA_IRTOOL generates dataset by using IR tool
%[A, b, x1, x0, err_x, dev, prob_info, k0, par0, xxOR, U, sig, V] = ...
%       generate_data_IRtool(n, problem, example, level, svd_analysis, plotting, varargin)
%
% problem       : 1,2,3,4
% example       : 1,2,3....p; p is dependent ob problem
% level         : double, noise level |w|/|Ax0|
% svd_anaylsis  : boolean
% plotting      : boolean, requires svd_analysis
%
%
% VARARGIN
% 'type'        : double problem type
% 'level'       : illposedness severity level
%
% PROBLEMS
%   1 : BLURRING
%           example : 1,...,9
%           'type', : 1,...,6
%           'level' : 1,2,3 (blurring level)
%   2 : SEISMIC
%           example : 1,...,8
%           'type', : 1,2
%           'level' : 1,2,3... makes overdetermined
%   3 : X-RAY TOMOGRAPHIC
%           example : 1,...,8
%           'type', : 1,2
%           no level input
%   4 : SPHERICAL
%           example : 1,...,8
%           no type  input
%           no level input
%
%

%% generate data
fprintf('\nData generation...\n'); tic;
[A, b0, x0, prob_info] = generate_IRtool_iko(n, problem, example, varargin{:});
err_x0 = @(xx)(sqrt(sum((x0- xx).^2, 1))/norm(x0));

%initial guess
[n,d]       = size(A);
x1          = zeros(d,1);

% %% DELETE DELETE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
% As          = full(A);
% [U,S,V]     = svd(As,'econ');
% x0          = V*ones(d,1);
% b0          = A*x0;
%% Add noise

w       = randn(n,numel(level));
dev     = level*norm(b0)./sqrt(sum(w.^2,1));
w       = dev(:)'.*w;
b       = b0 + w;





fprintf('%2.1e sec elapsed\n\n', toc);
%% svd analysis
if(svd_analysis == true)
    fprintf('SVD Analysis...\n'); tic;
    if(issparse(A))
        As = full(A);
    end
    %SVD
    [U,S,V]     = svd(As,'econ');
    sig         = diag(S);
    
    N           = numel(level);
    
    %%% TSVD
    k0          = zeros(N,1);
    for i = 1:N
        [xxOR{1}(:,i), k0(i)]   = LS_truncated_iko(U,sig,V, b(:,i), x1, x0 );
    end
    if(N ==1)
        err_x = @(xx)(sqrt(sum((xxOR{1} - xx).^2, 1))/norm(xxOR{1}));
    else
        err_x = @(xx, dev)(sqrt(sum((xxOR{1}(:,dev) - xx).^2, 1))/norm(xxOR{1}(:,dev)));
    end
    
    %%% Tikhonov
    par0        = zeros(N,1);
    if(N==1)
        [xxOR{2}, par0]   = LS_ridge_iko(U,sig,V, b, x1, @(x)err_x(x));
        fprintf('Noise deviaiton     : %1.2e\n', dev);
        fprintf('Effective rank      : %d\n', k0);
        fprintf('Oracle reg. par     : %1.1e\n', par0);
        fprintf('Oracle Error (log)  : %1.1f,    (dec) : %1.3f\n', log10(err_x(xxOR{2})), err_x(xxOR{2}));
    else
        for i =1:N
            [xxOR{2}(:,i), par0(i)]   = LS_ridge_iko(U,sig,V, b(:,i), x1, @(x)err_x(x, i));
        end
    end
    
    fprintf('\n%2.1e sec elapsed\n\n', toc);
end


%% PLOT
if(~svd_analysis && plotting)
    fprintf('!!!!!!!!!! PLOTTING REQURE SVD ANALYSIS !!!!!!!!!!')
elseif(svd_analysis && plotting)
    N = length(dev);
    if N ==1
        N = [];
    end
    figure('Name', 'Theoric Plot');
    
    %     subplot(1,3,1);
    %     %coherence
    %     imagesc(normA(A), [0 1]); colorbar;
    %     title('coherence of A: $<\hat{a_i}, \hat{a_j}>$', 'Interpreter', 'latex');
    %     xlabel('column index i'); ylabel('column index j')
    %     axis square
    %     set(gca,'fontsize',16)
    
    %singular values
    subplot(1,2,1);
    plot(log10(sig), 'linewidth', 3); hold on;
    plot(log10(abs(U'*b(:,[N, 1]))), 'linewidth', .5); hold on;
    plot(log10(dev([N, 1]).*ones(min(n,d),1)), 'linewidth', .5); hold on;
    title('Energy of Measurement'); xlabel('index i');
    if(numel(dev) == 1)
        legend('\sigma_i', '|u_i^Ty|', '\sigma_\omega');
    else
        legend('\sigma_i', 'L:|u_i^Ty|', 'S:|u_i^Ty|', 'L:\sigma_\omega', 'S:\sigma_\omega', 'Location', 'northeastoutside');
    end
    set(gca,'fontsize',16)
    grid on;
    
    subplot(1,2,2);
    plot(log10(sig), 'linewidth', 3); hold on;
    plot(log10(abs((V'*x0))), 'linewidth', .5)
    plot(log10(abs((U'*w(:, [N, 1]))./sig)), 'linewidth', .5);
    title('Eneryg of Iput and Noise'); xlabel('index i');
    if(numel(dev) == 1)
        legend('\sigma_i', '|v_i^Tx_0|', '|u_i^T\omega/\sigma_i|');
    else
        legend('\sigma_i', '|v_i^Tx_0|', 'L: |u_i^T\omega/\sigma_i|', 'S: |u_i^T\omega/\sigma_i|', 'Location', 'northeastoutside');
    end
    set(gca,'fontsize',16)
    grid on;
    
end

end

function [ o ] = normA( A )
%normA calculate correlation between columns or rows
[n,d]   = size(A);
if(n>=d)
    A = A ./ repmat(sqrt(sum(A.^2, 1)), n,1);
    o=A.'*A;
else
    A = A ./ repmat(sqrt(sum(A.^2, 2)), 1,d);
    o=A*A.';
end
end




