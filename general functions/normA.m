function [ o ] = normA( A, PLOT)
%normA shows correlation between columns of A
%   [ o ] = normA( A, PLOT)
%   if want to plot make PLOT = true
%
A = A ./ repmat(sqrt(sum(A.^2, 1)), size(A,1),1);
o=A'*A;
if(PLOT)
    figure;
    imagesc(o, [0, 1]);
    colorbar
    title('coherence of A: $<\hat{a_i}, \hat{a_j}>$', 'Interpreter', 'latex');
    xlabel('column index i'); ylabel('column index j')
    axis square
    set(gca,'fontsize',16)
end
end