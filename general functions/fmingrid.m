function [min_par, min_val, fun_par, grid] = fmingrid(fun, xmin, xmax, num1, num2)

%% Parameter Grid 1. iteration
grid1       = logpsce(log10(xmin), log10(xmax), num1);
fun_par1    = zeros(1,num1);

for l = 1:num1
    fun_par1(l) = fun(grid1(l));
end
ind1        = argmin(fun_par1);

%% Parameter Grid 2. iteration
grid2       = logspace(log10(grid1(ind1))-1, log10(grid1(ind1))+1, num2);
fun_par2    = zeros(1,num2);

for l = 1:G
    fun_par2(l) = fun(grid2(l));
end
[min_val, ind2] = min(fun_par2);

%% solution
min_par = grid2(ind2);

%% 3. output
if(nargout > 2)
    %collect
    fun_par     = [fun_par1, fun_par2];
    grid        = [grid1, grid2];
    %sort
    [grid, P]   = sort(grid, 'ascend');
    fun_par     = fun_par(P);
end

end