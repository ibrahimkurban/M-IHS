function [min_par, min_val, fun_par, grid] = grid_search(fun, grid1, num_grid2)
%%GRID_SEARCH finds the minimum point of a function which is entered as
%%function handle. Function is evaluated firslty over grid points given as
%%input. Then at minimum point a new grid is defined around min point.
%%Afterward min point is updated as the minimum of the function over
%%second grid.
%
%   Author  : Ibrahim Kurban Ozaslan
%   Date    : 14.07.2018
%   v1      : only 2 iterate, 2.grid is around linspace(min/10, min*10, num_grid2)
%
%   [min_par, min_val, fun_par, grid] = grid_search(fun, grid1, num_grid2)
%
%   fun         : fun object @ paramter i.e. grid
%   grid1       : given grid
%   num_grid2   : # of grid point in grid 2
%
%% Parameter Grid 1. iteration
G           = numel(grid1);
fun_par1    = zeros(1,G);

for l = 1:G
    fun_par1(l) = fun(grid1(l));
end
[min_val, ind1] = min(fun_par1);

if(nargin> 2)
%% Parameter Grid 2. iteration
grid2       = linspace(grid1(ind1)/10, grid1(ind1)*10, num_grid2);
G           = num_grid2;
fun_par2    = zeros(1,G);

for l = 1:G
    fun_par2(l) = fun(grid2(l));
end
[min_val, ind2] = min(fun_par2);

%% solution
min_par = grid2(ind2);

else
    fun_par2 = [];
    grid2   = [];
    min_par     = grid1(ind1);
end
%% 3. output
if(nargout > 2)
    %collect
    fun_par = [fun_par1, fun_par2];
    grid = [grid1, grid2];
    %sort
    [grid, P] = sort(grid, 'ascend');
    fun_par = fun_par(P);
end