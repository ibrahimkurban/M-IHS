function plotf(x, y, varargin)
figure;
grid on;
hold on

if(exist('y', 'var'))
    plot(x,y);
else
    plot(x)
end

if(size(varargin,2) == 1)
    title(varargin{1,1})
elseif((size(varargin,2) == 2))
    xlabel(varargin{1,1});
    ylabel(varargin{1, 2});
elseif((size(varargin,2) == 3))
    title(varargin{1,1})
    xlabel(varargin{1,2});
    ylabel(varargin{1,3});
end
