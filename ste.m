function [ y ] = ste( x,dim )
% ste Computes the standard error of a set of values. If x is a vector,
% y is the ste of x. If x is a matrix, y is a raw vector with the ste
% of each column of x.

if sum(isnan(x(:))) > 0,
    warning('Data contains NaN values. Use nanste to ignore them');
end
if nargin == 1,
    dim = 1;
end
if size(x,1) == 1 && size(x,2) > 1,
    x = x';
end

y = std(x,0,dim)./sqrt(size(x,dim));

if isempty(y)
    y = NaN;
end
end

