function C = supervertcat(varargin)
% C = supervertcat(A1, A2, ...)

[i j v] = cellfun(@find, varargin, 'uni',0);
r = [0 cumsum(cellfun('size',varargin,1))];
[i j v] = arrayfun(@(k) deal(r(k)+i{k}(:),j{k}(:),v{k}(:)), 1:nargin,'uni',0);
C = accumarray([cat(1,i{:}) cat(1,j{:})], cat(1,v{:}));
end