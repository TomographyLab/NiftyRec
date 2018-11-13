 function data = poisson(xm, seed)
%|function data = poisson(xm, seed, [options])
%|
%| option
%|	'factor'	see poisson2.m
%|
%| Generate Poisson random vector with mean xm.
%| For small, use poisson1.m
%| For large, use poisson2.m
%| see num. rec. C, P. 222
%|
%| Copyright 1997-4-29, Jeff Fessler, University of Michigan

%if nargin < 1, help(mfilename), error(mfilename), end
%if nargin == 1 & streq(xm, 'test'), poisson_test, return, end

if ~exist('seed','var')
	seed = [];
end

arg.factor = 0.85;
%arg = vararg_pair(arg, varargin);

data	= xm;
xm	= xm(:);

small = xm < 12;
data( small) = poisson1(xm( small), seed);
data(~small) = poisson2(xm(~small), 'seed', seed, 'factor', arg.factor);


