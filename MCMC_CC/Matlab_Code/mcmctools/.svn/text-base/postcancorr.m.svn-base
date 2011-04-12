function OUT = postcancorr(OPTS,varargin)
%POSTCANCORR Calculate canonical correlation between log Bayes and parameters.
%
%


I = [];

if (nargin >= 2)
	I = varargin{1};
end

FILTER = [];
PARAMS = {};

if (nargin >= 3)
	FILTER = varargin{2};
	PARAMS = {varargin{3:end}};
end

OUT =  postfun(OPTS,I,'postcancorr_obj',FILTER,PARAMS{:});

