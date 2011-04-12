function OUT = postmean(OPTS,varargin)
%POSTMEAN Calculate mean of samples from posterior distribution.
%  POSTMEAN(OPTS) mean is taken over all samples in OPTS.OutDir
%  POSTMEAN(OPTS,N) mean is taken over the last N samples in OPTS.OUTDIR.
%  POSTMEAN(OPTS,I) mean is taken over SAMPLE(I) where I is a vector of 
%  indices in the range 1 up to total number of samples in OPTS.OutDir.
%


Iempty = logical(0);
Iscalar = logical(0);

I = [];
if (nargin >= 2)
	I = varargin{1};
else
	Iempty = logical(1);
end

FILTER = [];
PARAMS = {};

if (nargin >= 3)
	FILTER = varargin{2};
	PARAMS = {varargin{3:end}};
end

[SX,NX] =  postfun(OPTS,I,'postmean_obj',FILTER,PARAMS{:});

OUT = SX/NX;
