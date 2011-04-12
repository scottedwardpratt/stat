function OUT = postmean(OPTS,varargin)
%POSTMEAN Calculate mean of samples from posterior distribution.
%  POSTMEAN(OPTS) mean is taken over all samples in OPTS.OutDir
%  POSTMEAN(OPTS,N) mean is taken over the last N samples in OPTS.OUTDIR.
%  POSTMEAN(OPTS,I) mean is taken over SAMPLE(I) where I is a vector of 
%  indices in the range 1 up to total number of samples in OPTS.OutDir.
%  POSTMEAN(OPTS,I,FILTER,P1,P2,...) sends data through FILTER prior to 
%  including in calculation of mean.  FILTER is a user-supplied m-file 
%  with the calling syntax FILTER(X,t,P1,P2,...) where the ith column of X 
%  corresponds to the ith parameter t is a vector indicating the iteration 
%  corresponding to the row in X and P1,P2,... are optional user-supplied 
%  parameters.
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

OUT =  postfun(OPTS,I,'postmean_obj',FILTER,PARAMS{:});

