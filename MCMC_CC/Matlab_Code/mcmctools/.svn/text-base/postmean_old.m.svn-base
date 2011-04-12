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
if (nargin == 2)
	I = varargin{1};
else
	Iempty = logical(1);
end

if (prod(size(I)) == 1)
	Iscalar = logical(1);
end

S = dir(OPTS.OutDir);

if (length(S) == 0)
	error('Can''t get guesses -- no OUT files in OutDir.');
end

N = strvcat(S.name);
ind_out = find(strcmp(cellstr(N(:,1:4)),'out_'));
%Get the numbers...
[N,ind] = sort(str2num(N(ind_out,5:end-4)));
ind_out = ind_out(ind);

t = 0;
N = 0;
SX = 0;
SXsq = 0;



for i = 1:length(ind_out)
	X = getfield(load(fullfile(OPTS.OutDir,S(ind_out(i)).name),'X'),'X');
	t = [t(end)+1:t(end)+size(X,1)];
	if ~Iempty
		if Iscalar
			ind = (t>I);
		else
			ind = ismember(t,I);
			I = setdiff(I,t(ind));
		end
	else %Use all samples.
		ind = ones(1,length(t));
	end
	
	N = N + nnz(ind);
	SX = SX + sum(X(ind,:),1);
	
end

OUT = SX./N;
