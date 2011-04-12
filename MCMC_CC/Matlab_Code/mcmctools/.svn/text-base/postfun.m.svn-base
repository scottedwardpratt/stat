function OUT = postfun(OPTS,I,FUN,varargin)
%POSTFUN Calculate a function over the MCMC chain stored in OPTS.OutDir.
%   POSTFUN(OPTS,N,FUN,P1,P2,...)
%   FUN should having calling syntax FUN(X,t,FLAG,P1,P2,...)
%   where X is a matrix of parameters (one column per parameter) and each row
%   corresponds to iteration t of the MCMC chain.  FLAG is 'start' when the
%   first file is read, '' while intermediate files are being read 
%   FLAG = 'report' when the last file is read.  The value of
%   FLAG can be used to determine when to carry out certain calculations
%   or may be ignored.
%   
%   The reason this is called POSTFUN is that it is tyically evoked after
%   the completion of an MCMC simulation and is used to examine statistically
%   or graphically the parameter chains stored in the OutDir.
%   
%

Iscalar = logical(0);
Iempty = isempty(I);


if (prod(size(I)) == 1)
	Iscalar = logical(1);
end

S = dir(fullfile(OPTS.OutDir,'out_*.mat'));

if (length(S) == 0)
	error('Can''t get guesses -- no OUT files in OutDir.');
end

N = strvcat(S.name);

if isempty(N) 
	error('Directory is empty.');
elseif (size(N,2) < 4)
	error('No output files in directory.');
end

% ind_out = find(strcmp(cellstr(N(:,1:4)),'out_*.mat'));
% 
% %Get the numbers...
% if isempty(ind_out)
% 	error('No output files in directory.');
% end
ind_out = [1:length(S)];
[N,ind] = sort(str2num(N(:,5:end-4)));


ind_out = ind_out(ind);

t = 0;
N = 0;

FILTER = '';

if (nargin >= 4)
	FILTER = varargin{1};
	PARAMS = {varargin{2:end}};
end


for i = 1:length(ind_out)
	if (i == length(ind_out))
		FLAG = 'report';
	elseif (i == 1)
		FLAG = 'start';
	else
		FLAG = '';
	end
	
	X = getfield(load(fullfile(OPTS.OutDir,S(ind_out(i)).name),'X'),'X');
	LOGBF = getfield(load(fullfile(OPTS.OutDir,S(ind_out(i)).name),'LOGBF'),'LOGBF');
	t = [t(end)+1:t(end)+size(X,1)];
	
	if ~Iempty
		if Iscalar
			ind = (t>I);
		else
			ind = ismember(t,I);
			I = setdiff(I,t(ind));
		end
	else %Use all samples.
		ind = logical(ones(1,size(X,1)));
	end
	
	N = N + nnz(ind);
	OUT = feval(FUN,X(ind,:),LOGBF(ind),t(ind),FILTER,FLAG,PARAMS{:});
end


clear(FUN);
