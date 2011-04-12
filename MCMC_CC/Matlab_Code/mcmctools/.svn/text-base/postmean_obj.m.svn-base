function OUT = postmean_obj(X,LOGBF,t,FILTER,FLAG,varargin)
%POSTMEAN_OBJ Objective function for POSTMEAN.
%

persistent SUMX N;

if isempty(SUMX)
	SUMX = 0;
end

if isempty(N)
	N = 0;
end

if ~isempty(FILTER)
	X = feval(FILTER,X,t,varargin{:});
end


if isstruct(X)
	VALS = struct2cell(X);
	VALS = cell2mat(VALS)';
	SUMX = SUMX + sum(VALS,1);
	N = N + size(VALS,1);
else
	SUMX = SUMX + sum(X,1);
	N = N + size(X,1);
end


if strcmp(FLAG,'report')
	OUT = SUMX./N;
else
	OUT = 0;
end
