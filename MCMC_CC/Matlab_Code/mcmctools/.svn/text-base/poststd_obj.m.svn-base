function OUT = poststd_obj(X,LOGBF,t,FILTER,FLAG,varargin)
%POSTSTD_OBJ Objective function for POSTSTD.
%

persistent SUMX SUMXsq N;


if isempty(SUMX)
	SUMX = 0;
end

if isempty(SUMXsq)
	SUMXsq = 0;
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
	SUMXsq = SUMXsq + sum(VALS.^2,1);
	N = N + size(VALS,1);
else
	N = N + size(X,1);
	SUMX = SUMX + sum(X,1);
	SUMXsq = SUMXsq + sum(X.^2,1);
end


if strcmp(FLAG,'report')
	OUT = (SUMXsq - ((SUMX.^2)./N))./(N-1);
	OUT = sqrt(OUT);
else
	OUT = 0;
end


