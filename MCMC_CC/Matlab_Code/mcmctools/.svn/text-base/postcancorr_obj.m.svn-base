function OUT = postcancorr_obj(X,LOGBF,t,FILTER,FLAG,varargin)
%POSTCANCORR_OBJ Objective function for POSTCANCORR.
%

persistent LOGBF_TOT X_TOT;


if ~isempty(FILTER)
	X = feval(FILTER,X,t,varargin{:});
end

if isstruct(X)
	VALS = struct2cell(X);
	VALS = cell2mat(VALS)';
	X_TOT = cat(1,X_TOT,VALS);
	LOGBF_TOT = cat(1,LOGBF_TOT,LOGBF);
else
	X_TOT = cat(1,X_TOT,X);
	LOGBF_TOT = cat(1,LOGBF_TOT,LOGBF);
end


if strcmp(FLAG,'report')
	ind = ~isnan(LOGBF_TOT);
	X_TOT = X_TOT(ind,:);
	LOGBF_TOT = LOGBF_TOT(ind,:);
	
	
	if (size(X_TOT,1) <= size(X_TOT,2))
		warning('on');	
		warning('Number of usable data points less than number of parameters.');
		warning('off');
	end
	
	[A,B] = canoncorr(X_TOT,LOGBF_TOT);
	OUT = A;
else
	OUT = 0;
end
