function [THETA0] = Demo_PreMCMC(ACTUAL, THETA0)

global CONFIG;
global LOGS;

for i=1:length(CONFIG.LOGPARAM)
	THETA0.(CONFIG.LOGPARAM{i}) = log(THETA0.(CONFIG.LOGPARAM{i}));
end

fid = fopen(CONFIG.LIKELIHOOD.DATAFILE, 'w');

DATA = normrnd(ACTUAL.mean,ACTUAL.sigma,CONFIG.PRE.NumDataPoints,1);

fprintf(fid, '%d\n', DATA);

fclose(fid);

end

