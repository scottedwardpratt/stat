function [THETA0] = cosmo_pre(THETA0, ACTUAL)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global CONFIG;

if(CONFIG.OPTIONS.DEBUG.HandleMessages)
    HandleMessages('iterate');
end

%Optimize the log of those that we want to be non-negative...
for i=1:length(CONFIG.LOGPARAM)
	THETA0.(CONFIG.LOGPARAM{i}) = log(THETA0.(CONFIG.LOGPARAM{i}));
end


%%%%%% SYNTHETIC DATA %%%%%%%%%
CMD = '/usr/local/bin/cosmosurvey';
F_ACTUAL = fieldnames(ACTUAL);
F_STATIC = fieldnames(CONFIG.STATIC);
for i = 1:length(F_ACTUAL)
    CMD = [CMD, ' -' F_ACTUAL{i} ' ' num2str(ACTUAL.(F_ACTUAL{i}), '%.15f')];
end
for j = 1:length(F_STATIC)
    if(any(strcmp(F_STATIC{j},{'nz','nl'})) == 1)
        CMD = [CMD, ' -' F_STATIC{j} ' ' num2str(CONFIG.STATIC.(F_STATIC{j}), '%f')];
    else
        CMD = [CMD, ' -' F_STATIC{j} ' ' num2str(CONFIG.STATIC.(F_STATIC{j}), '%.15f')];
    end
end
CMD = [CMD ' > ' CONFIG.LIKELIHOOD.DATAFILE];

[MSG,result] = system(CMD);
if(MSG ~= 0)
    disp(result);
    error('ERROR');
end;

NoisyData(CONFIG.LIKELIHOOD.DATAFILE, CONFIG.STATIC.nl, CONFIG.STATIC.nz, .1, 1000);

%%%%%%%%%%
%Add in warnings to ensure valid statistics at the end of the trace
if(CONFIG.OPTIONS.STATS.statsamplesize > CONFIG.OPTIONS.MCMC.MaxT/2)
    warnmsg = 'Statistics sample size is more than half of the total iterations. Consider resizing.';
    if(CONFIG.OPTIONS.DEBUG.HandleMessages)
        HandleMessages(4, warnmsg);
    end
    warning(warnmsg);
end
if(CONFIG.OPTIONS.STATS.statsamplesize > CONFIG.OPTIONS.MCMC.MaxT)
    errmsg = 'Statistics sample size is larger than total number of iterations.';
    if(CONFIG.OPTIONS.DEBUG.HandleMessages)
        HandleMessages(3, errmsg);
    end
    error(errmsg);
end
%%%%%%%%%%
%Testing initial parameters

NEWPARAMS =THETA0;

for i=1:length(CONFIG.LOGPARAM)
	NEWPARAMS.(CONFIG.LOGPARAM{i}) = exp(NEWPARAMS.(CONFIG.LOGPARAM{i}));
end

N= CONFIG.STATIC.nz*CONFIG.STATIC.nl; %Total number of data points

V = zeros(1,N);


for k = 1:length(CONFIG.CHOL)
    
    if(mod(str2num(CONFIG.CHOL{k}(2)),1)~= 0)
        if(CONFIG.OPTIONS.DEBUG.HandleMessages)
            HandleMessages(3, 'Error interpreting cholesky paramters.',['-' CONFIG.CHOL{k}]);
        end
        error('Error interpreting cholesky parameters.');
    end
    index = str2num(CONFIG.CHOL{k}(2));
    
    V(index) = NEWPARAMS.(CONFIG.CHOL{k});
end
        
TempL = BuildCholeskyFinal(V);

clear TempL; %since the factorization is stored in the queue, don't even bother using up the memory

end

