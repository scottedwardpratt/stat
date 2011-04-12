clear all;
close all;

global CONFIG;
global LOGS;

%disp('Starting MCMC search....');

%Initialize the CONFIG variable to default values
Configure('set');

%...then change things specific to your run.
CONFIG.OPTIONS.MCMC.('MaxT')     = 500;
CONFIG.OPTIONS.MCMC.('WriteOut') = 100;
CONFIG.OPTIONS.MCMC.('OutDir') = 'cosmo_log_demo';
CONFIG.OPTIONS.MCMC.('ReCalc') = false;
CONFIG.OPTIONS.MCMC.('UseLogLike') = false;
CONFIG.OPTIONS.LIKELIHOOD.DetMethod = '';
CONFIG.OPTIONS.LIKELIHOOD.Qsize = 3;
CONFIG.OPTIONS.LIKELIHOOD.Cusping = true;
CONFIG.OPTIONS.STATS.statsamplesize = 1;
CONFIG.FUNCTIONS.Likelihood = 'cosmoLikeLihood';
CONFIG.FUNCTIONS.Prior = 'cosmo_prior';
CONFIG.FUNCTIONS.Proposal = 'cosmo_proposal';
CONFIG.FUNCTIONS.PreMCMC = 'cosmo_pre';
CONFIG.FUNCTIONS.PostMCMC = 'cosmo_post';
CONFIG.CHOL = {'x11'};
CONFIG.LIKELIHOOD.DATAFILE = 'MCMCsyntheticdata.d';
CONFIG.OPTIONS.DEBUG.DebugFactor = false;
CONFIG.OPTIONS.DEBUG.ReportLike = false;
CONFIG.OPTIONS.DEBUG.DebugLikelihood = false;
			
THETA0 = struct('w',-1.0,...
                's8',1.8,...
                'oo',0.30,...
                'g',0.15,...
                'x11', 5);
            
ACTUAL = struct('w',-1.0,...
                's8',0.8,...
                'oo',0.3,...
                'g',0.550);
            
CONFIG.STATIC = struct('z2',2.5,...
                'fl',1.6e-13,...
                'do',12.70393163,...
                'z1',0.01,...
                'l2',46.0,...
                'l1',44.0,...
                'nz',10,...
                'nl',10);

F = fieldnames(THETA0);
VALS = num2cell(repmat(0.01,size(F)));


CONFIG.PARAMETERS = [F',fieldnames(CONFIG.STATIC)'];
CONFIG.PROPOSAL.MIXING_STDDEV = cell2struct(VALS,F);
CONFIG.LOGPARAM = setdiff(F,{'w', 'x21'});

%finally, check to make sure your changes were valid
Configure('check');

THETA0 = feval(CONFIG.FUNCTIONS.PreMCMC, THETA0, ACTUAL);

AUG = [];
%Perform the search!
Theta = mcmc(THETA0,AUG);

feval(CONFIG.FUNCTIONS.PostMCMC);
