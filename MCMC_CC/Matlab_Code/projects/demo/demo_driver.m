clear all;
close all;

global CONFIG;
global LOGS;
Configure('set');

CONFIG.LIKELIHOOD.DATAFILE = 'DemoMCMCDatafile.d';
CONFIG.PRE.NumDataPoints = 100;

ACTUAL = struct('mean',5,'sigma',1);
THETA0 = struct('mean',3,'sigma',1);

CONFIG.PARAMETERS = {'mean','sigma'};
CONFIG.LOGPARAM = {'sigma'};

VALS = num2cell(repmat(0.01,size(CONFIG.PARAMETERS)));
STDDEV = cell2struct(VALS,CONFIG.PARAMETERS);
CONFIG.PROPOSAL.MIXING_STDDEV = STDDEV;

CONFIG.FUNCTIONS.PreMCMC = 'Demo_PreMCMC';
CONFIG.FUNCTIONS.Prior = 'demo_prior';
CONFIG.FUNCTIONS.Proposal = 'demo_prop';

Configure('check');


THETA0 = feval(CONFIG.FUNCTIONS.PreMCMC, ACTUAL, THETA0);

AUG = [];
Theta = mcmc(THETA0, AUG);
feval(CONFIG.FUNCTIONS.PostMCMC);