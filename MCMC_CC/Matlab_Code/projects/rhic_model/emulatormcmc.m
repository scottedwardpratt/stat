clear all;
close all

global CONFIG;
global DATA;
global RANGES;
global LOGS;

Configure('set');

CONFIG.OPTIONS.MCMC.('MaxT')     = 500;
CONFIG.OPTIONS.MCMC.('WriteOut') = 100;
CONFIG.OPTIONS.MCMC.('OutDir') = 'cosmo_log_demo';
CONFIG.OPTIONS.MCMC.('ReCalc') = false;
CONFIG.OPTIONS.MCMC.('UseLogLike') = false;
CONFIG.FUNCTIONS.Likelihood = 'emulator_test_likelihood';
CONFIG.FUNCTIONS.Prior = 'emulator_prior';
CONFIG.FUNCTIONS.Proposal = 'emulator_proposal';
CONFIG.FUNCTIONS.PreMCMC = 'emulator_pre';
CONFIG.FUNCTIONS.PostMCMC = '';
CONFIG.OPTIONS.DEBUG.DebugFactor = false;
CONFIG.OPTIONS.DEBUG.ReportLike = false;
CONFIG.OPTIONS.DEBUG.DebugLikelihood = false;
CONFIG.OPTIONS.DEBUG.HandleMessages = false;

THETA0 = struct('HYDRO_T0', 0.674776,...
                'GLAUBER_WNBIN_RATIO', 0.936926,...
                'GLAUBER_K_TAU', 2.45906,...
                'HYDRO_INIT_NS', 1.52115,...
                'HYDRO_INIT_FLOW', 1.03372,...
                'HYDRO_SVRATIO', 0.207722,...
                'EQOFST_BUMP_HEIGHT', 0.493751);
            
ACTUAL = struct('HYDRO_T0', (1.1+0.6)/2,...
                'GLAUBER_WNBIN_RATIO', 0.8,...
                'GLAUBER_K_TAU', (2.4+3.3)/2,...
                'HYDRO_INIT_NS', 1.25,...
                'HYDRO_INIT_FLOW', .8,...
                'HYDRO_SVRATIO', 0.14,...
                'EQOFST_BUMP_HEIGHT', 0.5);
DATA = [];

RANGES = struct();
            
CONFIG.STATIC = struct();

F = fieldnames(THETA0);
VALS = num2cell(repmat(0.01,size(F)));

CONFIG.PARAMETERS = [F',fieldnames(CONFIG.STATIC)'];
CONFIG.PROPOSAL.MIXING_STDDEV = cell2struct(VALS,F);
%finally, check to make sure your changes were valid
Configure('check');

[DATA, RANGES] = feval(CONFIG.FUNCTIONS.PreMCMC,ACTUAL);

AUG = [];
%Perform the search!
Theta = mcmc(THETA0,AUG);

feval(CONFIG.FUNCTIONS.PostMCMC);


