#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include "parameterset.h"
#include "distribution.h"
#include "parametermap.h"
#include "random.h"
#include "visualization.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "emulator.h"
#include <time.h>
#include <sys/stat.h>

using namespace std;

class ParameterSet;
class ParameterSetList;
class ProposalDistribution;
class EmulatorHandler;
class LikelihoodDistribution;
class PriorDistribution;
class VizHandler;

class MCMCConfiguration{
public:
	MCMCConfiguration(string run_file);
	MCMCConfiguration(string run_file, string configuration);
	~MCMCConfiguration();
	void Initialize();
	
	bool LOGLIKE;
	bool LOGPRIOR;
	bool LOGPROPOSAL;
	bool VIZTRACE;
	bool APPEND_TRACE;
	parameterMap parmap;
	string dir_name;
	string parameterfile;
	string configname;
	string parameter_file_name;
	vector<bool> LogParam;
	string EmulatorParams;
	vector<string> ParamNames;
	string Observables;
	vector<string> ObservablesNames;
	
	CRandom *randnum;
	LikelihoodDistribution *Likelihood;
	ProposalDistribution *Proposal;
	PriorDistribution *Prior;
	ParameterSetList *DummyList;
	double Max_Ranges[10],Min_Ranges[10]; // 
};

class MCMCRun{
public:
	parameterMap local_parmap;
	MCMCRun(MCMCConfiguration *mcmc_config);
	MCMCRun(MCMCConfiguration *mcmc_config, ParameterSet Theta0);
	~MCMCRun();
	double Run();
	
	ParameterSet *BestParameterSetPtr;
	vector<double> ParamValues;
	double bestlikelihood;
	double Likelihood_Current;
	
	int  WRITEOUT;
	int  MAXITERATIONS;
	bool RANDOM_THETA0;
	bool RESCALED_TRACE;
	
	string tracedir;
	MCMCConfiguration *mcmcconfig;
	VizHandler *Visualizer;
	ParameterSetList *ThetaList;
	int Accept_Count;
	int Viz_Count;
};

#endif