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
//#include "emulator.h"
#include "quad.h"
#include <time.h>
#include <sys/stat.h>
#include "EmuPlusPlus.h"

using namespace std;

class ParameterSet;
class ParameterSetList;
class ProposalDistribution;
class EmulatorHandler;
class LikelihoodDistribution;
class PriorDistribution;
class VizHandler;

class MCMC{
public:
	MCMC(string run_file, string configuration);
	~MCMC();
	void GetRanges(); //This loads the names and ranges of the emulator parameters
	void FirstPass(); //This prepares the objects necessary for the MCMC Run
	void LoadObservables(); //This loads the names of the emulator observables
	void Run(); //This performs the MAXITERATIONS sampling

	bool LOGLIKE;
	bool LOGPRIOR;
	bool LOGPROPOSAL;
	bool CREATE_TRACE;
	bool APPEND_TRACE;
	bool SUPPRESS_ERRORS;
	bool RESCALED_TRACE;
	bool PRESCALED_PARAMS;
	int  WRITEOUT;
	int  BURN_IN;
	int  MAXITERATIONS;
	bool RANDOM_THETA0;
	bool VIZTRACE;
	bool QUIET;
	bool WAIT_TO_PRINT_DENSITY_PLOTS;
	bool MOVING_WINDOW;
	bool DENSITY_PLOTS;
	string MODEL;
	parameterMap parmap;
	string dir_name;
	string parameterfile;
	string configname;
	string parameter_file_name;
	vector<bool> LogParam;
	string EmulatorParams;
	vector<string> ParamNames;
	vector<string> ObservablesNames;
	
	CRandom *randnum;
	LikelihoodDistribution *Likelihood;
	ProposalDistribution *Proposal;
	PriorDistribution *Prior;
	double Max_Ranges[30],Min_Ranges[30];
	
	//ParameterSet *BestParameterSetPtr;
	vector<double> BestParameterSet;
	vector<double> ParamValues;
	vector<double> Theta;
	vector<double> Proposed_Theta;
	vector<double> Scaled_Theta; //These are scaled to [0,1] for the trace density plots
	double bestlikelihood;
	double Likelihood_Current;
	
	string tracedir;
	VizHandler *Visualizer;
	int WriteOutCounter;
	//ParameterSetList *ThetaList;
	vector<vector<double>> ThetaList;
	vector<vector<double>> Scaled_ThetaList;
	int Accept_Count;
	int Viz_Count;
};

#endif