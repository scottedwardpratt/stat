#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include "parameter.h"
#include "distribution.h"
#include "parametermap.h"
#include "random.h"
#include "visualization.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "emulator.h"
#include <time.h>

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
	parameterMap parmap;
	ParameterSetList *ThetaList;
	MCMC(string run_file);
	~MCMC();
	void Run();
	string dir_name;
	string runnickname;
	bool LOGLIKE;
	bool LOGPRIOR;
	bool LOGPROPOSAL;
	bool VIZTRACE;
	int  WRITEOUT;
	//Made public to test against Matlab, change back later.
	LikelihoodDistribution *Likelihood;
private:
	int MAXITERATIONS, Accept_Count, Viz_Count;
	string parameter_file_name;
	CRandom *randnum;
	ProposalDistribution *Proposal;
	PriorDistribution *Prior;
	VizHandler *Visualizer;
};

#endif