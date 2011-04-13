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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

class ParameterSet; class ParameterSetList; class ProposalDistribution;

class MCMC{
public:
	parameterMap parmap;
	ParameterSetList *ThetaList;
	MCMC(string run_file);
	~MCMC();
	void Run();
	
private:
	int MAXITERATIONS, WRITEOUT, Accept_Count;
	string dir_name, parameter_file_name;
	
	CRandom *randnum;
	bool LOGLIKE;
	bool LOGPRIOR;
	bool LOGPROPOSAL;
	

	// LikelihoodDistribution *Likelihood;
	ProposalDistribution *Proposal;
	// PriorDistribution *Prior;
};

#endif