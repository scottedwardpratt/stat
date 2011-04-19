#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "mcmc.h" 
#include "parametermap.h"

using namespace std;

class EmulatorHandler;
class MCMC;
class ParameterSet;

class Distribution{
public:
	Distribution(MCMC *mcmc_in);
	Distribution();
	
protected:
	MCMC * mcmc;
	bool SepMap;
	bool TIMING;
	parameterMap * parmap;
};

class ProposalDistribution:private Distribution {
public:
	ProposalDistribution(MCMC *mcmc_in);
	ParameterSet Iterate();
	double Evaluate(ParameterSet Theta);
private:
	bool SymmetricProposal;
	vector<double> MixingStdDev;
	double *Min_Ranges;
	double *Max_Ranges;
	gsl_rng * randy;
	int FindParam(string param_name);
};

class PriorDistribution:private Distribution {
public:
	PriorDistribution(MCMC *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class LikelihoodDistribution:private Distribution {
public:
	LikelihoodDistribution(MCMC *mcmc_in);
	~LikelihoodDistribution();
	double Evaluate(ParameterSet Theta);
private:
	vector<double> GetData();
	vector<double> DATA;
	bool UseEmulator;
	EmulatorHandler * emulator;
};
#endif