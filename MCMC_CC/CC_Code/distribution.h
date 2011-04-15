#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "mcmc.h" 

using namespace std;

class EmulatorHandler;

class Distribution{
public:
	MCMC * mcmc;
	Distribution(MCMC *mcmc_in);
	Distribution();
	
	virtual double Evaluate(ParameterSet Theta) {return(0);}
private:
	
};

class ProposalDistribution: public Distribution {
public:
	ProposalDistribution(MCMC *mcmc_in);
	ParameterSet Iterate();
	double Evaluate(ParameterSet Theta);
private:
	bool SymmetricProposal;
	vector<double> MixingStdDev;
	double *Ranges[2];
	gsl_rng *randy;
	int FindParam(string param_name);
};

class PriorDistribution:public Distribution {
public:
	PriorDistribution(MCMC *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class LikelihoodDistribution:public Distribution {
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