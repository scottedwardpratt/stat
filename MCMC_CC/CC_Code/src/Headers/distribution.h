#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "mcmc.h" 
#include "parametermap.h"

using namespace std;

class EmulatorHandler;
class MCMCConfiguration;
class MCMCRun;
class ParameterSet;

class Distribution{
public:
	Distribution(MCMCConfiguration *mcmc_in);
	Distribution();

protected:
	MCMCConfiguration * mcmc;
	bool SepMap;
	bool TIMING;
	bool VERBOSE;
	bool DEBUG;
	parameterMap * parmap;
	
	double Normal(double x, double mu, double sigma);
	double Gaussian(double x, double mu, double sigma);
	double Log_MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double LogNormal(double x, double mu, double sigma);
};

class ProposalDistribution:private Distribution {
public:
	ProposalDistribution(MCMCConfiguration *mcmc_in);
	ParameterSet Iterate(ParameterSet current);
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
	PriorDistribution(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class LikelihoodDistribution:private Distribution {
public:
	LikelihoodDistribution(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution();
	double Evaluate(ParameterSet Theta);
private:
	vector<double> GetData();
	vector<double> DATA;
	bool UseEmulator;
	ofstream emulator_test;
	EmulatorHandler * emulator;
};
#endif