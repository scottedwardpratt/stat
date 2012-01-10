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
	Distribution();
	virtual ~Distribution();

	//protected:
	MCMCConfiguration * mcmc;
	gsl_rng * randy;
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

/** ---------------------------------------- */

class ProposalDistribution:public Distribution {
public:
	ProposalDistribution(MCMCConfiguration *mcmc_in);
	ParameterSet Iterate(ParameterSet current);
	virtual double Evaluate(ParameterSet Theta);
	//protected:
	bool SymmetricProposal;
	bool Rescaled_Method;
	vector<double> MixingStdDev;
	double *Min_Ranges;
	double *Max_Ranges;
	int FindParam(string param_name);
};

class LikelihoodDistribution:public Distribution{
public:
	LikelihoodDistribution();
	virtual ~LikelihoodDistribution();
	virtual double Evaluate(ParameterSet Theta);
	//protected:
	virtual vector<double> GetData();
	vector<double> DATA;
	bool UseEmulator;
	ofstream emulator_test;
	EmulatorHandler * emulator;
};

class PriorDistribution:public Distribution {
public:
	PriorDistribution();
	virtual ~PriorDistribution();
	virtual double Evaluate(ParameterSet Theta);
};

/** ---------------------------------------- */

class PriorDistribution_RHIC:public PriorDistribution{
public:
	PriorDistribution_RHIC(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class PriorDistribution_Cosmo:public PriorDistribution {
public:
	PriorDistribution_Cosmo(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class LikelihoodDistribution_RHIC:public LikelihoodDistribution{
public:
	LikelihoodDistribution_RHIC(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution_RHIC();
	double Evaluate(ParameterSet Theta);
	//private:
	vector<double> GetData();
	vector<double> GetRealData(); //This is just a temporary name, this will replace "GetData()"
	vector<double> DATA;
	bool UseEmulator;
	ofstream emulator_test;
	EmulatorHandler * emulator;
};

class LikelihoodDistribution_Cosmo:public LikelihoodDistribution {
public:
	LikelihoodDistribution_Cosmo(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution_Cosmo();
	double Evaluate(ParameterSet Theta);
	//private:
	vector<double> GetData();
	vector<double> DATA;
	vector<int> intDATA;
	bool UseEmulator;
	ofstream emulator_test;
	EmulatorHandler * emulator;
};

#endif