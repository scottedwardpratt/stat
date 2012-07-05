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
class QuadHandler;
class MCMCConfiguration;
class MCMCRun;
class ParameterSet;
class emulator;

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
	double IntegratedNormal(double x, double mu, double sigma, double data_sigma);
	double Gaussian(double x, double mu, double sigma);
	double Gaussian(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double Gaussian(gsl_vector x, gsl_vector mu, gsl_matrix sigma, gsl_matrix data_sigma);
	double Log_MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma);
	double LogNormal(double x, double mu, double sigma);
};

/** ---------------------------------------- */

class ProposalDistribution:public Distribution {
public:
	ProposalDistribution(MCMCConfiguration *mcmc_in);
	ParameterSet Iterate(ParameterSet current, float scale);
	virtual double Evaluate(ParameterSet Theta1, ParameterSet Theta2, float scale);
	//protected:
	bool SymmetricProposal;
	bool Rescaled_Method;
	double SCALE;
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
	//emulator * My_emu;
	//EmulatorHandler * emulator;
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

class PriorDistribution_RHIC_PCA:public PriorDistribution{
public:
	PriorDistribution_RHIC_PCA(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class PriorDistribution_Cosmo:public PriorDistribution {
public:
	PriorDistribution_Cosmo(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class PriorDistribution_Test:public PriorDistribution {
public:
	PriorDistribution_Test(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
};

class LikelihoodDistribution_RHIC:public LikelihoodDistribution{
public:
	LikelihoodDistribution_RHIC(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution_RHIC();
	double Evaluate(ParameterSet Theta);
	double *Datamean;
	double *Dataerror;
	//private:
	vector<double> GetFakeData();
	vector<double> GetRealData(); 
	vector<double> GetRealError(); 
	vector<double> DATA;
	vector<double> ERROR;
	bool UseEmulator;
	bool FAKE_DATA;
	emulator * My_emu;
	//parameterMap * parmap;
	ofstream emulator_test;
	int FindParam(string param_name, vector<string> PNames);
	parameterMap observablesparmap;
};

class LikelihoodDistribution_RHIC_PCA:public LikelihoodDistribution{
public:
	LikelihoodDistribution_RHIC_PCA(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution_RHIC_PCA();
	double Evaluate(ParameterSet Theta);
	double *Datamean;
	double *Dataerror;
	//private:
	vector<double> GetRealData(); 
	vector<double> DATA;
	vector<double> ERROR;
	bool UseEmulator;
	ofstream emulator_test;
	QuadHandler * quad;
	parameterMap observablesparmap;
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
	//emulator * My_emu;
};

class LikelihoodDistribution_Test:public LikelihoodDistribution {
public:
	LikelihoodDistribution_Test(MCMCConfiguration *mcmc_in);
	~LikelihoodDistribution_Test();
	double Evaluate(ParameterSet Theta);
	//private:
	vector<double> GetData();
	vector<double> DATA;
	bool UseEmulator;
	ofstream emulator_test;
};

#endif