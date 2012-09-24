#ifndef __RHIC_H__
#define __RHIC_H__

#include "distribution.h"

class PriorDistribution_RHIC:public PriorDistribution{
public:
	PriorDistribution_RHIC(MCMCConfiguration *mcmc_in);
	double Evaluate(ParameterSet Theta);
	string PRIOR;
	bool SCALED;
	vector<double> GAUSSIAN_MEANS;
	vector<double> GAUSSIAN_STDVS;
	vector<double> STEP_MEANS;
	vector<string> STEP_SIDE;
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

#endif