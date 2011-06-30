#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>

using namespace std;

LikelihoodDistribution::LikelihoodDistribution(MCMCConfiguration *mcmc_in):Distribution(mcmc_in){
	SepMap = parameter::getB(mcmc->parmap, "LIKELIHOOD_PARAMETER_MAP", false);
	
	if(SepMap){
		string parmapfile = mcmc->parameterfile + "/likelihood.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
	}else{
		parmap = &(mcmc->parmap);
	}
	
	// cout << "Like param map made." << endl;
	
	UseEmulator = parameter::getB(*parmap, "USE_EMULATOR", false);
	TIMING = parameter::getB(*parmap, "TIMING", false) || parameter::getB(*parmap, "TIME_LIKELIHOOD", false);
	VERBOSE = parameter::getB(*parmap, "VERBOSE", false) || parameter::getB(*parmap, "VERBOSE_LIKELIHOOD", false);
	
	// cout << "params declared." << endl;
	if(UseEmulator){
		cout << "Turn off USE_EMULATOR" << endl;
		exit(-1);
		emulator = new EmulatorHandler(parmap, mcmc_in);
	}

	DATA = GetData();
	
}

LikelihoodDistribution::~LikelihoodDistribution(){
	delete emulator;
}

double LikelihoodDistribution::Evaluate(ParameterSet Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood = 0.0;
	stringstream ss;
	ifstream inputfile;
	
	if(TIMING){
		begintime = clock();
	}
	
	for(int i = 0; i < Theta.Values.size(); i++){
		cout << "Comparing " << Theta.Values[i] << " against " << DATA[i] << endl;
		double templike = log(gsl_ran_gaussian_pdf(Theta.Values[i]-DATA[i], .1));
		cout << "Likelihood: " << templike << endl;
		likelihood += templike;
	}
	
	if(!(mcmc->LOGLIKE)){
		cout << "Exponentiating." << endl;
		likelihood = exp(likelihood);
	}
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
	}
	
	cout << "Likelihood: " << likelihood << endl;
	
	return likelihood;
}

vector<double> LikelihoodDistribution::GetData(){
	vector<double> datameans;
	stringstream ss;
	parameterMap actualparmap;
	ifstream inputfile;
	
	string actual_filename = mcmc->parameterfile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	return temp_values;
}
#endif