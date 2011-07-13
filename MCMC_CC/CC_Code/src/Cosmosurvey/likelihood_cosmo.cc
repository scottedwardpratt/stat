#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>

using namespace std;

LikelihoodDistribution_Cosmo::LikelihoodDistribution_Cosmo(MCMCConfiguration *mcmc_in){
	mcmc=mcmc_in;
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

LikelihoodDistribution_Cosmo::~LikelihoodDistribution_Cosmo(){
	delete emulator;
}

double LikelihoodDistribution_Cosmo::Evaluate(ParameterSet Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood = 0.0;
	stringstream ss;
	ifstream inputfile;
	
	if(TIMING){
		begintime = clock();
	}
	
	ss << "cosmosurvey";
	for(int i = 0; i < Theta.Names.size(); i++){
		ss << " -" << Theta.Names[i] << " " << Theta.Values[i];
	}
	ss << " -nz 10 -nl 10 -ob .0406 > output.dat" << endl;
	
	cout << ss.str() << endl;
	cout << "Waiting on cosmosurvey...";
	cout.flush();
	int result = system((ss.str()).c_str());
	cout << "Done." << endl;
	
	inputfile.open("output.dat");
	
	while(!inputfile.eof()){
		double temp;
		inputfile >> temp;
		ModelMeans.push_back(temp);
	}
	
	
	for(int i = 1; i < ModelMeans.size(); i++){
		likelihood += log(gsl_ran_poisson_pdf(static_cast<unsigned int>(ModelMeans[i] + 0.5), DATA[i]));
	}
	printf("LogLikelihood=%g\b",likelihood);
	
	if(!(mcmc->LOGLIKE)){
		likelihood = exp(likelihood);
	}
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
	}
	
	return likelihood;
}

vector<double> LikelihoodDistribution_Cosmo::GetData(){
	vector<double> datameans;
	stringstream ss;
	parameterMap actualparmap;
	ifstream inputfile;
	
	string actual_filename = mcmc->parameterfile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	ParameterSet ActualParams;
	ActualParams.Initialize(temp_names, temp_values);
	
	ss << "cosmosurvey";
	for(int i = 0; i < temp_names.size(); i++){
		ss << " -" << temp_names[i] << " " << temp_values[i];
	}
	ss << " -nz 10 -nl 10 -ob .0406 > data_output.dat" << endl;
	
	cout << ss.str() << endl;
	cout << "Waiting on cosmosurvey...";
	cout.flush();
	int result = system((ss.str()).c_str());
	cout << "Done." << endl;
	
	inputfile.open("data_output.dat");
	
	int counter = 0;
	while(!inputfile.eof()){
		counter++;
		double temp;
		inputfile >> temp;
		cout << "Data point " << counter << ": " << temp << endl;
		datameans.push_back(temp);
	}
	return datameans;
}
#endif