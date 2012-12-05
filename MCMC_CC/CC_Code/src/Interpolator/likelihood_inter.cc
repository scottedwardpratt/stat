#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>
#include "rhicstat.h"

using namespace std;

LikelihoodDistribution_Interpolator::LikelihoodDistribution_Interpolator(MCMC *mcmc_in){
	mcmc=mcmc_in;
	SepMap = parameter::getB(mcmc->parmap, "LIKELIHOOD_PARAMETER_MAP", false);
	if(SepMap){
		string parmapfile = mcmc->parameterfile + "/likelihood.param";
		parmap = new parameterMap;
		//parameterMap *parmap;
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
		cout << "Emulator is being loaded from: " << mcmc->dir_name << endl;
		My_emu = new CRHICStat(mcmc->dir_name);
		//CRHICStat *My_emu = new CRHICStat(mcmc->dir_name);
	}
	else{
		cout << "The UseEmulator flag is set to false (or not set). We can't do anything without an emulator" << endl;
		exit(1);
	}

	cout << "We are using Scott's emulator, so all parameters are being varied in a scaled space from -root(3) to +root(3)" << endl;
	for(int i=0; i < mcmc->ParamNames.size(); i++){
		double r3 = sqrt(3);
		mcmc->Min_Ranges[i]=-r3;
		mcmc->Max_Ranges[i]=r3;
	}
}

LikelihoodDistribution_Interpolator::~LikelihoodDistribution_Interpolator(){
	//delete My_emu;
}

double LikelihoodDistribution_Interpolator::Evaluate(vector<double> Theta){
	clock_t begintime;
	double likelihood;
	double *x = &Theta[0];

	if(VERBOSE){
		cout << "Theta: ";
		for(int i = 0; i < Theta.size(); i++){
			cout << Theta[i] << " ";
		}
		cout << endl;
		/*for(int i = 0; i < Theta.size(); i++){
			cout << x[i] << " ";
		}
		cout << endl;*/
	}
	
	if(TIMING){
		begintime = clock();
	}

	if(UseEmulator){
		//GetLL takes an array of doubles, this should work, but if for some reason the values in Theta are not stored contiguosly, this may barf
		//likelihood = My_emu->GetLL(x);
		//cout << "Using 'x': " << likelihood << endl;
		likelihood = My_emu->GetLL(&Theta[0]);
		//cout << "Using '&Theta[0]': " << likelihood << endl;
	}
	else{
		//determine another way to fill the vectors
	}
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
	}

	return likelihood;
}
#endif
