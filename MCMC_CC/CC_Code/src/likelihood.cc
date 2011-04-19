#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>

using namespace std;

LikelihoodDistribution::LikelihoodDistribution(MCMC *mcmc_in):Distribution(mcmc_in){
	SepMap = parameter::getB(mcmc->parmap, "LIKELIHOOD_PARAMETER_MAP", false);
	
	if(SepMap){
		string parmapfile = mcmc->dir_name + "/mcmc/parameters/likelihood.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
	}else{
		parmap = &(mcmc->parmap);
	}
	
	UseEmulator = parameter::getB(*parmap, "USE_EMULATOR", false);
	TIMING = parameter::getB(*parmap, "TIMING", false) || parameter::getB(*parmap, "TIME_LIKELIHOOD", false);
	
	if(UseEmulator){
		emulator = new EmulatorHandler(parmap, mcmc_in);
	}
	else{
		exit(1);
	}

	DATA = GetData();
}

LikelihoodDistribution::~LikelihoodDistribution(){
	delete emulator;
}

double LikelihoodDistribution::Evaluate(ParameterSet Theta){
	clock_t begintime;
	if(TIMING){
		begintime = clock();
	}
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood;
	int foobar;
	
	if(UseEmulator){
		emulator->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
	}
	else{
		//determine another way to fill the vectors
	}
	
	//Initialize GSL containers
	gsl_matrix * sigma = gsl_matrix_calloc(ModelErrors.size(), ModelErrors.size());
	gsl_vector * diff = gsl_vector_alloc(ModelErrors.size());
	gsl_vector * temp = gsl_vector_alloc(ModelErrors.size());
	// cout << "Done allocating gsl containers." << endl;
	
	//Read in appropriate elements
	for(int i = 0; i<ModelErrors.size(); i++){
		gsl_matrix_set(sigma, i,i,Theta.GetValue("SIGMA"));
		gsl_vector_set(diff, i, ModelMeans[i]-DATA[i]);
	}
	
	//invert matrix using cholesky decomposition
	foobar = gsl_linalg_cholesky_decomp(sigma);
	foobar = gsl_linalg_cholesky_invert(sigma);
	
	//multiply matrix and left vector together using CBLAS routines
	gsl_blas_dgemv(CblasNoTrans,1.0, sigma, diff, 0.0, temp);
	
	gsl_blas_ddot(diff, temp, &likelihood);
	
	likelihood = (-1.0/2.0)*likelihood;
	
	if(!(mcmc->LOGLIKE)){
		likelihood = exp(likelihood);
	}
	
	//deallocate GSL containers.
	gsl_vector_free(diff);
	gsl_vector_free(temp);
	gsl_matrix_free(sigma);
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
		// cout << CLOCKS_PER_SEC << endl;
	}
	
	return likelihood;
}

vector<double> LikelihoodDistribution::GetData(){
	vector<double> datameans;
	vector<double> dataerror;
	
	parameterMap actualparmap;
	
	string actual_filename = mcmc->dir_name+"/mcmc/parameters/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	ParameterSet ActualParams(mcmc->ThetaList);
	ActualParams.Initialize(temp_names, temp_values);
	emulator->QueryEmulator(ActualParams, datameans, dataerror);
	return datameans;
}
#endif