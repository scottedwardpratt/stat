#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

LikelihoodDistribution::LikelihoodDistribution(MCMC *mcmc_in):Distribution(mcmc_in){
	emulator = new EmulatorHandler(mcmc);
	UseEmulator = parameter::getB(mcmc->parmap, "LIKELIHOOD_USE_EMULATOR", false);
	
	DATA = GetData();
}

LikelihoodDistribution::~LikelihoodDistribution(){
	delete emulator;
}

double LikelihoodDistribution::Evaluate(ParameterSet Theta){
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
	//gsl_vector * diff2 = gsl_vector_alloc(ModelErrors.size());
	gsl_vector * temp = gsl_vector_alloc(ModelErrors.size());
	
	//Read in appropriate elements
	for(int i = 0; i<ModelErrors.size(); i++){
		gsl_matrix_set(sigma, i,i,ModelErrors[i]);
		gsl_vector_set(diff, i, ModelMeans[i]-DATA[i]);
	}
	
	//invert matrix using cholesky decomposition
	foobar = gsl_linalg_cholesky_decomp(sigma);
	foobar = gsl_linalg_cholesky_invert(sigma);
	
	//multiply matrix and left vector together using CBLAS routines
	gsl_blas_dgemv(CblasNoTrans,1.0, sigma, diff, 0.0, temp);
	
	likelihood = (-1/2)*gsl_vector_mul(diff, temp);
	
	if(!(mcmc->LOGLIKE)){
		likelihood = exp(likelihood);
	}
	
	//deallocate GSL containers.
	gsl_vector_free(diff);
	gsl_vector_free(temp);
	gsl_matrix_free(sigma);
	
	return likelihood;
}

//Fix this later.
vector<double> LikelihoodDistribution::GetData(){
	vector<double> data (50, 0);
	
	return data;
}
#endif