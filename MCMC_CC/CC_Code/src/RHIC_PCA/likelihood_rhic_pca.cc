#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>

using namespace std;

LikelihoodDistribution_RHIC_PCA::LikelihoodDistribution_RHIC_PCA(MCMCConfiguration *mcmc_in){
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
		quad = new QuadHandler(parmap, mcmc_in);
	}
	else{
		exit(1);
	}

	DATA = GetRealData();	
}

LikelihoodDistribution_RHIC_PCA::~LikelihoodDistribution_RHIC_PCA(){
	delete quad;
}

double LikelihoodDistribution_RHIC_PCA::Evaluate(ParameterSet Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood;
	
	if(TIMING){
		begintime = clock();
	}
	
	if(UseEmulator){
		quad->QueryQuad(Theta, ModelMeans, ModelErrors); //fills vectors with quad output
	}
	else{
		//determine another way to fill the vectors
	}
	
	//Initialize GSL containers
	int N = ModelErrors.size();
	//cout << "ModelErrors.size() = " << N << endl;
	gsl_matrix * sigma = gsl_matrix_calloc(N,N);
	//gsl_matrix * sigma_data = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		gsl_matrix_set(sigma,i,i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, DATA[i]);
		//cout << "ModelErrors " << ModelErrors[i] << " ModelMeans " << ModelMeans[i] << " Data " << DATA[i] << endl;
		//cout << "i: " << i << " Data: " << DATA[i] <<  " Mean: " << ModelMeans[i] << " Error: " << ModelErrors[i] << endl;
	}
	
	likelihood = Log_MVNormal(*model, *mu, *sigma);
	//likelihood = Gaussian(*model, *mu, *sigma);
	//likelihood = Gaussian(*model, *mu, *sigma, *sigma_data); //This is the integrated likelihood
	
	if(!(mcmc->LOGLIKE)){ //If you don't want the loglikelihood, and we've used Log_MVN, we have to exponentiate.
		likelihood = exp(likelihood);
	}
	/*if(mcmc->LOGLIKE){ //If you do want the loglikelihood, and we've used Gaussian, we have to take the log.
		likelihood = log(likelihood);
	}*/
	
	if(VERBOSE){
		double sum = 0.0;
		
		for(int i = 0; i< N; i++){
			sum += (gsl_vector_get(model, i) - gsl_vector_get(mu, i));
		}
		sum = sum/(double)N;
		cout << "Average difference between outputs:" << sum << endl;
	}
	
	//deallocate GSL containers.
	gsl_vector_free(model);
	gsl_vector_free(mu);
	gsl_matrix_free(sigma);
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
	}
	
	//cout << "PCA 0: " << ModelMeans[0] << endl;
	
	return likelihood;
}

vector<double> LikelihoodDistribution_RHIC_PCA::GetRealData(){
	vector<double> datameans;
	float myints[] = {403.723,467.153,743.79,1042.58,5.28,4.81,5.47,180.682,460.801,739.43,1025.55,0.0891618,0.0271948,0.0498033,4.27,3.99,4.53};
	datameans.assign (myints,myints+17);
	//float myints[] = {0.658581,0.545024,1.4456,-2.03515,1.22693,0.731763};
	//datameans.assign (myints,myints+6);
	return datameans;
}

#endif