#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>

using namespace std;

LikelihoodDistribution_Test::LikelihoodDistribution_Test(MCMC *mcmc_in){
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
		cout << "This is a guassian test. The UseEmulator flag will be ignored." << endl;
	}
	else{
		exit(1);
	}

	DATA = GetData();

	//testing the outputs of the emulator at various points	// 
	emulator_test.open("PCA0.dat");
	emulator_test.close();
	
}

LikelihoodDistribution_Test::~LikelihoodDistribution_Test(){
}

double LikelihoodDistribution_Test::Evaluate(vector<double> Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood;
	
	if(TIMING){
		begintime = clock();
	}
	

	//Fill vectors:
	if(UseEmulator){
		// This is a gaussian test, so there shouldn't be an 'emulator' per say.
		//emulator->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
	}
	else{
		int errors[]={0,0,0,0}; //For now
		ModelErrors.assign (errors,errors+4);
		int means[4];
		for(int i=0;i<4;i++){
			means[i]=Theta[i];
		}
		ModelMeans.assign(means,means+4);
	}
	
	
	//Initialize GSL containers
	int N = ModelErrors.size();
	//gsl_matrix * sigma = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		// gsl_matrix_set(sigma, i,i,Theta.GetValue("SIGMA"));
		//gsl_matrix_set(sigma,i,i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, DATA[i]);
	}
	
	//likelihood = Log_MVNormal(*model, *mu, *sigma);
	
	for(int i=0;i<4;i++){
		likelihood+=(ModelMeans[i]-DATA[i])*(ModelMeans[i]-DATA[i])/2; //This is a gaussian, but it would be "log(exp(stuff))", so I just left it "stuff"
	}

	if(!(mcmc->LOGLIKE)){
		likelihood = exp(likelihood);
	}
	
	if(VERBOSE){
		/*double sum = 0.0;
		
		for(int i = 0; i< N; i++){
			sum += (gsl_vector_get(model, i) - gsl_vector_get(mu, i));
		}
		sum = sum/(double)N;
		cout << "Average difference between outputs:" << sum << endl;*/
	}
	
	//deallocate GSL containers.
	gsl_vector_free(model);
	gsl_vector_free(mu);
	//gsl_matrix_free(sigma);
	
	if(TIMING){
		cout << "Likelihood evaluation took " << (clock()-begintime)*1000/CLOCKS_PER_SEC << " ms." << endl;
	}
	
	cout << "PCA 0: " << ModelMeans[0] << endl;
	
	emulator_test.open("PCA0.dat", ios_base::app);
	emulator_test << ModelMeans[0] << endl;
	emulator_test.close();
	// emulator_test << ModelMeans[0] << endl;
	
	return likelihood;
}

vector<double> LikelihoodDistribution_Test::GetData(){
	//Four separate gaussians with means of 1, 2, 3, and 4 and standard deviations of 2. The means will be the input, 
	//and three values on the gaussian will be the output.

	vector<double> datameans;
	int mymeans[]={1,2,3,4};
	datameans.assign (mymeans,mymeans+4);
	//vector<double> dataerror;
	/*
	parameterMap actualparmap;
	
	string actual_filename = mcmc->parameterfile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	vector<double> ActualParams;
	ActualParams.Initialize(temp_names, temp_values);
	
	emulator->QueryEmulator(ActualParams, datameans, dataerror);*/

	return datameans;
}
#endif