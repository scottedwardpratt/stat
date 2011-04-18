#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

LikelihoodDistribution::LikelihoodDistribution(MCMC *mcmc_in):Distribution(mcmc_in){
	cout << "Likelihood:Start" << endl;
	SepMap = parameter::getB(mcmc->parmap, "LIKELIHOOD_PARAMETER_MAP", false);
	
	if(SepMap){
		// cout << "Seperate parameter map for likelihood." << endl;
		string parmapfile = mcmc->dir_name + "/mcmc/parameters/likelihood.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
		//parameter::ReadParsFromFile(parmap, parameter_file_name);
	}else{
		// cout << "Using MCMC parameter map for likelihood." << endl;
		parmap = &(mcmc->parmap);
	}
	
	UseEmulator = parameter::getB(*parmap, "USE_EMULATOR", false);
	
	if(UseEmulator){
		// cout << "About to create new emulator handler." << endl;
		emulator = new EmulatorHandler(parmap, mcmc_in);
	}
	else{
		// cout << "Unable to use necessary emulator." << endl;
		exit(1);
	}

	DATA = GetData();
	cout << "Likelihood: Done." << endl;
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
	cout << "Reading in data." << endl;
	vector<double> datameans;
	vector<double> dataerror;
	//vector<string> stringtemp;
	//vector<double> doubtemp;
	
	parameterMap actualparmap;
	
	string actual_filename = mcmc->dir_name+"/mcmc/parameters/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	//cout << "Reading in actual parameters from " << actual_filename << endl;
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	//cout << "Parameter names read in." << endl;
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	//cout << "Parameter values read in." << endl;
	
	ParameterSet ActualParams(mcmc->ThetaList);
	//cout << "New parameterset generated." << endl;
	ActualParams.Initialize(temp_names, temp_values);
	//cout << "Parameters initialized." << endl;
	emulator->QueryEmulator(ActualParams, datameans, dataerror);
	cout << "Done reading in data." << endl;
	return datameans;
}
#endif