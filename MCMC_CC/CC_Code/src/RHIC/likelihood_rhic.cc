#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>

using namespace std;

LikelihoodDistribution_RHIC::LikelihoodDistribution_RHIC(MCMCConfiguration *mcmc_in){
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
		emulator = new EmulatorHandler(parmap, mcmc_in);
	}
	else{
		exit(1);
	}

	DATA = GetData();

	//testing the outputs of the emulator at various points	// 
	emulator_test.open("PCA0.dat");
	emulator_test.close();
	
}

LikelihoodDistribution_RHIC::~LikelihoodDistribution_RHIC(){
	delete emulator;
}

double LikelihoodDistribution_RHIC::Evaluate(ParameterSet Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood;
	
	if(TIMING){
		begintime = clock();
	}
	
	if(UseEmulator){
		emulator->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
	}
	else{
		//determine another way to fill the vectors
	}
	
	//Initialize GSL containers
	int N = ModelErrors.size();
	gsl_matrix * sigma = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		// gsl_matrix_set(sigma, i,i,Theta.GetValue("SIGMA"));
		gsl_matrix_set(sigma,i,i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, DATA[i]);
	}
	
	likelihood = Log_MVNormal(*model, *mu, *sigma);
	
	if(!(mcmc->LOGLIKE)){
		likelihood = exp(likelihood);
	}
	
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
	
	cout << "PCA 0: " << ModelMeans[0] << endl;
	
	emulator_test.open("PCA0.dat", ios_base::app);
	emulator_test << ModelMeans[0] << endl;
	emulator_test.close();
	// emulator_test << ModelMeans[0] << endl;
	
	return likelihood;
}

vector<double> LikelihoodDistribution_RHIC::GetData(){
	vector<double> datameans;
	vector<double> dataerror;
	
	parameterMap actualparmap;
	
	string actual_filename = mcmc->parameterfile + "/actual.param";
	parameter::ReadParsFromFile(actualparmap, actual_filename);
	
	vector<string> temp_names = parameter::getVS(actualparmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(actualparmap, "VALUES", "");
	
	ParameterSet ActualParams;
	ActualParams.Initialize(temp_names, temp_values);
	
	emulator->QueryEmulator(ActualParams, datameans, dataerror);

	return datameans;
}

vector<double> LikelihoodDistribution_RHIC::GetRealData(){
/*	vector<double> datameans;
	vector<double> dataerror;
	
	parameterMap observablesparmap; //This is the observables in the emulator ouput

	string EmulatorObservables=emulator->Observables; //The name of the observables we are working with
	
	string observables_filename = mcmc->parameterfile + "/" + EmulatorObservables + ".param"; //I would like to just call 
	//on the observables list in "data-prep" but I think the format is wrong. If I am going to do it the way it is right
	//now, I will need to add a line to setup.sh to generate the observabes list.
	parameter::ReadParsFromFile(observablesparmap, observables_filename);
	
	//I will need to open up the real data, then I want to run over it with a for loop. Each pass will take one of the
	//emulate observables and find a match, then use that to fill the datameans and dataerrors vectors.
	//Becuase of the for loop method, I don't know if it is best to actually store the emulated observables in a
	//parameter array. Perhaps I should store the data in a parameter array and just do a readline loop for the
	//emulator observables.

	// *********************/
	// Here we go:
	vector<double> datameans;
	parameterMap observablesparmap;
	string EmulatorObservables=emulator->Observables;
	string observables_filename = "data-prep/" + EmulatorObservables + ".param";
	string data_filename = "pathname for real data";
	parameter::ReadParsFromFile(observablesparmap, data_filename);

	ifstream emulated_observables (observables_filename.c_str());

	while(!emulated_observables.eof()){
		string line;
		double mean, error;
		getline(emulated_observables, line, '\n');
		//emulated_observables.getline(line,100);
		if(line.compare(0,1,"#") != 0 && !line.empty() ){ //if the line is not a comment
			// How this readline will work will depend on how the data file is formated
			mean = parameter::getD(observablesparmap,line,-999);
			datameans.push_back(mean);
		}
	}
	
	return datameans;
}
#endif