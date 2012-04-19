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

	DATA = GetRealData();
	ERROR = GetRealError();

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
	gsl_matrix * sigma_data = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		//gsl_matrix_set(sigma, i,i,Theta.GetValue("SIGMA"));
		//ModelErrors[i]=0.1; //Suppressing the error of the emulator
		if (ModelErrors[i] < 1) { ModelErrors[i]=1; }
		//Bad things can happen in the likelihood calculations if the errors get too large or small
		//ModelErrors[i]=ModelErrors[i]*2;
//		if (ModelErrors[i] < 0.5*fabs(DATA[i]-ModelMeans[i])) { ModelErrors[i]=0.5*fabs(DATA[i]-ModelMeans[i]); }
//		if (ModelErrors[i] > 5*fabs(DATA[i]-ModelMeans[i])) { ModelErrors[i]=5*fabs(DATA[i]-ModelMeans[i]); }
		gsl_matrix_set(sigma,i,i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, DATA[i]);
		gsl_matrix_set(sigma_data,i,i,ERROR[i]);
		//cout << "i: " << i << " Data: " << DATA[i] <<  " Mean: " << ModelMeans[i] << " Error: " << ModelErrors[i] << endl;
	}
	
	//likelihood = Log_MVNormal(*model, *mu, *sigma);
	//likelihood = Gaussian(*model, *mu, *sigma);
	likelihood = Gaussian(*model, *mu, *sigma, *sigma_data); //This is the integrated likelihood
	
	/*if(!(mcmc->LOGLIKE)){ //If you don't want the loglikelihood, and we've used Log_MVN, we have to exponentiate.
		likelihood = exp(likelihood);
	}*/
	if(mcmc->LOGLIKE){ //If you do want the loglikelihood, and we've used Gaussian, we have to take the log.
		likelihood = log(likelihood);
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

vector<double> LikelihoodDistribution_RHIC::GetFakeData(){
	//This makes some fake results by querying the emulator with the parameters from actual.param
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
	vector<double> datameans;
	string EmulatorObservables=emulator->Observables;
	//cout << "Observables: " << EmulatorObservables << endl;
	string observables_filename = mcmc->dir_name + "/" + EmulatorObservables + ".datnames";
	cout << "Observables being read from: " << observables_filename << endl;
	string data_filename = mcmc->dir_name + "/exp_data/results.dat";
	//string data_filename = mcmc->dir_name + "/analysis/exp_data/" + "cent20to30" + "/results.dat";
	cout << "Results being read from: " << data_filename << endl;
	
	fstream data;
	string type, param_name;
	int count = 0;
	
	parameter::ReadParsFromFile(observablesparmap, observables_filename.c_str());
	vector<string> PNames;
	PNames = parameter::getVS(observablesparmap,"NAMES","blah blah");
	/*ifstream emulated_observables (observables_filename.c_str());

	while(!emulated_observables.eof()){
		string line;
		getline(emulated_observables, line, '\n');
		if(line.compare(0,1,"#") != 0 && !line.empty() ){ //if the line is not a comment
			PNames.push_back(line);
		}
	}
	emulated_observables.close();
*/
	int numparams = PNames.size();
	cout << "There are " << numparams  << " parameters used in the emulator." << endl;
	Datamean=new double[numparams];
	Dataerror=new double[numparams];
	vector<double> temp (numparams, .01);
	double dump;
	
	data.open(data_filename.c_str(), fstream::in);
	if(data){
		while(data >> type){
			if(strcmp(type.c_str(), "double") == 0){
				data >> param_name;
				int index = FindParam(param_name, PNames);
				if(index != -1){ //returns -1 if not found
					data >> Datamean[index]; //Mean
					data >> Dataerror[index]; //Error
					count++;
					cout << param_name << " index: " << index << " " << Datamean[index] << endl;
				}else{
					cout << "Not using observable: " << param_name << endl;
					data >> dump; //we aren't using the observable, so we need to get the data out of the stream
					data >> dump;
					//exit(1);
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(data, temp, '\n');
				}
				else{
					cout << "Unrecognized variable type " << type << endl;
					exit(1);
				}
			}
		}
		data.close();
	}else{
		cout << "Warning: Unable to open data file in model directory." << endl;
		 exit(1);
	}

	if(count!=numparams){
		cout << "Not all emulated observables found!" << endl;
		exit(1);
	}

	datameans.assign(Datamean,Datamean+numparams);

	if(VERBOSE){
		for(int i = 0; i<numparams; i++){
			cout << "Data array: " << endl;
			cout << "i: " << i << " Data: " << datameans[i] << endl;
		}
	}
	return datameans;
}

vector<double> LikelihoodDistribution_RHIC::GetRealError(){
	vector<double> dataerrors;
	string EmulatorObservables=emulator->Observables;
	string observables_filename = mcmc->dir_name + "/" + EmulatorObservables + ".datnames";
	string data_filename = mcmc->dir_name + "/exp_data/results.dat";
	
	fstream data;
	string type, param_name;
	int count = 0;
	
	parameter::ReadParsFromFile(observablesparmap, observables_filename.c_str());
	vector<string> PNames;
	PNames = parameter::getVS(observablesparmap,"NAMES","blah blah");

	int numparams = PNames.size();
	Datamean=new double[numparams];
	Dataerror=new double[numparams];
	vector<double> temp (numparams, .01);
	double dump;
	
	data.open(data_filename.c_str(), fstream::in);
	if(data){
		while(data >> type){
			if(strcmp(type.c_str(), "double") == 0){
				data >> param_name;
				int index = FindParam(param_name, PNames);
				if(index != -1){ //returns -1 if not found
					data >> Datamean[index]; //Mean
					data >> Dataerror[index]; //Error
					count++;
					cout << param_name << " index: " << index << " " << Datamean[index] << endl;
				}else{
					data >> dump; //we aren't using the observable, so we need to get the data out of the stream
					data >> dump;
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(data, temp, '\n');
				}
				else{
					cout << "Unrecognized variable type " << type << endl;
					exit(1);
				}
			}
		}
		data.close();
	}else{
		cout << "Warning: Unable to open data file in model directory." << endl;
		 exit(1);
	}

	if(count!=numparams){
		cout << "Not all emulated observables found!" << endl;
		exit(1);
	}

	dataerrors.assign(Dataerror,Dataerror+numparams);

	if(VERBOSE){
		for(int i = 0; i<numparams; i++){
			cout << "Error array: " << endl;
			cout << "i: " << i << " Error: " << dataerrors[i] << endl;
		}
	}
	return dataerrors;
}

int LikelihoodDistribution_RHIC::FindParam(string name, vector<string> PNames){ 
	int out = -1;
	int i = 0;
	bool Found = false;
	
	while(i < PNames.size()){
		// cout << "FindParam: Comparing " << name << " to " << PNames[i] << endl;
		if(strcmp(PNames[i].c_str(), name.c_str()) == 0){
			if(!Found){
				out = i;
				Found = true;
			}else{ //A matching parameter has already been found, multiple parameters with the same name.
				cout << PNames[out] << endl;
				cout << PNames[i] << endl;
				cout << "In ProposalName::FindParam; Duplicate parameter names found. Please change parameter names." << endl;
				exit(1);
			}
		}
		i++;
	}
	return out;
}
#endif