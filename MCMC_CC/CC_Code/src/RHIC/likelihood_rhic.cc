#ifndef __LIKELIHOOD_CC__
#define __LIKELIHOOD_CC__

#include "mcmc.h"
#include <time.h>
#include "EmuPlusPlus.h"

using namespace std;

LikelihoodDistribution_RHIC::LikelihoodDistribution_RHIC(MCMC *mcmc_in){
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
	FAKE_DATA = parameter::getB(*parmap, "FAKE_DATA", false);

	// cout << "params declared." << endl;
	if(UseEmulator){
		cout << "Emulator is being loaded from: " << mcmc->dir_name + "/Emulator.statefile" << endl;
		My_emu = new emulator(mcmc->dir_name + "/Emulator.statefile");
		//emulator My_emu(mcmc->dir_name + "/Emulator.statefile");
		cout << "Emulator loaded. Test (number of params): _" << My_emu->number_params << "_" << endl;
	}
	else{
		cout << "The UseEmulator flag is set to false (or not set). We can't do anything without an emulator" << endl;
		exit(1);
	}

	if(FAKE_DATA){
		DATA = GetFakeData();
	}else{
		DATA = GetRealData();
	}
	//ERROR = GetRealError();

	//testing the outputs of the emulator at various points	// 
	//emulator_test.open("PCA0.dat");
	//emulator_test.close();
	
}

LikelihoodDistribution_RHIC::~LikelihoodDistribution_RHIC(){
	//delete My_emu;
}

double LikelihoodDistribution_RHIC::Evaluate(vector<double> Theta){
	clock_t begintime;
	vector<double> ModelMeans;
	vector<double> ModelErrors;
	double likelihood;
	
	if(TIMING){
		begintime = clock();
	}

	if(UseEmulator){
		My_emu->QueryEmulator(Theta, ModelMeans, ModelErrors); //fills vectors with emulator output
	}
	else{
		//determine another way to fill the vectors
	}

	//Initialize GSL containers
	int N = ModelErrors.size();
	gsl_matrix * sigma = gsl_matrix_calloc(N,N);
	//gsl_matrix * sigma_data = gsl_matrix_calloc(N,N);
	gsl_vector * model = gsl_vector_alloc(N);
	gsl_vector * mu = gsl_vector_alloc(N);
	// cout << "Done allocating gsl containers." << endl;
	
	if(VERBOSE){
		cout << "Theta: ";
		for(int i = 0; i < Theta.size(); i++){
			cout << Theta[i] << " ";
		}
		cout << endl << "Observable, Model means, Model errors, Data" << endl;
		for(int i = 0; i<N; i++){
			cout << mcmc->ObservablesNames[i] << " " << ModelMeans[i] << " " << ModelErrors[i] << " " << DATA[i] << endl;
		}
	}
	//Read in appropriate elements
	for(int i = 0; i<N; i++){
		//cout << " Data: " << DATA[i] << " Emu: " << ModelMeans[i] << " +/-: " << ModelErrors[i] << endl;
		if(mcmc->SUPPRESS_ERRORS){
			//ModelErrors[i]=ModelMeans[i]*0.1; //What a reasonable error is depends on the observable
			ModelErrors[i]=1;
		}
		gsl_matrix_set(sigma, i, i,ModelErrors[i]);
		gsl_vector_set(model, i,ModelMeans[i]);
		gsl_vector_set(mu, i, DATA[i]);
		//cout << "i: " << i << " Data: " << DATA[i] <<  " Mean: " << ModelMeans[i] << " Error: " << ModelErrors[i] << endl;
	}
	
	likelihood = Log_MVNormal(*model, *mu, *sigma);
	//likelihood = Gaussian(*model, *mu, *sigma);
	//likelihood = Gaussian(*model, *mu, *sigma, *sigma_data); //This is the integrated likelihood
	
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
	
	//emulator_test.open("PCA0.dat", ios_base::app);
	//emulator_test << ModelMeans[0] << endl;
	//emulator_test.close();
	//emulator_test << ModelMeans[0] << endl;

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
	
	vector<double> ActualParams;
	ActualParams = temp_values;
	
	My_emu->QueryEmulator(ActualParams, datameans, dataerror);

	cout << "We are using FAKE DATA!!!!!!!!" << endl;
	cout << "The parameter values read in from actual.param are:" << endl;
	for(int i = 0; i < ActualParams.size(); i++){
		cout << ActualParams[i] << " ";
	}
	cout << endl << "Thses have given us parameter values of:" << endl;
	for(int i = 0; i<datameans.size(); i++){
		cout << mcmc->ObservablesNames[i] << " " << datameans[i] << endl;
	}

	return datameans;
}

vector<double> LikelihoodDistribution_RHIC::GetRealData(){
	vector<double> datameans;
	
	string data_filename = mcmc->dir_name + "/exp_data/results.dat";
	fstream data;
	string type, obsv_name;
	int count=0;
	vector<string> PNames;
	PNames = mcmc->ObservablesNames;

	int numparams = PNames.size();
	cout << "There are " << numparams  << " observables used in the emulator." << endl;
	Datamean=new double[numparams];
	Dataerror=new double[numparams];
	vector<double> temp (numparams, .01);
	double dump;
	
	data.open(data_filename.c_str(), fstream::in);
	if(data){
		while(data >> type){
			if(strcmp(type.c_str(), "double") == 0){
				data >> obsv_name;
				int index = FindParam(obsv_name, PNames);
				if(index != -1){ //returns -1 if not found
					data >> Datamean[index]; //Mean
					data >> Dataerror[index]; //Error
					count++;
					cout << obsv_name << " index: " << index << " " << Datamean[index] << endl;
				}else{
					cout << "Not using observable: " << obsv_name << endl;
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
		cout << "Not all emulated observables found! count=" << count << " numparams= " << numparams << endl;
		exit(1);
	}

	datameans.assign(Datamean,Datamean+numparams);

	if(VERBOSE){
		cout << "Data array: " << endl;
		for(int i = 0; i<numparams; i++){
			cout << mcmc->ObservablesNames[i] << " " << datameans[i] << endl;
		}
	}
	return datameans;
}

vector<double> LikelihoodDistribution_RHIC::GetRealError(){
	vector<double> dataerrors;
	string data_filename = mcmc->dir_name + "/exp_data/results.dat";
	
	fstream data;
	string type, obsv_name;
	int count = 0;
	
	vector<string> PNames;
	PNames = mcmc->ObservablesNames;

	int numparams = PNames.size();
	Datamean=new double[numparams];
	Dataerror=new double[numparams];
	vector<double> temp (numparams, .01);
	double dump;
	
	data.open(data_filename.c_str(), fstream::in);
	if(data){
		while(data >> type){
			if(strcmp(type.c_str(), "double") == 0){
				data >> obsv_name;
				int index = FindParam(obsv_name, PNames);
				if(index != -1){ //returns -1 if not found
					data >> Datamean[index]; //Mean
					data >> Dataerror[index]; //Error
					count++;
					cout << obsv_name << " index: " << index << " " << Datamean[index] << endl;
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
