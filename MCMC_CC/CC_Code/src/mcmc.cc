#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMCConfiguration::MCMCConfiguration(string run_file){
	configname = "default";
	dir_name = run_file;
	parameterfile = run_file+"/mcmc/parameters/" + configname;
	cout << "In config: " << parameterfile << endl;
	parameter_file_name = parameterfile + "/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	
	EmulatorParams = parameter::getS(parmap, "EMULATOR_PARAMETERS", "");
	ParamNames = parameter::getVS(parmap, "PARAMETER_NAMES", "");
	vector<string> temp_logparam = parameter::getVS(parmap, "LOG_PARAMETERS", "");
	
	
	for(int i =0; i<temp_logparam.size(); i++){
		if(strcmp(temp_logparam[i].c_str(), "true") == 0 || strcmp(temp_logparam[i].c_str(), "True") == 0){
			LogParam.push_back(true);
		}else if(strcmp(temp_logparam[i].c_str(), "false") == 0 || strcmp(temp_logparam[i].c_str(), "False") == 0){
			LogParam.push_back(false);
		}else{
			cout << "Unrecognized LogParam value " << temp_logparam[i] << endl;
			exit(1);
		}
	}
	
	randnum = new CRandom(1234);
	Likelihood = new LikelihoodDistribution(this);
	Proposal = new ProposalDistribution(this);
	Prior = new PriorDistribution(this);
}

MCMCConfiguration::MCMCConfiguration(string run_file, string configuration){
	configname = configuration;
	dir_name = run_file;
	parameterfile = run_file+"/mcmc/parameters/" + configuration;
	cout << "In config: " << parameterfile << endl;
	parameter_file_name = parameterfile + "/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	
	EmulatorParams = parameter::getS(parmap, "EMULATOR_PARAMETERS", "");
	ParamNames = parameter::getVS(parmap, "PARAMETER_NAMES", "");
	vector<string> temp_logparam = parameter::getVS(parmap, "LOG_PARAMETERS", "");
	
	
	for(int i =0; i<temp_logparam.size(); i++){
		if(strcmp(temp_logparam[i].c_str(), "true") == 0 || strcmp(temp_logparam[i].c_str(), "True") == 0){
			LogParam.push_back(true);
		}else if(strcmp(temp_logparam[i].c_str(), "false") == 0 || strcmp(temp_logparam[i].c_str(), "False") == 0){
			LogParam.push_back(false);
		}else{
			cout << "Unrecognized LogParam value " << temp_logparam[i] << endl;
			exit(1);
		}
	}
	
	randnum = new CRandom(1234);
	Likelihood = new LikelihoodDistribution(this);
	Proposal = new ProposalDistribution(this);
	Prior = new PriorDistribution(this);
}

MCMCConfiguration::~MCMCConfiguration(){
	
}

MCMCRun::MCMCRun(MCMCConfiguration *mcmc_config){
	mcmcconfig = mcmc_config;
	local_parmap = mcmcconfig->parmap;
	tracedir = mcmcconfig->dir_name + "/mcmc/trace/" + mcmcconfig->configname;
	
	MAXITERATIONS = parameter::getI(local_parmap, "MAX_ITERATIONS", 500);
	WRITEOUT = parameter::getI(local_parmap, "WRITEOUT", 100);
	VIZTRACE = parameter::getB(local_parmap, "VISUALIZE_TRACE", true);
	APPEND_TRACE = parameter::getB(local_parmap, "APPEND_TRACE", false);
	
	ThetaList = new ParameterSetList(this);
	
	if(VIZTRACE){
		Visualizer = new VizHandler(this);
		Viz_Count = parameter::getI(local_parmap, "VIZ_COUNT", floor(MAXITERATIONS/200));
		Visualizer->UpdateTraceFig();
	}
	
	if(APPEND_TRACE){
		string addon = "";
		bool Done = false;
		int filecount = 0;
		while(!Done){
			struct stat st;
			stringstream ss;
			string tempfile = tracedir + addon;
			if(stat(tempfile.c_str(), &st)==0){
				//directory exists.
				cout << tempfile << " exists, trying next option..." << endl;
				filecount++;
				ss << "_" << filecount;
				addon = ss.str();
				ss.str(string());
			}else{
				//doesn't exist
				Done = true;
				tracedir = tempfile;
			}
		}
	}else{
		cout << "Deleting prior trace data." << endl;
		string cmd = "rm " + tracedir + "/output*.dat " + tracedir + "/trace.dat";
		system(cmd.c_str());
	}
	
	string command = "mkdir -p "+ tracedir;
	
	system(command.c_str());
	printf("Iteration\tAlpha\tResult\n");
}

MCMCRun::MCMCRun(MCMCConfiguration *mcmc_config, ParameterSet Theta0){
	mcmcconfig = mcmc_config;
	local_parmap = mcmcconfig->parmap;
	tracedir = mcmcconfig->dir_name + "/mcmc/trace/" + mcmcconfig->configname;
	
	MAXITERATIONS = parameter::getI(local_parmap, "MAX_ITERATIONS", 500);
	WRITEOUT = parameter::getI(local_parmap, "WRITEOUT", 100);
	VIZTRACE = parameter::getB(local_parmap, "VISUALIZE_TRACE", true);
	APPEND_TRACE = parameter::getB(local_parmap, "APPEND_TRACE", false);
	
	ThetaList = new ParameterSetList(this, Theta0);
	
	if(VIZTRACE){
		Visualizer = new VizHandler(this);
		Viz_Count = parameter::getI(local_parmap, "VIZ_COUNT", floor(MAXITERATIONS/200));
		Visualizer->UpdateTraceFig();
	}
	
	if(APPEND_TRACE){
		string addon = "";
		bool Done = false;
		int filecount = 0;
		while(!Done){
			struct stat st;
			stringstream ss;
			string tempfile = tracedir + addon;
			if(stat(tempfile.c_str(), &st)==0){
				//directory exists.
				cout << tempfile << " exists, trying next option..." << endl;
				filecount++;
				ss << "_" << filecount;
				addon = ss.str();
				ss.str(string());
			}else{
				//doesn't exist
				Done = true;
				tracedir = tempfile;
			}
		}
	}else{
		cout << "Deleting prior trace data." << endl;
		string cmd = "rm " + tracedir + "/output*.dat " + tracedir + "/trace.dat";
		system(cmd.c_str());
	}
	
	string command = "mkdir -p "+ tracedir;
	
	system(command.c_str());
	printf("Iteration\tAlpha\tResult\n");
}

MCMCRun::~MCMCRun(){

}

void MCMCRun::Run(){
	double Likelihood_Current,Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	
	double LOGBF,alpha;
	ParameterSet * ThetaZeroPtr = ThetaList->Theta[0];
	ParameterSet CurrentParameters = *ThetaZeroPtr;
	
	if(!ThetaZeroPtr){
		cout << "Error getting zeroth iteration parameters." << endl;
	}
	
	Likelihood_Current = mcmcconfig->Likelihood->Evaluate(*ThetaZeroPtr);
	cout << "Likelihood theta0 eval." << endl;
	Proposal_Current = mcmcconfig->Proposal->Evaluate(*ThetaZeroPtr);
	cout << "Propr theta0 " << endl;
	Prior_Current = mcmcconfig->Prior->Evaluate(*ThetaZeroPtr);
	cout << "Prior" << endl;
	Accept_Count = 0;
	for(int i =1; i<=MAXITERATIONS; i++){
		LOGBF = 0;
		ParameterSet Temp_Theta = mcmcconfig->Proposal->Iterate(CurrentParameters);
		Likelihood_New = mcmcconfig->Likelihood->Evaluate(Temp_Theta);
		Prior_New = mcmcconfig->Prior->Evaluate(Temp_Theta);
		Proposal_New = mcmcconfig->Proposal->Evaluate(Temp_Theta);
		

		if(mcmcconfig->LOGLIKE){
			LOGBF += (Likelihood_New-Likelihood_Current);
		}
		else{
			LOGBF +=log(Likelihood_New/Likelihood_Current);
		}
		if(mcmcconfig->LOGPRIOR){
			LOGBF += (Prior_New-Prior_Current);
		}else
		{
			LOGBF +=log(Prior_New/Prior_Current);
		}
		if(mcmcconfig->LOGPROPOSAL){
			LOGBF += (Proposal_New-Proposal_Current);
		}else
		{
			LOGBF +=log(Proposal_New/Proposal_Current);
		}
		
		alpha = min(1.0,exp(LOGBF));
		printf("%5d\t%.5f\t",i,alpha);
		if(alpha > mcmcconfig->randnum->ran()){ //Accept the proposed set.
			printf("Accept\n");
			Accept_Count++;
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			Proposal_Current = Proposal_New;
			CurrentParameters = Temp_Theta;
		}else{
			printf("Reject\n");
		}
		
		ThetaList->Add(CurrentParameters);
		
		if(VIZTRACE){
			if((i+1) % Viz_Count == 0){
				Visualizer->UpdateTraceFig();
			}
		}
		if((i+1) % WRITEOUT == 0){
			cout << "Writing out." << endl;
			if(VIZTRACE){
				Visualizer->UpdateTraceFig();
			}
			ThetaList->WriteOut();
		}
	}
	ThetaList->WriteOut();
	if(VIZTRACE){
		Visualizer->FinalTrace();
	}
	cout << "Accepts: " << Accept_Count << endl;
	cout << "Iterations: " << MAXITERATIONS << endl;
	cout << "Acceptance ratio: " << (double)Accept_Count/(double)MAXITERATIONS << endl;
}
#endif