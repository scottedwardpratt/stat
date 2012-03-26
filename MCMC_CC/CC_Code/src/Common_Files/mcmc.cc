#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMCConfiguration::MCMCConfiguration(string info_dir){
	configname = "default";
	dir_name = info_dir;
	parameterfile = info_dir+"/parameters/" + configname;
	cout << "In config: " << parameterfile << endl;
	parameter_file_name = parameterfile + "/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	
	EmulatorParams = parameter::getS(parmap, "EMULATOR_PARAMETERS", "");
	ParamNames = parameter::getVS(parmap, "PARAMETER_NAMES", "blah blah blah");
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
	Proposal = new ProposalDistribution(this);
	
	if(parameter::getS(parmap,"MODEL","NOMODEL")=="CosmoSurvey"){
		Likelihood = new LikelihoodDistribution_Cosmo(this);
		Prior = new PriorDistribution_Cosmo(this);
	}
	else if(parameter::getS(parmap,"MODEL","NOMODEL")=="RHIC"){
	 Likelihood = new LikelihoodDistribution_RHIC(this);
	 Prior = new PriorDistribution_RHIC(this);
	 }
	else if((parameter::getS(parmap,"MODEL","NOMODEL")=="TEST")||(parameter::getS(parmap,"MODEL","NOMODEL")=="Test")){
	 Likelihood = new LikelihoodDistribution_Test(this);
	 Prior = new PriorDistribution_Test(this);
	 }
	else{
		printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
		exit(1);
	}
}

MCMCConfiguration::MCMCConfiguration(string info_dir, string configuration){
	configname = configuration;
	dir_name = info_dir;
	parameterfile = info_dir+"/parameters/" + configuration;
	cout << "In config: " << parameterfile << endl;
	parameter_file_name = parameterfile + "/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	
	EmulatorParams = parameter::getS(parmap, "EMULATOR_PARAMETERS", "");
	ParamNames = parameter::getVS(parmap, "PARAMETER_NAMES", "blah blah blah");
	for(int i=0;i<ParamNames.size();i++) printf("%s ",ParamNames[i].c_str());
	printf("\n");
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
	
	// cout << "stuff done." << endl;
	randnum = new CRandom(1234);
	if(parameter::getS(parmap,"MODEL","NOMODEL")=="CosmoSurvey"){
		Likelihood = new LikelihoodDistribution_Cosmo(this);
		Prior = new PriorDistribution_Cosmo(this);
	}
	else if(parameter::getS(parmap,"MODEL","NOMODEL")=="RHIC"){
	 Likelihood = new LikelihoodDistribution_RHIC(this);
	 Prior = new PriorDistribution_RHIC(this);
	 }
	else if((parameter::getS(parmap,"MODEL","NOMODEL")=="TEST")||(parameter::getS(parmap,"MODEL","NOMODEL")=="Test")){
	 Likelihood = new LikelihoodDistribution_Test(this);
	 Prior = new PriorDistribution_Test(this);
	 }
	else{
		printf("Must define parameter MODEL in parameter file, or yours is unrecognized\n");
		exit(1);
	}
	/*Likelihood = new LikelihoodDistribution(this);
	 cout << "Like done." << endl; 
	Prior = new PriorDistribution(this);
	 cout << "Prior done." << endl;*/
	
	Proposal = new ProposalDistribution(this);
	//cout << "Proposal done." << endl;
	
}

MCMCConfiguration::~MCMCConfiguration(){
	
}

MCMCRun::MCMCRun(MCMCConfiguration *mcmc_config){
	mcmcconfig = mcmc_config;
	local_parmap = mcmcconfig->parmap;
	tracedir = mcmcconfig->dir_name + "/trace/" + mcmcconfig->configname;
	
	MAXITERATIONS = parameter::getI(local_parmap, "MAX_ITERATIONS", 500);
	WRITEOUT = parameter::getI(local_parmap, "WRITEOUT", 100);
	VIZTRACE = parameter::getB(local_parmap, "VISUALIZE_TRACE", true);
	APPEND_TRACE = parameter::getB(local_parmap, "APPEND_TRACE", false);
	
	ThetaList = new ParameterSetList(this);
	//Likelihood_Current=0;

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
	tracedir = mcmcconfig->dir_name + "/trace/" + mcmcconfig->configname;
	
	MAXITERATIONS = parameter::getI(local_parmap, "MAX_ITERATIONS", 500);
	WRITEOUT = parameter::getI(local_parmap, "WRITEOUT", 100);
	VIZTRACE = parameter::getB(local_parmap, "VISUALIZE_TRACE", true);
	APPEND_TRACE = parameter::getB(local_parmap, "APPEND_TRACE", false);
	
	ThetaList = new ParameterSetList(this, Theta0);
	//Likelihood_Current=0;
	
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

/** This runs MAXITERATIONS samplings */
double MCMCRun::Run(){
	double Likelihood_Current, Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	
	double LOGBF,alpha;
	ParameterSet * ThetaZeroPtr = ThetaList->Theta[0];
	ParameterSet CurrentParameters = *ThetaZeroPtr;
	
	if(!ThetaZeroPtr){
		cout << "Error getting zeroth iteration parameters." << endl;
	}
	
	Likelihood_Current = mcmcconfig->Likelihood->Evaluate(*ThetaZeroPtr);
	Proposal_Current = mcmcconfig->Proposal->Evaluate(*ThetaZeroPtr);
	Prior_Current = mcmcconfig->Prior->Evaluate(*ThetaZeroPtr);

	Accept_Count = 0;
	for(int i =1; i<=MAXITERATIONS; i++){
		LOGBF = 0;
		ParameterSet Temp_Theta = mcmcconfig->Proposal->Iterate(CurrentParameters);
		Likelihood_New = mcmcconfig->Likelihood->Evaluate(Temp_Theta);
		if(i==1){
			bestlikelihood=Likelihood_New;
			BestParameterSetPtr=&CurrentParameters;
		}
		Prior_New = mcmcconfig->Prior->Evaluate(Temp_Theta);
		Proposal_New = mcmcconfig->Proposal->Evaluate(Temp_Theta);
		
		// cout << "Likelihood of proposed set: " << Likelihood_New << endl;
		// cout << "Likelihood of current set: " << Likelihood_Current << endl;
		// cout << "Prior of proposed set: " << Prior_New << endl;
		// cout << "Prior of current set: " << Prior_Current << endl;
		// cout << "Proposal of proposed set: " << Proposal_New;
		// cout << " Proposal of current set: " << Proposal_Current << endl;
		
		if(mcmcconfig->LOGLIKE){
			LOGBF += (Likelihood_New-Likelihood_Current);
			printf(" ll_new=%g, ll_current=%g\n",Likelihood_New,Likelihood_Current);
		}
		else{
			LOGBF +=log(Likelihood_New/Likelihood_Current);
		}
		/*
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
		 */
		
		alpha = min(1.0,exp(LOGBF));
		//alpha = min(0.97,exp(LOGBF));
		// cout << "exp(LOGBF): " << exp(LOGBF) << endl;
		printf("%5d\talpha=%6.5f\t",i,alpha);
		if(alpha > (mcmcconfig->randnum->ran())){ //Accept the proposed set.
		//if(exp(LOGBF) > 1){ //Accept the proposed set.
			printf("Accept\n");
			Accept_Count++;
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			Proposal_Current = Proposal_New;
			CurrentParameters = Temp_Theta;
			if(Likelihood_Current>bestlikelihood && i>1){
				bestlikelihood=Likelihood_New;
				BestParameterSetPtr=&CurrentParameters;
				printf("XXXXXXXXX YIPPEE!! Best parameters so far, likelihood=%g\n",bestlikelihood);
			}
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
	ThetaList->MakeTrace();
	if(VIZTRACE){
		Visualizer->FinalTrace();
	}
	double ratio = (double)Accept_Count/(double)MAXITERATIONS;
	cout << "Accepts: " << Accept_Count << endl;
	cout << "Iterations: " << MAXITERATIONS << endl;
	cout << "Acceptance ratio: " << ratio << endl;
	printf("-------- Best Parameter Set, likelihood=%g -------------\n",bestlikelihood);
	BestParameterSetPtr->Print();
	
	return ratio;
}

#endif