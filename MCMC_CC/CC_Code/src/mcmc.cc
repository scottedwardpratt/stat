#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMC::MCMC(string run_file){
	dir_name = run_file;
	parameter_file_name = run_file+"/mcmc/parameters/mcmc.param";
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	runnickname = parameter::getS(parmap, "NAME", "mcmcrun");
	MAXITERATIONS = parameter::getI(parmap, "MAX_ITERATIONS", 500);
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	WRITEOUT = parameter::getI(parmap, "WRITEOUT", 100);
	VIZTRACE = parameter::getB(parmap, "VISUALIZE_TRACE", true);
	APPEND_TRACE = parameter::getB(parmap, "APPEND_TRACE", false);
	
	randnum = new CRandom(1234);
	ThetaList = new ParameterSetList(this);
	
	Likelihood = new LikelihoodDistribution(this);
	Proposal = new ProposalDistribution(this);
	Prior = new PriorDistribution(this);
	if(VIZTRACE){
		cout << "Visualizing trace." << endl;
		Viz_Count = parameter::getI(parmap, "VIZ_COUNT", floor(MAXITERATIONS/200));
		Visualizer = new VizHandler(this);
		Visualizer->UpdateTraceFig();
	}
	
	Accept_Count = 0;
	
	tracedir = run_file+"/mcmc/trace/"+ runnickname;
	if(APPEND_TRACE){
		string addon;
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

void MCMC::Run(){
	double Likelihood_Current,Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	
	double LOGBF,alpha;
	ParameterSet * ThetaZeroPtr = ThetaList->Theta[0];
	ParameterSet CurrentParameters = *ThetaZeroPtr;
	
	if(!ThetaZeroPtr){
		cout << "Error getting zeroth iteration parameters." << endl;
	}
	
	Likelihood_Current = Likelihood->Evaluate(*ThetaZeroPtr);
	Proposal_Current = Proposal->Evaluate(*ThetaZeroPtr);
	Prior_Current = Prior->Evaluate(*ThetaZeroPtr);
	
	Accept_Count = 0;
	for(int i =1; i<=MAXITERATIONS; i++){
		LOGBF = 0;
		ParameterSet Temp_Theta = Proposal->Iterate(CurrentParameters);
		Likelihood_New = Likelihood->Evaluate(Temp_Theta);
		Prior_New = Prior->Evaluate(Temp_Theta);
		Proposal_New= Proposal->Evaluate(Temp_Theta);
		
		if(LOGLIKE){
			LOGBF += (Likelihood_New-Likelihood_Current);
		}
		else{
			LOGBF +=log(Likelihood_New/Likelihood_Current);
		}
		if(LOGPRIOR){
			LOGBF += (Prior_New-Prior_Current);
		}else
		{
			LOGBF +=log(Prior_New/Prior_Current);
		}
		if(LOGPROPOSAL){
			LOGBF += (Proposal_New-Proposal_Current);
		}else
		{
			LOGBF +=log(Proposal_New/Proposal_Current);
		}
		
		alpha = min(1.0,exp(LOGBF));
		printf("%5d\t%.5f\t",i,alpha);
		if(alpha > randnum->ran()){ //Accept the proposed set.
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
	Visualizer->FinalTrace();
	cout << "Accepts: " << Accept_Count << endl;
	cout << "Iterations: " << MAXITERATIONS << endl;
	cout << "Acceptance ratio: " << (double)Accept_Count/(double)MAXITERATIONS << endl;
}
#endif