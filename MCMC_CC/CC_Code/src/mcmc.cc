#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMC::MCMC(string run_file){
	dir_name = run_file;
	parameter_file_name = run_file+"/mcmc/parameters/mcmc.param";
	
	cout << "Reading in " << parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name.c_str());
	MAXITERATIONS = parameter::getI(parmap, "MAX_ITERATIONS", 500);
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	WRITEOUT = parameter::getI(parmap, "WRITEOUT", 100);
	
	randnum = new CRandom(1234);
	ThetaList = new ParameterSetList(this);
	
	Likelihood = new LikelihoodDistribution(this);
	Proposal = new ProposalDistribution(this);
	Prior = new PriorDistribution(this);
	
	// cout << "Distributions Declared." << endl;
	Accept_Count = 0;
	
	string command = "mkdir -p "+run_file+"/mcmc/trace";
	system(command.c_str());
	printf("Iteration\tAlpha\tResult\n");
	
	// cout << "MCMC Constructor: Done" << endl;
}

void MCMC::Run(){
	double Likelihood_Current,Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	
	double LOGBF,alpha;
	ParameterSet * ThetaZeroPtr = ThetaList->Theta[0];
	
	if(!ThetaZeroPtr){
		cout << "Error getting zeroth iteration parameters." << endl;
	}
	
	Likelihood_Current = Likelihood->Evaluate(*ThetaZeroPtr);
	Proposal_Current = Proposal->Evaluate(*ThetaZeroPtr);
	Prior_Current = Prior->Evaluate(*ThetaZeroPtr);
	
	for(int i =1; i<=MAXITERATIONS; i++){
		cout << "MCMC step " << i << endl;
		LOGBF = 0;
		ParameterSet Temp_Theta = Proposal->Iterate();
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
		// 
		// cout << "New: " << Likelihood_New << endl;
		// cout << "Current: " << Likelihood_Current << endl;
		// cout << log(Likelihood_New/Likelihood_Current) << endl;
		
		alpha = min(1.0,exp(LOGBF));
		printf("%5d\t%5g\t",i,alpha);
		if(alpha > randnum->ran()){ //Accept the proposed set.
			printf("Accept\n");
			Accept_Count++;
			ThetaList->Add(Temp_Theta);
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			Proposal_Current = Proposal_New;
		}else{
			printf("Reject\n");
		}
	}
}
#endif