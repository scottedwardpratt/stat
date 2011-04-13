#ifndef __MCMC_CC__
#define __MCMC_CC__

#include "mcmc.h"

using namespace std;

MCMC::MCMC(string run_file){
	dir_name = run_file;
	parameter_file_name = run_file+"/parameters/mcmc.param";
	
	cout << "Reading in " <<parameter_file_name << endl;
	
	parameter::ReadParsFromFile(parmap, parameter_file_name);
	MAXITERATIONS = parameter::getI(parmap, "MAX_ITERATIONS", 500);
	LOGLIKE = parameter::getB(parmap, "LOGLIKE", true);
	LOGPRIOR = parameter::getB(parmap, "LOGPRIOR", true);
	LOGPROPOSAL = parameter::getB(parmap, "LOGPROPOSAL", true);
	WRITEOUT = parameter::getI(parmap, "WRITEOUT", 100);
	
	//Likelihood = new LikelihoodDistribution;
	Proposal = new ProposalDistribution(this);
	//Prior = new PriorDistribution;
	
	Accept_Count = 0;
	
	randnum = new CRandom(1234);
	ThetaList = new ParameterSetList(this);
	
	string command = "mkdir -p "+run_file+"/output/mcmc";
	system(command.c_str());
	printf("Iteration\tAlpha\tResult\n");
}

void MCMC::Run(){
	double Likelihood_Current,Likelihood_New;
	double Prior_Current, Prior_New;
	double Proposal_Current, Proposal_New;
	
	double LOGBF,alpha;
	
	// Likelihood_Current = Likelihood.Evaluate(ThetaList->Theta[0]);
	// Prior_Current = Prior.Evaluate(ThetaList->Theta[0]);
	Proposal_Current = Proposal->Evaluate(ThetaList->Theta[0]);
	
	
	for(int i =1; i<=MAXITERATIONS; i++){
		LOGBF = 0;
		ParameterSet Temp_Theta = Proposal->Iterate();
		// Likelihood_New = Likeliood.Evaluate(Temp_Theta);
		// Prior_New = Prior.Evaluate(Temp_Theta);
		Proposal_New= Proposal->Evaluate(Temp_Theta);
		
		if(LOGLIKE){
			LOGBF += Likelihood_New-Likelihood_Current;
		}else{
			LOGBF +=log(Likelihood_New/Likelihood_Current);
		}
		if(LOGPRIOR){
			LOGBF += Prior_New-Prior_Current;
		}else{
			LOGBF +=log(Prior_New/Prior_Current);
		}
		if(LOGPROPOSAL){
			LOGBF += Proposal_New-Proposal_Current;
		}else{
			LOGBF +=log(Proposal_New/Proposal_Current);
		}
		
		alpha = min(1.0,exp(LOGBF));
		printf("%5d\t%5g\t",i,alpha);
		if(alpha > randnum->ran()){ //Accept the proposed set.
			Accept_Count++;
			ThetaList->Add(Temp_Theta);
			Likelihood_Current = Likelihood_New;
			Prior_Current = Prior_New;
			Proposal_Current = Proposal_New;
			printf("Accept\n");
		}else{
			printf("Reject\n");
		}
	}
}
#endif