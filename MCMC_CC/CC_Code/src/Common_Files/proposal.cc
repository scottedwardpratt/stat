#ifndef __PROPOSAL_CC__
#define __PROPOSAL_CC__ 

#include "mcmc.h"

using namespace std;

ProposalDistribution::ProposalDistribution(MCMCConfiguration * mcmc_in){
	mcmc=mcmc_in;
	string type, param_name;
	int count = 0;
	int numparams = mcmc->ParamNames.size();
	vector<double> temp (numparams, .01);
	
	SepMap = parameter::getB(mcmc->parmap, "PROPOSAL_PARAMETER_MAP", false);
	
	if(SepMap){
		string parmapfile = mcmc->parameterfile + "/proposal.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
	}else{
		parmap = &(mcmc->parmap);
	}
	
	Rescaled_Method = parameter::getB(*parmap, "RESCALED_PROPOSAL", true);
	MixingStdDev = parameter::getV(*parmap, "MIXING_STD_DEV", "0");
	SymmetricProposal = parameter::getB(*parmap, "SYMMETRIC_PROPOSAL", true);
	SCALE = parameter::getD(*parmap, "SCALE", 1.0);
	
	const gsl_rng_type * rngtype;
	rngtype = gsl_rng_default;
	gsl_rng_env_setup();
	randy = gsl_rng_alloc(rngtype);
	gsl_rng_set(randy, time(NULL));

	
}

int ProposalDistribution::FindParam(string name){
	vector<string> PNames = mcmc->ParamNames;
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

ParameterSet ProposalDistribution::Iterate(ParameterSet current, float scale){
	if(SymmetricProposal){
		//We use the scale set in the parameter file
		ParameterSet proposed = current;
		
		for(int i=0; i<proposed.Values.size(); i++){
			proposed.Values[i] = (current.Values[i] - mcmc->Min_Ranges[i])/(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]); //scale to between 0 and 1
			//proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, SCALE*MixingStdDev[i]/sqrt((double)proposed.Names.size()));
			proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, SCALE*MixingStdDev[i]);
			proposed.Values[i] = proposed.Values[i] - floor(proposed.Values[i]);
			proposed.Values[i] = (proposed.Values[i]*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]))+mcmc->Min_Ranges[i];
		}	

		return proposed;
	} else {
		// We use whatever scale we just got passed from the rest of the code
		ParameterSet proposed = current;
		
		for(int i=0; i<proposed.Values.size(); i++){
			proposed.Values[i] = (current.Values[i] - mcmc->Min_Ranges[i])/(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]); //scale to between 0 and 1
			//proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, scale*MixingStdDev[i]/sqrt((double)proposed.Names.size()));
			proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, SCALE*scale*MixingStdDev[i]);
			proposed.Values[i] = proposed.Values[i] - floor(proposed.Values[i]);
			proposed.Values[i] = (proposed.Values[i]*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]))+mcmc->Min_Ranges[i];
		}	

		return proposed;
	}
}

double ProposalDistribution::Evaluate(ParameterSet Theta1, ParameterSet Theta2, float scale){
	// At the moment we are using a gaussian proposal distribution, so the scale is the standard deviation
	double probability;
	double exponent = 0, prefactor = 1;
	
	if(SymmetricProposal){
		// If it's symmetric this doesn't matter
		probability = 1.0;
	}else{
		// We need to actually calculate the proposal
		/*cout << "ERROR: ProposalDistribution::Evaluate:" << endl;
		cout << "Allowed for unsymmetric proposal without defining method to evaluate proposal." << endl;
		exit(-1);*/
		for(int i=0; i<Theta1.Values.size(); i++){
			exponent += -(Theta1.Values[i]-Theta2.Values[i])*(Theta1.Values[i]-Theta2.Values[i])/(2*SCALE*SCALE*scale*scale*MixingStdDev[i]*MixingStdDev[i]);
			prefactor = prefactor/(SCALE*scale*sqrt(2*M_PI));
		}
		/*cout << scale << endl;
		cout << exponent << " " << prefactor << endl;
		probability = prefactor*exp(exponent);*/
	}
	
	return probability;
}
#endif