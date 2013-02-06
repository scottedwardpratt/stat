#ifndef __PROPOSAL_CC__
#define __PROPOSAL_CC__ 

#include "mcmc.h"

using namespace std;

ProposalDistribution::ProposalDistribution(MCMC * mcmc_in){
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
	PREFACTOR = parameter::getD(*parmap, "PREFACTOR", 1.0);
	SCALE = parameter::getD(*parmap, "SCALE", 1.0);
	MIN = parameter::getD(*parmap, "MIN", 0.0);
	MAX = parameter::getD(*parmap, "MAX", 1.0);
	
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

vector<double> ProposalDistribution::Iterate(vector<double> current, float scale){
	vector<double> proposed = current;
	if(SymmetricProposal){
		//We use the scale set in the parameter file
		
		for(int i=0; i<proposed.size(); i++){
			proposed[i] = (current[i] - mcmc->Min_Ranges[i])/(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]); //scale to between 0 and 1
			proposed[i] = proposed[i] + PREFACTOR*gsl_ran_gaussian(randy, SCALE*MixingStdDev[i]);
			proposed[i] = proposed[i] - floor(proposed[i]); //If it isn't in the range of [0,1], this will force it to be
			proposed[i] = (proposed[i]*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]))+mcmc->Min_Ranges[i]; //Scale back to original range
		}	

	} else {
		// We use whatever scale we just got passed from the rest of the code and the range set in the parameter file
		scale = (scale*(MAX-MIN))+MIN;
		
		for(int i=0; i<proposed.size(); i++){
			proposed[i] = (current[i] - mcmc->Min_Ranges[i])/(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]); //scale to between 0 and 1
			proposed[i] = proposed[i] + PREFACTOR*gsl_ran_gaussian(randy, scale*MixingStdDev[i]);
			proposed[i] = proposed[i] - floor(proposed[i]); //If it isn't in the range of [0,1], this will force it to be
			proposed[i] = (proposed[i]*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]))+mcmc->Min_Ranges[i]; //Scale back to original range
		}	
	}

	/*double diff = 0;
	for(int i=0; i<proposed.size(); i++){
		diff += (proposed[i]-current[i])*(proposed[i]-current[i])/((mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i])*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i]));
	}
	cout << "Average difference between points was: " << sqrt(diff)/proposed.size() << endl;*/

	return proposed;
}

double ProposalDistribution::Evaluate(vector<double> Theta1, vector<double> Theta2, float scale){
	// At the moment we are using a gaussian proposal distribution, so the scale is the standard deviation
	double probability;
	double exponent = 0, prefactor = 1;
	
	if(SymmetricProposal){
		// If it's symmetric this doesn't matter
		probability = 1.0;
	} else {
		scale = (scale*(MAX-MIN))+MIN;
		for(int i=0; i<Theta1.size(); i++){
			exponent += -(Theta1[i]-Theta2[i])*(Theta1[i]-Theta2[i])/(2*scale*scale*MixingStdDev[i]*MixingStdDev[i]);
			prefactor = prefactor/(scale*sqrt(2*M_PI));
		}
		/*cout << exponent << " " << prefactor << endl;*/
		probability = prefactor*exp(exponent);
	}
	
	return probability;
}
#endif