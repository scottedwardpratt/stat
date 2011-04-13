#ifndef __PROPOSAL_CC__
#define __PROPOSAL_CC__ 

#include "mcmc.h"

using namespace std;

ProposalDistribution::ProposalDistribution(MCMC * mcmc_in): Distribution(mcmc_in){
	fstream ranges;
	string type, param_name;
	int count = 0;
	
	randy=gsl_rng_alloc(gsl_rng_ranlxd1);
	int numparams = mcmc->ThetaList->ParamNames.size();
	vector<double> temp (numparams, .01);
	MixingStdDev = parameter::getV(mcmc->parmap, "MIXING_STD_DEV", 0);
	Ranges[1]=new double[numparams];
	Ranges[2]=new double[numparams];
	
	SymmetricProposal = parameter::getB(mcmc->parmap, "SYMMETRIC_PROPOSAL", true);
	
	ranges.open("ranges.dat", fstream::in);
	
	if(ranges){
		while(!ranges.eof()){
			ranges >> type;
			if(strcmp(type.c_str(), "double") != 0){
				ranges >> param_name;
				int index = FindParam(param_name);
				if(index != -1){ //returns -1 if not found
					ranges >> Ranges[1][index]; //minimum
					ranges >> Ranges[2][index]; //maximum
					if(Ranges[1][index] > Ranges[2][index]){ //Flip them
						double temp2 = Ranges[1][index];
						Ranges[1][index] = Ranges[2][index];
						Ranges[2][index] = temp2;
					}
					count++;
				}else{
					cout << "Unrecognized parameter name " << param_name << endl;
					exit(1);
				}
			}else{
				cout << "Unrecognized variable type " << type << endl;
				exit(1);
			}
		}
		if(count != numparams){
			cout << "Error: number of parameters in ranges data and total number of parameters are different." << endl;
			exit(1);
		}
	}else{
		cout << "Unable to open ranges.dat" << endl;
		exit(1);
	}
}

int ProposalDistribution::FindParam(string name){
	vector<string> PNames = mcmc->ThetaList->ParamNames;
	int out = -1;
	int i = 0;
	bool Found = false;
	
	while(i < PNames.size()){
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
	}
	return out;
}

ParameterSet ProposalDistribution::Iterate(){
	ParameterSet current = mcmc->ThetaList->CurrentParameters();
	ParameterSet proposed = current;
	
	for(int i=0; i<proposed.Values.size(); i++){
		do{
			proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, MixingStdDev[i]);
		}while(proposed.Values[i]>Ranges[1][i]&&proposed.Values[i]<Ranges[2][i]);
	}
	
	return proposed;
}

double ProposalDistribution::Evaluate(ParameterSet Theta){
	double probability;
	
	if(SymmetricProposal){
		probability = 1;
	}else{
		//do something else.
	}
	
	return probability;
}
#endif