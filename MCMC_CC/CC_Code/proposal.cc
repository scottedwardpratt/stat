#ifndef __PROPOSAL_CC__
#define __PROPOSAL_CC__ 

#include "mcmc.h"

using namespace std;

ProposalDistribution::ProposalDistribution(MCMC * mcmc_in){
	fstream ranges;
	string type, param_name;
	int count = 0;
	
	mcmc = mcmc_in;
	int numparams = mcmc->ThetaList->ParamNames.size();
	vector<double> temp (numparams, .01);
	MixingStdDev = parameter::getV(mcmc->parmap, "MIXING_STD_DEV", 0);
	Ranges[1]=new double[numparams];
	Ranges[2]=new double[numparams];
	
	
	
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



#endif