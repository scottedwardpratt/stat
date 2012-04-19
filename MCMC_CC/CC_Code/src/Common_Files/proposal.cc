#ifndef __PROPOSAL_CC__
#define __PROPOSAL_CC__ 

#include "mcmc.h"

using namespace std;

ProposalDistribution::ProposalDistribution(MCMCConfiguration * mcmc_in){
	mcmc=mcmc_in;
	fstream ranges;
	string type, param_name;
	int count = 0;
	bool Ranges = false;
	int numparams = mcmc->ParamNames.size();
	Min_Ranges=new double[numparams];
	Max_Ranges=new double[numparams];
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
	
	//Determine the acceptable ranges of the parameters
	string filename = mcmc->dir_name + "/ranges.dat";
	ranges.open(filename.c_str(), fstream::in);
	
	if(ranges){
		Ranges = true;
		while(ranges >> type){
			if(strcmp(type.c_str(), "double") == 0){
				ranges >> param_name;
				int index = FindParam(param_name);
				if(index != -1){ //returns -1 if not found
					ranges >> Min_Ranges[index]; //minimum
					ranges >> Max_Ranges[index]; //maximum
					if(Min_Ranges[index] > Max_Ranges[index]){ //Flip them
						double temp2 = Min_Ranges[index];
						Min_Ranges[index] = Max_Ranges[index];
						Max_Ranges[index] = temp2;
					}
					count++;
				}else{
					cout << "Unrecognized parameter name " << param_name << endl;
					exit(1);
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(ranges, temp, '\n');
				}
				else{
					cout << "Unrecognized variable type " << type << endl;
					exit(1);
				}
			}
		}
		ranges.close();
	}else{
		cout << "Warning: Unable to open ranges.dat in model directory." << endl;
		// exit(1);
	}
	if(count != numparams){
		cout << "Error: number of parameters in ranges data and total number of parameters are different." << endl;
		cout << "Count: " << count << " NumParams: " << numparams << endl;
		exit(1);
	}
	/* //This seems to be a copy of what is above, so I have commented it out. JN 4/19/12
	count=0;
	string emulator_ranges = mcmc->dir_name + "/ranges.dat";
	ranges.open(emulator_ranges.c_str(), fstream::in);
	if(ranges){
		Ranges = true;
		while(ranges >> type){
			if(strcmp(type.c_str(), "double") == 0){
				ranges >> param_name;
				int index = FindParam(param_name);
				if(index != -1){ //returns -1 if not found
					ranges >> Min_Ranges[index]; //minimum
					ranges >> Max_Ranges[index]; //maximum
					if(Min_Ranges[index] > Max_Ranges[index]){ //Flip them
						double temp2 = Min_Ranges[index];
						Min_Ranges[index] = Max_Ranges[index];
						Max_Ranges[index] = temp2;
					}
					count++;
				}else{
					cout << "Unrecognized parameter name " << param_name << endl;
					exit(1);
				}
			}else{
				if(strncmp(type.c_str(), "#", 1) == 0){
					string temp;
					getline(ranges, temp, '\n');
				}
				else{
					cout << "Unrecognized variable type " << type << endl;
					exit(1);
				}
			}
		}
		ranges.close();
	}
	else{
		cout << "Warning:Unable to open /ranges.dat" << endl;
	}
	*/
	if(!Ranges){
		cout << "Unable to find any ranges.dat file." << endl;
		exit(-1);
	}
	if(count != numparams){
		cout << "Error: number of parameters in ranges data and total number of parameters are different." << endl;
		cout << "Count: " << count << " NumParams: " << numparams << endl;
		exit(1);
	}

	
	
	//gsl_rng * r;
	
	/*
	const gsl_rng_type * rngtype;
	rngtype = gsl_rng_default;
	gsl_rng_env_setup();
	randy = gsl_rng_alloc(rngtype);
	gsl_rng_set(randy, time(NULL));
*/
	
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

ParameterSet ProposalDistribution::Iterate(ParameterSet current){
	ParameterSet proposed = current;
	
	if(!Rescaled_Method){
		if(proposed.Values.size() != proposed.Names.size() || proposed.Values.size() != MixingStdDev.size()){
			cout << "Error: Parameter names/values aren't same size, or no mixing standard deviation for a parameter." << endl;
		}

		for(int i=0; i<proposed.Values.size(); i++){
			do{
				proposed.Values[i] = current.Values[i];
				if(i == proposed.Names.size() || i == MixingStdDev.size()){
					cout << "Error: iterating over non-existant parameter name or standard deviation." << endl;
					cout << "(Seg-Fault catch)" << endl;
					exit(-1);
				}
				/* 
				 added factor of 0.025 below for testing
				 */
				//proposed.Values[i] = proposed.Values[i] + 0.025*gsl_ran_gaussian(randy, MixingStdDev[i]);
				//proposed.Values[i] = proposed.Values[i] + 0.5*gsl_ran_gaussian(randy, MixingStdDev[i]);
				//proposed.Values[i] = proposed.Values[i] + 2*gsl_ran_gaussian(randy, MixingStdDev[i]);
				proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, MixingStdDev[i]);
			}while((proposed.Values[i] < Min_Ranges[i]) || (proposed.Values[i]>Max_Ranges[i]));
		}
	}
	else{
		for(int i=0; i<proposed.Values.size(); i++){
			do{
				proposed.Values[i] = current.Values[i];
				proposed.Values[i] = (proposed.Values[i] - Min_Ranges[i])/(Max_Ranges[i]-Min_Ranges[i]); //scale to between 0 and 1
				
				proposed.Values[i] = proposed.Values[i] + gsl_ran_gaussian(randy, MixingStdDev[0]/sqrt((double)proposed.Names.size()));
			}while((proposed.Values[i] < 0.0) || proposed.Values[i] > 1.0);
			
			proposed.Values[i] = (proposed.Values[i]*(Max_Ranges[i]-Min_Ranges[i]))+Min_Ranges[i];
		}
	}
	
	return proposed;
}

double ProposalDistribution::Evaluate(ParameterSet Theta){
	double probability = 0.01;
	
	if(SymmetricProposal){
		probability = 1.0;
	}else{
		cout << "ERROR: ProposalDistribution::Evaluate:" << endl;
		cout << "Allowed for unsymmetric proposal without defining method to evaluate proposal." << endl;
		exit(-1);
	}
	
	return probability;
}
#endif