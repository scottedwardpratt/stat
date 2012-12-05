#ifndef __PARAMETER_CC__
#define __PARAMETER_CC__ 

#include "parameterset.h"
#include <sstream>
#include <time.h>
#include <fstream>

using namespace std;

ParameterSet::ParameterSet(ParameterSetList *list){
	Used = false;
	paramlist = list;
}

ParameterSet::ParameterSet(){
	Used = false;
	paramlist = NULL;
}

void ParameterSet::Initialize(vector<string> names, vector<double> values){
	Names = names;
	Values = values;
	Used = true;
}

void ParameterSet::Initialize(ParameterSet ParamSetIn){
	Names = ParamSetIn.Names;
	Values = ParamSetIn.Values;
	Used = true;
}

void ParameterSet::Print(){
	if(Used){
		cout << "This parameter set is alive." << endl;
	}else{
		cout << "This parameter set isn't alive." << endl;
	}
	cout << "This parameter set contains " << Names.size() << " parameters." << endl;
	
	cout << "And "<< Values.size() << " values" << endl;
	
	for(int i=0; i < Names.size(); i++){
		cout << Names[i] << ":\t" << Values[i] << endl;
	}
	
	if(paramlist != NULL){
		if(paramlist->mcmc != NULL){
			cout << "Belongs to: " << paramlist->mcmc->tracedir << endl;
		}else{
			cout << "Parameter list doesn't have a MCMC run associated with it." << endl;
		}
	}else{
		cout << "No parameter list!" << endl;
	}
}

void ParameterSet::Reset(){
	Names.clear();
	Values.clear();
	Used = false;
	InTrace = false;
}

int ParameterSet::GetIndex(string ParamName){
	int out = -1;
	int i = 0;
	bool Found = false;
	
	while(i < Names.size()){
		//cout << "FindParam: Comparing " << name << " to " << PNames[i] << endl;
		if(strcmp(Names[i].c_str(), ParamName.c_str()) == 0){
			if(!Found){
				out = i;
				Found = true;
			}else{ //A matching parameter has already been found, multiple parameters with the same name.
				cout << Names[out] << endl;
				cout << Names[i] << endl;
				cout << "In ParameterSet::GetIndex; Duplicate parameter names found. Please change parameter names." << endl;
				exit(1);
			}
		}
		i++;
	}
	return out;
}

double ParameterSet::GetValue(string ParamName){
	int index = GetIndex(ParamName);
	
	return Values[index];
}

void ParameterSet::VizTrace(){
	if(!InTrace){
		InTrace = true;
	}
}

ParameterSetList::ParameterSetList(MCMC *mcmc_in){
	mcmc = mcmc_in;
	Theta = new ParameterSet*[mcmc->WRITEOUT+1];
	for(int i = 0; i < mcmc->WRITEOUT; i++){
		Theta[i] = new ParameterSet(this);
	}
	WriteOutCounter = 0;
	CurrentIteration = 0;
	HoldOver = new ParameterSet(this);
	
	GetTheta0FromFile();
}

ParameterSetList::ParameterSetList(MCMC *mcmc_in, int seed){ //This one creates a random set of theta0s
	mcmc = mcmc_in;
	Theta = new ParameterSet*[mcmc->WRITEOUT+1];
	for(int i = 0; i < mcmc->WRITEOUT; i++){
		Theta[i] = new ParameterSet(this);
	}
	WriteOutCounter = 0;
	CurrentIteration = 0;
	HoldOver = new ParameterSet(this);
	
	GetRandomTheta0(mcmc, seed);
}

ParameterSetList::ParameterSetList(MCMC *mcmc_in, ParameterSet Theta0){
	mcmc = mcmc_in;
	Theta = new ParameterSet*[mcmc->WRITEOUT+1];
	for(int i = 0; i < mcmc->WRITEOUT; i++){
		Theta[i] = new ParameterSet(this);
	}
	WriteOutCounter = 0;
	CurrentIteration = 0;
	HoldOver = new ParameterSet(this);
	
	GetTheta0FromFile();
	Theta[0]->Reset();
	Theta[0]->Initialize(Theta0);
}

void ParameterSetList::GetTheta0FromFile(){
	parameterMap parmap;
	string theta0_filename = mcmc->parameterfile + "/theta0.param";
	ParameterSet temp_set(this);
	parameter::ReadParsFromFile(parmap, theta0_filename);
	vector<string> temp_names = parameter::getVS(parmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(parmap, "VALUES", "");
	
	ParamNames = temp_names;
	temp_set.Names = temp_names;
	temp_set.Values = temp_values;
	
	Add(temp_set);
}

void ParameterSetList::GetRandomTheta0(MCMC *mcmc_in, int seed){ //This one creates random theta 0s
	mcmc = mcmc_in;
	srand(seed);
	parameterMap parmap;
	ParameterSet temp_set(this);
	vector<string> temp_names = mcmc->ParamNames;
	vector<double> temp_values;

	cout << "We are using random theta0 values. They are:" << endl;
	for(int i=0;i<mcmc->ParamNames.size();i++){
		temp_values.push_back((double(rand())/double(RAND_MAX))*(mcmc->Max_Ranges[i]-mcmc->Min_Ranges[i])+mcmc->Min_Ranges[i]);
		//cout << mcmc->ParamNames[i] << " " << temp_values[i] << endl;
	}

	ParamNames = temp_names;
	temp_set.Names = temp_names;
	temp_set.Values = temp_values;
	
	Add(temp_set);
}

void ParameterSetList::PrintDataToFile(){
	cout << "Printing data to file." << endl;
	stringstream ss;
	ss << WriteOutCounter+1;
	string filename = mcmc->tracedir + "/output"+ ss.str() +".dat";
	
	ofstream outputfile;
	if(Theta[0]->Used){
		outputfile.open(filename.c_str());

		cout << "Writing out to: " << filename << endl;

		if(outputfile){
			outputfile << "#ITERATION,";
			for(int i = 0; i<Theta[0]->Names.size(); i++){
				outputfile << Theta[0]->Names[i] << ",";
			}
			// outputfile << "LIKELIHOOD" << endl;

			for(int i =0; i < mcmc->WRITEOUT; i++){
				if(Theta[i]->Used){
					outputfile << i+WriteOutCounter*mcmc->WRITEOUT << ",";
					for(int j=0; j< Theta[i]->Values.size(); j++){
						if(Theta[i]){
							if(mcmc->RESCALED_TRACE){
								outputfile << (Theta[i]->Values[j]-mcmc->Min_Ranges[j])/(mcmc->Max_Ranges[j]-mcmc->Min_Ranges[j]);
								if(j!=Theta[i]->Values.size()-1){
									outputfile << ",";
								}
							}else{
								outputfile << Theta[i]->Values[j];
								if(j!=Theta[i]->Values.size()-1){
									outputfile << ",";
								}
							}
						}
						else{
							cout << "Error: Accessing empty element." << endl;
							exit(1);
						}
					}
					outputfile << endl;
					// outputfile << "," << LikelihoodArray[i] << endl;
				}
			}
			outputfile.close();
		}else{
			cout << "Error in writing output file." << endl;
			exit(1);
		}
		WriteOutCounter++;
		cout << "Done printing to file." << endl;
	}
}

void ParameterSetList::Add(ParameterSet Theta_In){
	if(CurrentIteration >= mcmc->WRITEOUT){
		cerr << "Error: Thetalist index out of bounds (Greater than WRITEOUT)." << endl;
		exit(1);
	}else{
		Theta[CurrentIteration]->Initialize(Theta_In);
		CurrentIteration++;
	}
}

void ParameterSetList::WriteOut(){
	HoldOver->Reset();
	HoldOver->Initialize(*Theta[--CurrentIteration]);
	PrintDataToFile();
	for(int i = 0; i < mcmc->WRITEOUT; i++){
		Theta[i]->Reset();
	}
	CurrentIteration = 0;
}

ParameterSet ParameterSetList::CurrentParameters(){
	ParameterSet tempset(this);
	if(CurrentIteration == 0){
		tempset = *HoldOver;
	}
	else{
		tempset = *Theta[CurrentIteration -1];
	}
	return tempset;
}

void ParameterSetList::MakeTrace(){
	stringstream ss;
	ss << "cat ";
	
	for(int i = 1; i <=ceil((double)(mcmc->MAXITERATIONS)/(double)(mcmc->WRITEOUT)); i++){
		cout << "Parsing output" << i << ".dat" << endl;
		ss << mcmc->tracedir << "/output" << i << ".dat ";
	}
	ss << "> " << mcmc->tracedir << "/trace.dat" << endl;
	
	string command = ss.str();
	system(command.c_str());
	
	ss.str(string());
}
#endif