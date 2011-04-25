#ifndef __PARAMETER_CC__
#define __PARAMETER_CC__ 

#include "parameter.h"
#include <sstream>

using namespace std;

ParameterSet::ParameterSet(ParameterSetList *list){
	Used = false;
	paramlist = list;
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
			cout << "Belongs to: " << paramlist->mcmc->dir_name << endl;
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
	GetTheta0();
	string command = "mkdir -p " + mcmc->dir_name + "/mcmc/trace";
	system(command.c_str());
	
	HoldOver = new ParameterSet(this);
}

void ParameterSetList::GetTheta0(){
	parameterMap parmap;
	string theta0_filename = mcmc->dir_name+"/mcmc/parameters/theta0.param";
	ParameterSet temp_set(this);
	parameter::ReadParsFromFile(parmap, theta0_filename);
	vector<string> temp_names = parameter::getVS(parmap, "NAMES", "");
	vector<double> temp_values = parameter::getV(parmap, "VALUES", "");
	vector<string> temp_logparam = parameter::getVS(parmap, "LOGPARAM", "");
	EmulatorParams = parameter::getS(parmap, "EMULATOR", "");
	
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
	
	ParamNames = temp_names;
	temp_set.Names = temp_names;
	temp_set.Values = temp_values;
	
	Add(temp_set);
}

void ParameterSetList::PrintDataToFile(){
	stringstream ss;
	ss << WriteOutCounter+1;
	string filename = mcmc->dir_name + "/mcmc/trace/" + mcmc->runnickname + "/output"+ ss.str() +".dat";
	
	ofstream outputfile;
	
	outputfile.open(filename.c_str());
	
	cout << "Writing out to: " << filename << endl;
	
	if(outputfile){
		outputfile << "#ITERATION";
		for(int i = 0; i<Theta[0]->Names.size(); i++){
			outputfile << Theta[0]->Names[i] << "\t";
		}
		outputfile << endl;
		
		for(int i =0; i < mcmc->WRITEOUT; i++){
			outputfile << i+WriteOutCounter*mcmc->WRITEOUT << "\t";
			for(int j=0; j< Theta[i]->Values.size(); j++){
				if(Theta[i]){
					outputfile << Theta[i]->Values[j] << "\t";
				}
				else{
					cout << "Error: Accessing empty element." << endl;
					exit(1);
				}
			}
			outputfile << endl;
		}
		outputfile.close();
	}else{
		cout << "Error in writing output file." << endl;
		exit(1);
	}
	WriteOutCounter++;
	cout << "Done printing to file." << endl;
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
	// cout << "WriteOut." << endl;
	HoldOver->Reset();
	// cout << "Reset." << endl;
	cout << CurrentIteration << endl;
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
#endif