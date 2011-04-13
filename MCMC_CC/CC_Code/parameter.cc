#ifndef __PARAMETER_CC__
#define __PARAMETER_CC__ 

using namespace std;

ParameterSet::ParameterSet(ParameterSetList *list){
	Used = false;
	paramlist = list;
}

ParameterSet::Initialize(vector<string> names, vector<double> values){
	Names = names;
	Values = values;
	Used = true;
}

ParameterSet::Initialize(ParameterSet ParamSetIn){
	Names = ParamSetIn.Names;
	Values = ParamSetIn.Values;
	Used = true;
}

ParameterSetList::ParameterSetList(MCMC *mcmc_in){
	mcmc = mcmc_in;
	Theta = new ParameterSet[mcmc->WRITEOUT](this);
	WriteOutCounter = 1;
	CurrentIteration = 0;
	system("mkdir -p output/mcmc");
	GetTheta0();
}

void ParameterSetList::GetTheta0(){
	parameterMap parmap;
	theta0_filename = mcmc->dirname+"/parameters/theta0.param";
	ParameterSet temp_set;
	parameter::ReadParsFromFile(parmap, theta0_filename);
	vector<string> temp_names = parameter::getVS(parmap, "PARAM_NAMES", NULL);
	vector<double> temp_values = parameter::getV(parmap, "PARAM_VALUES", NULL);
	vector<string> temp_logparam = parameter::getVS(parmap, "PARAM_LOGPARAM", NULL);
	
	for(int i =0; i<temp_logparam.size(); i++){
		if(strcmp(temp_logparam[i], "true") == 0 || strcmp(temp_logparam[i], "True") == 0){
			LogParam.push_back(true);
		}else if(strcmp(temp_logparam[i], "false") == 0 || strcmp(temp_logparam[i], "False") == 0){
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

void ParameterSetList::PrintData(){
	filename = "output/mcmc/ouput"+WriteOutCounter+".dat";
	FILE * outputfile = fopen(filename, "w");
	
	cout << "Writing out to: " << filename << endl;
	if(outputFile != NULL){
		fputs("#", outputfile);
		for(int i = 0; i<Theta[0].Names.size(); i++){
			fputs(Theta[0].Names[i], outputFile);
			fputs("\t", outputFile);
		}
		fputs("\n", outputFile);
		
		for(int j =0; j< Theta.size(); j++){
			for(int k= 0; k<Theta[j].Values.size(); k++){
				fputs(Theta[j].Values[k], outputFile);
				fputs("\t", outputFile);
			}
			fputs("\n", outputFile);
		}
		fclose(outputFile);
	}else{
		cout << "Error in writing output file." << endl;
		exit(1);
	}
	WriteOutCounter++;
}

void ParameterSetList::Add(ParameterSet Theta_In){
	Theta[CurrentIteration].Initialize(Theta_In);
	CurrentIteration++;
	
	if(CurrentIteration > mcmc->WRITEOUT){
		PrintData();
		CurrentIteration = 1;
		ParameterSet Theta_0 = Theta[0];
		delete[] Theta;
		Theta = new ParameterSet[mcmc->WRITEOUT](this);
		Theta[0] = Theta_0;
	}
}

ParameterSet CurrentParameters(){
	return Theta[CurrentIteration];
}
#endif