#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: mcmc run_file\n");
		exit(-1);
	}
	
	string run_file = argv[1];
	MCMC *mcmcdefault = new MCMC(run_file, "default");
	// cout << "Config done." << endl;
	vector<string> Names = mcmcdefault->ParamNames;
	// MCMC* runlow;
	// MCMC* runbig;
	MCMC *run;
	ifstream input;
	ParameterSet * Theta0 = new ParameterSet();
	// vector<double> BigRatios;
	// vector<double> LowRatios;
	vector<double> Ratios;
	double ratio;
	
	// input.open("params.dat");
	input.open("params_nosigma.dat");
	
	if(input){
		for(int i = 0; i < 5; i++){
			Theta0->Reset();
			vector<double> Tempvalues;
			for(int index = 0; index< Names.size(); index++){
				double temp;
				input >> temp;
				Tempvalues.push_back(temp);
			}
			Theta0->Initialize(Names, Tempvalues);
			Theta0->Print();
			run = new MCMC(mcmcdefault, *Theta0);
			ratio = run->Run();
			Ratios.push_back(ratio);
		}
	}
	else{
		cout << "Unable to open parameter file." << endl;
	}
	
	cout << "Acceptance ratios: " << endl;
	
	for(int j = 0; j < Ratios.size(); j++){
		cout << Ratios[j] << endl;
	}
}