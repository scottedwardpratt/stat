#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: mcmc run_file\n");
		exit(-1);
	}
	/*
	double HYDRO_T0 0.6 1.1
	double GLAUBER_WNBIN_RATIO 0.6 1.0  
	double GLAUBER_K_TAU 2.4 3.3
	double HYDRO_INIT_NS 0.5 2.0
	double HYDRO_INIT_FLOW 0.4 1.2
	double HYDRO_SVRATIO 0.03 0.25
	double EQOFST_BUMP_HEIGHT 0.2 0.8
	*/

	string run_file = argv[1];
	MCMCConfiguration *mcmclowlambda = new MCMCConfiguration(run_file, "lowlambda");
	MCMCConfiguration *mcmcbiglambda = new MCMCConfiguration(run_file, "biglambda");
	vector<string> Names = mcmcbiglambda->ParamNames;
	MCMCRun* runlow;
	MCMCRun* runbig;
	ifstream input;
	ParameterSet * Theta0 = new ParameterSet();
	
	input.open("params.dat");
	
	if(input){
		for(int i = 0; i < 5; i++){
			Theta0->Reset();
			vector<double> Tempvalues;
			for(int index = 0; index< 8; index++){
				double temp;
				input >> temp;
				Tempvalues.push_back(temp);
			}
			Theta0->Initialize(Names, Tempvalues);
			Theta0->Print();
			runbig = new MCMCRun(mcmcbiglambda, *Theta0);
			runbig->Run();
		}
	}
	else{
		cout << "Unable to open parameter file." << endl;
	}
}