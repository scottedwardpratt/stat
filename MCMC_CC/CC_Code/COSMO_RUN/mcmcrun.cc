#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: mcmc run_file\n");
		exit(-1);
	}
	
	string run_file = argv[1];
	MCMC *mcmcdefault = new MCMC(run_file, "default");
	vector<string> Names = mcmcdefault->ParamNames;
	MCMC *run;
	ifstream input;
	// ParameterSet * Theta0 = new ParameterSet();
	vector<double> Ratios;
	double ratio;
			
	run = new MCMC(mcmcdefault);
	ratio = run->Run();

	cout << "Done successfully." << endl;
}