#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: mcmc run_file\n");
		exit(-1);
	}
	
	string run_file = argv[1];
	MCMC *mcmc = new MCMC(run_file, "default");
	mcmc->FirstPass();
	mcmc->Run();

	cout << "Done successfully." << endl;
}