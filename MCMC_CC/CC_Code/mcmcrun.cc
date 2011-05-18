#include "mcmc.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: mcmc run_file\n");
		exit(-1);
	}
	
	string run_file = argv[1];
	MCMCConfiguration *mcmcfoobar = new MCMCConfiguration(run_file);
	cout << "Done making config." << endl;
	MCMCRun *mcmcrun = new MCMCRun(mcmcfoobar);
	cout << "Done making run" << endl;
	
	return 0;
}