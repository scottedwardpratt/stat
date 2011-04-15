#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include "mcmc.h"

using namespace std;

class MCMC;

class EmulatorHandler{
public:
	EmulatorHandler(MCMC *mcmc_in);
	~EmulatorHandler();
	void QueryEmulator(ParameterSet Theta,vector<double> &Means, vector<double> &Errors);
private:
	MCMC *mcmc;
	string EmulatorScriptHome;
	string EmInputFile;
	string EmOutputFile;
	char * pPath;
};

#endif