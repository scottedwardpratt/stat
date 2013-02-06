#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include "mcmc.h"

using namespace std;

class MCMC;

class EmulatorHandler{
public:
	EmulatorHandler(parameterMap *parmap, MCMC * mcmc_in);
	~EmulatorHandler();
	void QueryEmulator(vector<double> Theta,vector<double> &Means, vector<double> &Errors);
private:
	string EmulatedParams;
	MCMC *mcmc;
	string EmulatorScriptHome;
	string EmInputFile;
	string EmOutputFile;
	string EmErrorFile;
	char * pPath;
};

#endif