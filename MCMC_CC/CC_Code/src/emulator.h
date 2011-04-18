#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include "mcmc.h"

using namespace std;

class MCMC;
class ParameterSet;

class EmulatorHandler{
public:
	EmulatorHandler(parameterMap *parmap, MCMC * mcmc_in);
	~EmulatorHandler();
	void QueryEmulator(ParameterSet Theta,vector<double> &Means, vector<double> &Errors);
private:
	string EmulatedParams;
	MCMC *mcmc;
	string EmulatorScriptHome;
	string EmInputFile;
	string EmOutputFile;
	char * pPath;
};

#endif