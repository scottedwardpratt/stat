#ifndef __EMULATOR_H__
#define __EMULATOR_H__

#include "mcmc.h"

using namespace std;

class MCMC;
class ParameterSet;

class EmulatorHandler{
public:
	EmulatorHandler(parameterMap *parmap, MCMCConfiguration * mcmc_in);
	~EmulatorHandler();
	void QueryEmulator(ParameterSet Theta,vector<double> &Means, vector<double> &Errors);
	string Observables;
private:
	string EmulatedParams;
	MCMCConfiguration *mcmc;
	string EmulatorScriptHome;
	string EmInputFile;
	string EmOutputFile;
	string EmErrorFile;
	char * pPath;
};

#endif