#ifndef __QUAD_H__
#define __QUAD_H__

#include "mcmc.h"

using namespace std;

class MCMC;
class ParameterSet;

class QuadHandler{
public:
	QuadHandler(parameterMap *parmap, MCMCConfiguration * mcmc_in);
	~QuadHandler();
	void QueryQuad(ParameterSet Theta,vector<double> &Means, vector<double> &Errors);
private:
	string EmulatedParams;
	MCMCConfiguration *mcmc;
	string QuadScriptHome;
	string EmInputFile;
	string EmOutputFile;
	string EmErrorFile;
	char * pPath;
};

#endif