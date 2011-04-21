#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

#include "mcmc.h"

using namespace std;

class VizHandler{
public:
	VizHandler(MCMC * mcmc_in);
	~VizHandler();
	void operator() (const string& command);
	void UpdateTraceFig();
	string header;
	
protected:
	int HighestItnReadIn;
	MCMC *mcmc;
	string *paramvalues;
	FILE * gnuplotpipe;
};
#endif