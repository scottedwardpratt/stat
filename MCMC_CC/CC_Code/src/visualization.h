#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

#include "mcmc.h"
#include <deque>

using namespace std;

class VizHandler{
public:
	VizHandler(MCMCRun * mcmc_in);
	~VizHandler();
	void operator() (const string& command);
	void UpdateTraceFig();
	void FinalTrace();
	string header;
private:
	string gnuplotterm;
	string gnuplotstyle;
	int ThetaListSize;
	int HighestItnReadIn;
	bool MovingWindow;
	MCMCRun *mcmc;
	string *paramvalues;
	deque<string> *DequeParameterValues;
	FILE * gnuplotpipe;
};
#endif