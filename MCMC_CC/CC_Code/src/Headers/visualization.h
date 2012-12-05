#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

#include "mcmc.h"
#include <deque>

using namespace std;

class VizHandler{
public:
	VizHandler(MCMC * mcmc_in);
	~VizHandler();
	void operator() (const string& command);
	void UpdateTraceFig();
	void FinalTrace();
private:
	string gnuplotterm;
	string gnuplotstyle;
	int ThetaListSize;
	int HighestItnReadIn;
	bool MovingWindow;
	bool DensityPlot;
	MCMC *mcmc;
	string *paramvalues;
	deque<string> *DequeParameterValues;
	FILE * gnuplotpipe;
	//FILE * gnuplotmultipipe;
	vector<string> DensityPlotFileNames;
	vector<string> DensityPlotCommands;
	string header;
	int Densities[7][7][100][100]; //This is the array which is holding the density plot data. In principle the first two indices should be NUM_PARAMS, and NUM_PARAMS-1, but I don't know how to do that here. The last two are the number of bins.
	int BINS;
};
#endif