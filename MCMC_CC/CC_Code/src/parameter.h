#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include "mcmc.h"

using namespace std;

class ParameterSetList;
class MCMCConfiguration;
class MCMCRun;
class ParameterSet;

class ParameterSetList{
public:
	ParameterSetList(MCMCRun *mcmc_in);
	ParameterSetList(MCMCRun *mcmc_in, ParameterSet Theta0);
	MCMCRun *mcmc;
	ParameterSet **Theta;
	int WriteOutCounter;
	vector<string> ParamNames;
	string EmulatorParams;
	int CurrentIteration;
	
	void Add(ParameterSet Theta_In);
	void GetTheta0FromFile();
	void PrintDataToFile();
	void WriteOut();
	ParameterSet CurrentParameters();
	ParameterSet * HoldOver;
};

class ParameterSet{
public:
	vector<string> Names;
	vector<double> Values;
	ParameterSetList *paramlist;
	
	ParameterSet(ParameterSetList *list);
	ParameterSet();
	void Initialize(ParameterSet ParamSetIn);
	void Initialize(vector<string> names, vector<double> values);
	void Print();
	void Reset();
	int GetIndex(string ParamName);
	double GetValue(string ParamName);
	void VizTrace();
	bool InTrace;
	bool Used;
};
#endif