#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include "mcmc.h"

using namespace std;

class ParameterSetList; class MCMC; class ParameterSet;

class ParameterSetList{
public:
	ParameterSetList(MCMC *mcmc_in);
	MCMC *mcmc;
	ParameterSet **Theta;
	int WriteOutCounter;
	vector<bool> LogParam;
	vector<string> ParamNames;
	string EmulatorParams;
	int CurrentIteration;
	
	void Add(ParameterSet Theta_In);
	void GetTheta0();
	void PrintDataToFile();
	ParameterSet CurrentParameters();
};

class ParameterSet{
public:
	vector<string> Names;
	vector<double> Values;
	ParameterSetList *paramlist;
	
	ParameterSet(ParameterSetList *list);
	void Initialize(ParameterSet ParamSetIn);
	void Initialize(vector<string> names, vector<double> values);
	void Print();
	void Reset();
private:
	bool Used;
};
#endif