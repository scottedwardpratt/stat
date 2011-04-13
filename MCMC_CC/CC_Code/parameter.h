#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include "mcmc.h"

using namespace std;

class ParameterSetList; class MCMC; class ParameterSet;

class ParameterSetList{
public:
	ParameterSetList(MCMC *mcmc_in);
	MCMC *mcmc;
	ParameterSet *Theta;
	int WriteOutCounter;
	vector<bool> LogParam;
	vector<string> ParamNames;
	int CurrentIteration;
	
	void Add(ParameterSet Theta_In);
	void GetTheta0();
	void PrintData();
	ParameterSet CurrentParameters();
};

class ParameterSet{
public:
	vector<string> Names;
	vector<double> Values;
	ParameterSet(ParameterSetList *list);
	void Initialize(vector<string> names, vector<double> values);
	
private:
	bool Used;
	ParameterSetList *paramlist;
};
#endif