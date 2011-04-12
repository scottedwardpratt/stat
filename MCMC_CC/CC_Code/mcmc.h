#ifndef __MCMC_H__
#define __MCMC_H__

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <iostream>
#include <fstream>
#include "parameter.h"
#include "distribution.h"

using namespace std;

//Hacked out of parametermap.h from Coralutils folder. Fix later.
typedef map<string,string> parameterMap;

namespace parameter {
  bool   getB(parameterMap ,string ,bool);
  int    getI(parameterMap ,string ,int);
  string getS(parameterMap ,string ,string);
  double getD(parameterMap ,string ,double);
  vector< double > getV(parameterMap, string, double);
  vector< string > getVS(parameterMap, string, string);
  vector< vector< double > > getM(parameterMap, string, double);
  void set(parameterMap&, string, double);
  void set(parameterMap&, string, int);
  void set(parameterMap&, string, bool);
  void set(parameterMap&, string, string);
  void set(parameterMap&, string, char*);
  void set(parameterMap&, string, vector< double >);
  void set(parameterMap&, string, vector< string >);
  void set(parameterMap&, string, vector< vector< double > >);
  void ReadParsFromFile(parameterMap&, const char *filename);
  void ReadParsFromFile(parameterMap&, string filename);
  void PrintPars(parameterMap&);
};

class ParameterSet; class ParameterSetList; class ProposalDistribution;

class MCMC{
public:
	parameterMap parmap;
	ParameterSetList *ThetaList;
	MCMC(string run_file);
	~MCMC();
	void Run();
	
private:
	int MAXITERATIONS, WRITEOUT, Accept_Count;
	string dir_name, parameter_file_name;
	
	bool LOGLIKE;
	bool LOGPRIOR;
	bool LOGPROPOSAL;
	

	// LikelihoodDistribution *Likelihood;
	//ProposalDistribution *Proposal;
	// PriorDistribution *Prior;
};

#endif