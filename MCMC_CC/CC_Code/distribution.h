#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__ 

using namespace std;

class Distribution{
public:
	double Evaluate(ParameterSet Theta);

};

class ProposalDistribution: Distribution {
public:
	MCMC * mcmc;
	vector<double> MixingStdDev;
	double *Ranges[2];
	
	ProposalDistribution(MCMC * mcmc_in);
	ParameterSet Iterate();
private:
	int FindParam(string param_name);
};
#endif