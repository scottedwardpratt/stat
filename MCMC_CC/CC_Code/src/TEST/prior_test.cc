#ifndef __PRIOR_CC__
#define __PRIOR_CC__

#include "distribution.h"

using namespace std;

PriorDistribution_Test::PriorDistribution_Test(MCMC * mcmc_in){
	mcmc=mcmc_in;
	SepMap = parameter::getB(mcmc->parmap, "PRIOR_PARAMETER_MAP", false);
	
	if(SepMap){
		string parmapfile = mcmc->dir_name + "/parameters/prior.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
		//parameter::ReadParsFromFile(parmap, parameter_file_name);
	}else{
		parmap = &(mcmc->parmap);
	}
}

double PriorDistribution_Test::Evaluate(ParameterSet Theta){
	double mean = parameter::getD(*parmap, "PRIOR_MEAN", -3.7372);
	double sigma = parameter::getD(*parmap, "PRIOR_SIGMA", 1.6845);
	return Normal(log(Theta.GetValue("SIGMA")), mean, sigma);
}
#endif