#ifndef __PRIOR_CC__
#define __PRIOR_CC__

#include "distribution.h"

using namespace std;

PriorDistribution_Cosmo::PriorDistribution_Cosmo(MCMC * mcmc_in){
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

double PriorDistribution_Cosmo::Evaluate(ParameterSet Theta){
	return 1.0;
}
#endif