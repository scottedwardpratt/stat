#ifndef __PRIOR_CC__
#define __PRIOR_CC__

#include "distribution.h"

using namespace std;

PriorDistribution::PriorDistribution(MCMC * mcmc_in):Distribution(mcmc_in){
	
}

double PriorDistribution::Evaluate(ParameterSet Theta){
	return 1.0;
}
#endif