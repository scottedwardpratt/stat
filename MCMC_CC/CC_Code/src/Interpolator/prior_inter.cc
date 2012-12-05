#ifndef __PRIOR_CC__
#define __PRIOR_CC__

#include "distribution.h"

using namespace std;

PriorDistribution_Interpolator::PriorDistribution_Interpolator(MCMC * mcmc_in){
	mcmc=mcmc_in;
	SepMap = parameter::getB(mcmc->parmap, "PRIOR_PARAMETER_MAP", false);
	
	if(SepMap){
		string parmapfile = mcmc->dir_name + "/defaultpars/prior.param";
		parmap = new parameterMap;
		parameter::ReadParsFromFile(*parmap, parmapfile);
		//parameter::ReadParsFromFile(*parmap, parameter_file_name);
	}else{
		parmap = &(mcmc->parmap);
	}

	PRIOR = parameter::getS(*parmap,"PRIOR","UNIFORM"); // Specify what type of prior function to use: UNIFORM, GAUSSIAN, STEP. //I need to add "MIXED"
	SCALED = parameter::getB(*parmap,"SCALED",true); //Specifies wethere the values are given 0 to 1, or min_val to max_val

	if( strcmp(PRIOR.c_str(),"GAUSSIAN")==0 ){
		// Read in parameters for Gaussians
		GAUSSIAN_MEANS = parameter::getV(*parmap, "GAUSSIAN_MEANS","");
		GAUSSIAN_STDVS = parameter::getV(*parmap, "GAUSSIAN_STDVS","");
		if(GAUSSIAN_MEANS.size()==0 || GAUSSIAN_STDVS.size()==0){
			cout << "Error in prior_Interpolator.cc: PriorDistribution_Interpolator::PriorDistribution_Interpolator(MCMC)" << endl;
			cout << "GAUSSIAN_MEANS or GAUSSIAN_STDVS not specified. Exiting" << endl;
			exit(1);
		}
		if(GAUSSIAN_MEANS.size() != GAUSSIAN_STDVS.size()){
			cout << "Error in prior_Interpolator.cc: PriorDistribution_Interpolator::PriorDistribution_Interpolator(MCMC)" << endl;
			cout << "Length of GAUSSIAN_MEANS and GAUSSIAN_STDVS are not the same." << endl;
			cout << "Lenght of GAUSSIAN_MEANS = " << GAUSSIAN_MEANS.size() << endl;
			cout << "Lenght of GAUSSIAN_STDVS = " << GAUSSIAN_STDVS.size() << endl;
			cout << "Exiting" << endl;
			exit(1);
		}
	}
}

double PriorDistribution_Interpolator::Evaluate(vector<double> Theta){
	/*double mean = parameter::getD(*parmap, "PRIOR_MEAN", -3.7372);
	double sigma = parameter::getD(*parmap, "PRIOR_SIGMA", 1.6845);
	return Normal(log(Theta.GetValue("SIGMA")), mean, sigma);*/
	if( strcmp(PRIOR.c_str(), "UNIFROM")==0 ){
		// If the prior is uniform, it doesn't matter what value we return as long as it is consistent
		return 1.0;
	}
	if( strcmp(PRIOR.c_str(), "GAUSSIAN")==0 ){
		// The return value needs to be caluclated form a multivariate gaussian
		int N = Theta.size();
		gsl_matrix * sigma = gsl_matrix_calloc(N,N);
		gsl_vector * theta = gsl_vector_alloc(N);
		gsl_vector * means = gsl_vector_alloc(N);
		for(int i = 0; i<N; i++){
			gsl_matrix_set(sigma, i, i,GAUSSIAN_STDVS[i]);
			gsl_vector_set(theta, i, Theta[i]);
			gsl_vector_set(means, i, GAUSSIAN_MEANS[i]);
		}
		return MVNormal(*theta,*means,*sigma);
	}
}
#endif
