#ifndef __PRIOR_CC__
#define __PRIOR_CC__

#include "distribution.h"

using namespace std;

PriorDistribution_RHIC::PriorDistribution_RHIC(MCMCConfiguration * mcmc_in){
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
			cout << "Error in prior_rhic.cc: PriorDistribution_RHIC::PriorDistribution_RHIC(MCMCConfiguration)" << endl;
			cout << "GAUSSIAN_MEANS or GAUSSIAN_STDVS not specified. Exiting" << endl;
			exit(1);
		}
		if(GAUSSIAN_MEANS.size() != GAUSSIAN_STDVS.size()){
			cout << "Error in prior_rhic.cc: PriorDistribution_RHIC::PriorDistribution_RHIC(MCMCConfiguration)" << endl;
			cout << "Length of GAUSSIAN_MEANS and GAUSSIAN_STDVS are not the same." << endl;
			cout << "Lenght of GAUSSIAN_MEANS = " << GAUSSIAN_MEANS.size() << endl;
			cout << "Lenght of GAUSSIAN_STDVS = " << GAUSSIAN_STDVS.size() << endl;
			cout << "Exiting" << endl;
			exit(1);
		}
	}
	if( strcmp(PRIOR.c_str(),"STEP")==0 ){
		/* JFN 9/24/12 9:20am- Just a note. At the moment I am coding a step function prior because Scott asked me to. But this is
		both a) going to cause runtime errors, and b) the wrong way to do this. If we have a step function prior, that should be
		actualized by adjusting the ranges.*/
		// Read in parameters for step functions
		STEP_MEANS = parameter::getV(*parmap, "STEP_MEANS", "");
		STEP_SIDE  = parameter::getVS(*parmap, "STEP_SIDE", "");
		if(STEP_MEANS.size()==0 || STEP_SIDE.size()==0){
			cout << "Error in prior_rhic.cc: PriorDistribution_RHIC::PriorDistribution_RHIC(MCMCConfiguration)" << endl;
			cout << "STEP_MEANS or STEP_SIDE not specified. Exiting" << endl;
			exit(1);
		}
		if(STEP_MEANS.size() != STEP_SIDE.size()){
			cout << "Error in prior_rhic.cc: PriorDistribution_RHIC::PriorDistribution_RHIC(MCMCConfiguration)" << endl;
			cout << "Length of STEP_MEANS and STEP_SIDE are not the same." << endl;
			cout << "Lenght of STEP_MEANS = " << STEP_MEANS.size() << endl;
			cout << "Lenght of STEP_SIDE = " << STEP_SIDE.size() << endl;
			cout << "Exiting" << endl;
			exit(1);
		}
	}
}

double PriorDistribution_RHIC::Evaluate(ParameterSet Theta){
	/*double mean = parameter::getD(*parmap, "PRIOR_MEAN", -3.7372);
	double sigma = parameter::getD(*parmap, "PRIOR_SIGMA", 1.6845);
	return Normal(log(Theta.GetValue("SIGMA")), mean, sigma);*/
	if( strcmp(PRIOR.c_str(), "UNIFROM")==0 ){
		// If the prior is uniform, it doesn't matter what value we return as long as it is consistent
		return 1.0;
	}
	if( strcmp(PRIOR.c_str(), "GAUSSIAN")==0 ){
		// The return value needs to be caluclated form a multivariate gaussian
		int N = Theta.Values.size();
		gsl_matrix * sigma = gsl_matrix_calloc(N,N);
		gsl_vector * theta = gsl_vector_alloc(N);
		gsl_vector * means = gsl_vector_alloc(N);
		for(int i = 0; i<N; i++){
			gsl_matrix_set(sigma, i, i,GAUSSIAN_STDVS[i]);
			gsl_vector_set(theta, i, Theta.Values[i]);
			gsl_vector_set(means, i, GAUSSIAN_MEANS[i]);
		}
		return MVNormal(*theta,*means,*sigma);
	}
	if( strcmp(PRIOR.c_str(), "STEP")==0 ){
		// The return value is either 0, or an arbitary consistent value
		for(int i=0; i < Theta.Values.size(); i++){
			// At the moment, scaled vs unscaled doesn't make any difference. In both they are treated as scaled. To be fixed
			if(SCALED){
				if(((Theta.Values[i] > STEP_MEANS[i]) && ( strcmp(STEP_SIDE[i].c_str(),"LOW"))) || ((Theta.Values[i] < STEP_MEANS[i]) && ( strcmp(STEP_SIDE[i].c_str(),"HIGH")))){
					return 0; //... this seems like it could cause problems
				}
			}
			if(!SCALED){
				if(((Theta.Values[i] > STEP_MEANS[i]) && ( strcmp(STEP_SIDE[i].c_str(),"LOW"))) || ((Theta.Values[i] < STEP_MEANS[i]) && ( strcmp(STEP_SIDE[i].c_str(),"HIGH")))){
					return 0; //... this seems like it could cause problems
				}
			}
		}
		return 1.0; // if the thetas have survived all of the step function checks
	}
}
#endif
