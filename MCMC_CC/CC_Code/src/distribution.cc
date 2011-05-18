#ifndef __DISTRIBUTION_CC__
#define __DISTRIBUTION_CC__

#include "mcmc.h"

using namespace std;

Distribution::Distribution(MCMCConfiguration *mcmc_in){
	mcmc = mcmc_in;
}

double Distribution::Normal(double x, double mu, double sigma){
	return (1/sqrt(2*M_PI*pow(sigma,2)))*exp(-pow((x-mu),2)/(2*pow(sigma,2)));
}

double Distribution::Gaussian(double x, double mu, double sigma){
	return Normal(x,mu,sigma);
}

double Distribution::Log_MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma){
	double OUT, det;
	int foobar,N;
	N = x.size;
	
	gsl_matrix * sigma_inv = gsl_matrix_calloc(N,N);
	gsl_matrix * tempsigma = gsl_matrix_alloc(N,N);
	gsl_permutation * p = gsl_permutation_alloc(N);
	gsl_vector * diff = gsl_vector_alloc(N);
	gsl_vector * temp = gsl_vector_alloc(N);
	
	gsl_vector_memcpy(diff, &x);
	gsl_matrix_memcpy(tempsigma, &sigma);
	gsl_vector_sub(diff, &mu);
	//invert matrix using LU decomposition
	gsl_linalg_LU_decomp(tempsigma, p, &foobar);
	gsl_linalg_LU_invert(tempsigma, p, sigma_inv);
	
	det = gsl_linalg_LU_det(tempsigma, foobar);
	
	//multiply matrix and left vector together using CBLAS routines
	gsl_blas_dsymv(CblasUpper,1.0, sigma_inv, diff, 0.0, temp);
	
	gsl_blas_ddot(diff, temp, &OUT);
	
	OUT = (-1.0/2.0)*OUT -(N/2)*log(2*M_PI) - (0.5)*log(det);
	
	gsl_matrix_free(sigma_inv);
	gsl_permutation_free(p);
	gsl_vector_free(diff);
	gsl_vector_free(temp);
	return OUT;
}

double Distribution::MVNormal(gsl_vector x, gsl_vector mu, gsl_matrix sigma){
	double like = Log_MVNormal(x,mu,sigma);
	return exp(Log_MVNormal(x,mu,sigma));
}

double Distribution::LogNormal(double x, double mu, double sigma){
	double OUT = (1.0/(x*sigma))*(1.0/sqrt(2.0*M_PI))*exp((-1.0/(2.0*pow(sigma,2)))*(pow((log(x)-mu), 2)));
	return OUT;
}

#endif