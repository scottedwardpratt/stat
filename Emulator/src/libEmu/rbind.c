#include "emulator.h"
#include "estimator.h"
#include "maximise.h"
#include "../multifit.h"
#include "../useful.h"

/**
 * just enough setup and teardown to call the emulator directly from R
 */



//! callEmulator from R
/**
 * so R can't handle passing anything other than arrays, 
 * as such anything which is canonically a 2d array (such as xmodel) will be passed 
 * in and out as a flat vector.
 */
void callEmulator(double* xmodel, int* nparams_in,  double* training, int *nmodelpts, int* nthetas_in, double* final_emulated_x, int *nemupts_in, double* final_emulated_y, double* final_emulated_variance , double* range_min_in, double* range_max_in ){
	int i, j;
	int nmodel_points = *nmodelpts;
	int nparams = *nparams_in;
	int nemupts = *nemupts_in;
	int nthetas = *nthetas_in;
	gsl_matrix *xmodel = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_matrix *new_x = gsl_matrix_alloc(nmodel_points, nparams);
	gsl_vector *training_vec = gsl_vector_alloc(nmodel_points);
	gsl_vector *emulated_variance = gsl_vector_alloc(nemupts);
	gsl_vector *emulated_y = gsl_vector_alloc(nemupts);
	gsl_vector *thetas = gsl_vector_alloc(nthetas);
	//\todo sethis up
	gsl_rng* random; 
	const gsl_rng_type *rng_type;

	// need these to call the routines which actually do the emulation
	emuResult theResult;
	eopts theOptions;

	// setup the rng
	T = gsl_rng_default;
	random_number = gsl_rng_alloc(T);
	gsl_rng_set(random, get_seed());

	// fill in xmodel 
	// you have to do this to new_x at the end
	for(j=0; j < nparams; j++){
		for(i = 0; i < nmodelpts; i++){
			gsl_matrix_set(xmodel, i, j, xmodel[i+j*nmodelpts]);
		}
	}
	// fill in the training vec
	for(i = 0; i < nmodelpts; i++){
		gsl_vector_set(training_vec, training[i]);
 

	// fill in the options
	theOptions.nmodel_points = nmodel_points;
	theOptions.nemu_points = nemupts;
	theOptions.nparams = nparams;
	theOptions.range_min = *range_min_in;
	theOptions.range_max = *range_max_in;
	theOptions.xmodel = xmodel;
	theOptions.training = training_vec;
	theOptions.thetas = thetas;

	// alloc the result bits
	theResult.nemu_points = nemupts;
	theResult.nparams = nparams;
	theResult.new_x = new_x;
	theResult.new_mean = emulated_y; 
	theResult.new_var = emualted_var;


	// actually do the emulation, this function 
	// can be found in ../multifit.c
	evaluateRegion(&theResult, &theOptions, random);
	

	// fill in emulated_y, emulated_variance
	for(i = 0; i < nemupts; i++){
		final_emulated_y[i] = gsl_vector_get(emulated_y, i);
		final_emulated_variance[i] = gsl_vector_get(emulated_variance, i);
	}
	// fill in final emulated_x
	for(j = 0; j < nparams; j++){
		for(i = 0; i < nmodelpts; i++){
			final_emulated_x[i+j*nmodelpts] = gsl_matrix_get(new_x, i, j);
		}
	}

	// tidy up
	gsl_matrix_free(xmodel);
	gsl_matrix_free(new_x);
	gsl_vector_free(emulated_variance);
	gsl_vector_free(emulated_y);
	gsl_vector_free(thetas);
	gsl_rng_free(random);
}
	
