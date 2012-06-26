/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>

#include "ExternalModel.h"
#include "SimpleMetropolisHastings.h"
#include "Trace.h"

/**
 * Test case for madai::ExternalModel and
 * madai::SimpleMetropolisHastings classes.
 */
int main(int argc, char ** argv) {
	if (argc < 2) {
		std::cerr <<
			"Useage:\n\texternal_model_plus_mhmcmc_test CONFIG_FILE\n\n"
			"where CONFIG_FILE is a suitable configuration for the\n"
			"madai::ExternalModel class.\n\n";
		return EXIT_FAILURE;
	}
	std::string filename(argv[1]);
	madai::ExternalModel external_model;
	if (filename == "-") {
		//FIXME document this.
		external_model.LoadConfigurationFile(std::cin);
	} else {
		external_model.LoadConfigurationFile(filename);
	}
	if (! external_model.good()) {
	  std::cerr << "Something is wrong with the external model\n";
	  return 1;
	}
	madai::SimpleMetropolisHastings simple_mcmc(&external_model);

	std::vector< madai::Parameter > const * parameters
		= &(external_model.GetParameters());
	for (unsigned int i = 0; i < parameters->size(); i++)
		simple_mcmc.ActivateParameter((*parameters)[i].m_Name);
	simple_mcmc.SetOutputScalarToOptimize(
		external_model.GetScalarOutputNames().at(0) );

	madai::Trace trace;
	unsigned int numberIter = 10000;
	for (unsigned int count = 0; count < numberIter; count ++)
		simple_mcmc.NextIteration(&trace);

	trace.writeHead(std::cout, external_model.GetParameters(), external_model.GetScalarOutputNames());
	trace.write(std::cout);
	return 0;
}
