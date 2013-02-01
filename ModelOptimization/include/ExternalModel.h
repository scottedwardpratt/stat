/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#ifndef __ExternalModel_h_
#define __ExternalModel_h_

#include "Parameter.h" // Parameter class 
#include "Model.h"

#include "process_pipe.h"
/* madai::process_pipe struct, madai::create_process_pipe fn.
 Note that process_pipe is full of POSIX-specific inter-process
 communication code.  We can and should fix that at some point.  Until we
 fix that, this class will only compile for a POSIX (Unix, Linux,
 BSD, Macosx) target.  */

#include <vector> // std::vector
#include <string> // std::string 
#include <iostream>


namespace madai {

	class ExternalModel : public Model {
	private:
		unsigned int number_of_parameters, number_of_outputs;
		std::vector< std::string > command_arguments;
		process_pipe process;
		std::string configFileName;
		typedef enum {
			UNINITIALIZED,
			READY,
			ERROR
		} internal_state;
		internal_state stateFlag;
		
	public:

		bool good() { return (this->stateFlag == READY); }
		ExternalModel();
		ExternalModel(const std::string & configFileName);
		virtual ~ExternalModel();

		/** 
		 * Loads a configuration from a file.  The format of the file is
		 * defined by this function.  We'll lock it down later.
		 */
		virtual ErrorType LoadConfigurationFile( const std::string fileName );
		virtual ErrorType LoadConfigurationFile( std::istream & configFile );

		/** Get the valid range for the parameter at parameterIndex. */
		virtual void GetRange(unsigned int parameterIndex, double range[2]) const;

		/** 
		 * Get the scalar outputs from the model evaluated at x.  If an
		 * error happens, the scalar output array will be left incomplete.
		 */
		virtual ErrorType GetScalarOutputs(
				const std::vector< double > & parameters,
				std::vector< double > & scalars ) const;

	  // Not implemented yet.  Should we do this numerically?
		/** Get both scalar values and the gradient of the parameters. */
		virtual ErrorType GetScalarAndGradientOutputs(
			const std::vector< double > & parameters,
			const std::vector< bool > & activeParameters,
			std::vector< double > & scalars,
      unsigned int outputIndex, std::vector< double > & gradient) const;
    
    // Not implemented yet.
    // Proposed function for interaction with the MCMC
    virtual ErrorType GetLikeAndPrior( const std::vector< double > & parameters,
                                       double & Like,
                                       double & Prior ) const;
	}; // end Model
} // end namespace madai

#endif // __Model_h_
