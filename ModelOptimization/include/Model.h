/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
	Cory Quammen <cquammen AT cs.unc.edu>
	Russell Taylor <taylorr AT cs.unc.edu>
	Scott Pratt <pratt AT nscl.msu.edu>
	Kevin Novak <novakkev AT msu.edu>
	Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#ifndef __Model_h_
#define __Model_h_

#include "Parameter.h"
#include "random.h"
#include <gsl/gsl_randist.h>
#include "parametermap.h"

#include <cfloat>
#include <vector>

namespace madai {

	class Model {
	public:
		typedef enum {
			NO_ERROR = 0,
			INVALID_OUTPUT_INDEX,
			INVALID_ACTIVE_PARAMETERS,
			FILE_NOT_FOUND_ERROR,
			OTHER_ERROR
		} ErrorType;
		
		Model() {};
		virtual ~Model() {};
		
		/** Loads a configuration from a file. **/
		virtual ErrorType LoadConfigurationFile( const std::string fileName ) = 0;
		
		/** Get the number of parameters. */
		virtual unsigned int GetNumberOfParameters() const
		{
			return static_cast<unsigned int>(m_Parameters.size());
		};
		
		/** Get names of the parameters. */
		virtual const std::vector< Parameter > & GetParameters() const
		{
			return m_Parameters;
		};
		
		/** Set the parameters. */
		virtual ErrorType SetParameters(const std::vector< Parameter > &){
			return NO_ERROR;
		}
		
		/** Get the number of scalar outputs. */
		virtual unsigned int GetNumberOfScalarOutputs() const
		{
			return static_cast<unsigned int>(m_ScalarOutputNames.size());
		};
		
		/** Get the names of the scalar outputs of the model. */
		virtual const std::vector< std::string > & GetScalarOutputNames() const
		{
			return m_ScalarOutputNames;
		};
		
		/** Get the valid range for the parameter at parameterIndex. */
		virtual void GetRange( unsigned int parameterIndex, double range[2] ) const
		{
			range[0] = this->m_Parameters.at(parameterIndex).m_MinimumPossibleValue;
			range[1] = this->m_Parameters.at(parameterIndex).m_MaximumPossibleValue;
		}
		
		/** Get the scalar outputs from the model evaluated at x. */
		virtual ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                                            std::vector< double > & scalars ) const = 0;
		
		/** Get both scalar values and the gradient of the parameters. */
		virtual ErrorType GetScalarAndGradientOutputs( const std::vector< double > & parameters,
                                                      const std::vector< bool > & activeParameters,
                                                      std::vector< double > & scalars,
                                                      unsigned int outputIndex, 
                                                      std::vector< double > & gradient) const = 0;
  
    std::string   m_DirectoryName;
    std::string   m_ParameterFile;
    std::string   m_ParameterFileName;
    bool          m_LogLike;
    bool          m_PrescaledParams;
    parameterMap  m_ParameterMap;
		
	protected:
		/** Subclasses must populate this vector with the names of the
		 model parameters. */
		std::vector< Parameter > m_Parameters;
		
		/** Subclasses must populate these vectors with the names of the
		 scalar outputs. */
		std::vector< std::string > m_ScalarOutputNames;
		
		/** Add a parameter. */
		void AddParameter( const std::string & name,
											double minimumPossibleValue = -DBL_MAX,
											double maximumPossibleValue =  DBL_MAX )
		{
			m_Parameters.push_back(
				Parameter(name, minimumPossibleValue, maximumPossibleValue) );
		}
		
		/** Add a scalar output name. */
		void AddScalarOutputName( const std::string & name )
		{
			m_ScalarOutputNames.push_back( name );
		}
		
	}; // end Model
	
	class CRHICQuadModel : public Model{
	public:
		
	protected:
		
		
		
	};
	
} // end namespace madai

#endif // __Model_h_
