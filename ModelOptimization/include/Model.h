#ifndef __Model_h_
#define __Model_h_

#include "Parameter.h"

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

  /** Loads a configuration from a file. */
  virtual ErrorType LoadConfigurationFile( const std::string fileName ) = 0;

  /** Get the number of parameters. */
  virtual unsigned int GetNumberOfParameters() const
  {
    return m_Parameters.size();
  };

  /** Get names of the parameters. */
  virtual const std::vector< Parameter > & GetParameters() const
  {
    return m_Parameters;
  };

  /** Get the number of scalar outputs. */
  virtual unsigned int GetNumberOfScalarOutputs() const
  {
    return m_ScalarOutputNames.size();
  };

  /** Get the names of the scalar outputs of the model. */
  virtual const std::vector< std::string > & GetScalarOutputNames() const
  {
    return m_ScalarOutputNames;
  };

  /** Get the valid range for the parameter at parameterIndex. */
  virtual void GetRange( unsigned int parameterIndex, double range[2] ) = 0;

  /** Get the scalar outputs from the model evaluated at x. */
  virtual ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                                      std::vector< double > & scalars ) const = 0;

  /** Get both scalar values and the gradient of the parameters. */
  virtual ErrorType GetScalarAndGradientOutputs( const std::vector< double > & parameters,
                                                 const std::vector< bool > & activeParameters,
                                                 std::vector< double > & scalars,
                                                 unsigned int outputIndex, std::vector< double > & gradient) const = 0;

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
    Parameter newParameter;
    newParameter.m_Name = name;
    newParameter.m_MinimumPossibleValue = minimumPossibleValue;
    newParameter.m_MaximumPossibleValue = maximumPossibleValue;

    m_Parameters.push_back( newParameter );
  }

  /** Add a scalar output name. */
  void AddScalarOutputName( const std::string & name )
  {
    m_ScalarOutputNames.push_back( name );
  }

}; // end Model

} // end namespace madai

#endif // __Model_h_
