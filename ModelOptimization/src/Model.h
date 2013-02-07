/*=========================================================================
 *
 *  Copyright The University of North Carolina at Chapel Hill
 *  All rights reserved.
 *
 *  Licensed under the MADAI Software License. You may obtain a copy of
 *  this license at
 *
 *         https://madai-public.cs.unc.edu/software/license/
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef __Model_h__
#define __Model_h__

#include "Parameter.h"
#include "parametermap.h"
#include "random.h"

#include <gsl/gsl_randist.h>

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
  }

  /** Get the names of the scalar outputs of the model. */
  virtual const std::vector< std::string > & GetScalarOutputNames() const
  {
    return m_ScalarOutputNames;
  }

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

  // Proposed function for interaction with the MCMC:
  virtual ErrorType GetLikeAndPrior( const std::vector< double > & parameters,
                                     double & LikeNew,
                                     double & PriorNew) const = 0;


  std::string   m_DirectoryName;
  std::string   m_ParameterFileName;
  bool          m_LogLike;
  parameterMap  m_ParameterMap;

protected:
  /** Subclasses must populate this vector with the names of the
   * model parameters. */
  std::vector< Parameter > m_Parameters;

  /** Subclasses must populate these vectors with the names of the
   * scalar outputs. */
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

} // end namespace madai

#endif // __Model_h