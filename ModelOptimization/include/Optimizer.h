#ifndef __Optimizer_h_
#define __Optimizer_h_

#include <set>

#include "Model.h"
#include "Trace.h"

namespace madai {

class Optimizer {
public:
  typedef enum {
    NO_ERROR = 0,
    INVALID_PARAMETER_INDEX_ERROR,
    INVALID_OUTPUT_SCALAR_INDEX_ERROR
  } ErrorType;

  Optimizer( const Model *model );
  virtual ~Optimizer();

  void ActivateParameter( const std::string & parameterName );
  void DeactivateParameter( const std::string & parameterName );

  /** Resets a parameter value. */
  virtual ErrorType SetParameterValue( const std::string & parameterName,
                                  double value );

  /** Sets the output scalar value to optimize. */
  ErrorType SetOutputScalarToOptimize( const std::string & scalarName );

  /** Compute the next set of parameters and the output scalar values,
   * and save them in the trace file. */
  virtual void NextIteration(Trace *trace) = 0;
  // scalars;
  // m_Model->GetScalarOutputs( currentPosition, scalars, gradient );
  // m_Trace->RecordData( scalars )

  // update position

  /** Get the current parameter values. */
  const std::vector< double > & GetCurrentParameters() const;

protected:
  const Model *m_Model;

  std::set< std::string > m_ActiveParameters;

  std::vector< double > m_CurrentParameters;

  std::string m_OutputScalarToOptimize;

  Optimizer() {}; // intentionally hidden

  /** Subclasses that need to reset internal state when a parameter
  value has been changed outside the operation of the optimization
  algorithm should override this method. */
  virtual void ParameterSetExternally() {};

  unsigned int GetParameterIndex( const std::string & parameterName );

}; // end Optimizer

} // end namespace madai

#endif // __Optimizer_h_
