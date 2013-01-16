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
#ifndef __Optimizer_h_
#define __Optimizer_h_

#include <set>

#include "Model.h"
#include "Trace.h"

namespace madai {

class Model;

class Optimizer {
public:
  typedef enum {
    NO_ERROR = 0,
    INVALID_PARAMETER_INDEX_ERROR,
    INVALID_OUTPUT_SCALAR_INDEX_ERROR
  } ErrorType;

  Optimizer( const Model *model );
  virtual ~Optimizer();
  const Model * GetModel() const;
  std::set< std::string > GetActiveParameters();

  void ActivateParameter( const std::string & parameterName );
  void DeactivateParameter( const std::string & parameterName );

  /** Resets a parameter value. */
  virtual ErrorType SetParameterValue( const std::string & parameterName,
                                       double value );

  /** Sets the output scalar value to optimize. */
  ErrorType SetOutputScalarToOptimize( const std::string & scalarName );
  std::string GetOutputScalarToOptimize();

	ErrorType SetOutputScalarToOptimizeIndex(unsigned int idx);
	unsigned int GetOutputScalarToOptimizeIndex() const;

  /** Compute the next set of parameters and the output scalar values,
   * and save them in the trace file. */
  
  virtual void NextIteration(Trace *trace) = 0;
  //{  /* suggested structure for this function */
  //std::vector< double > scalarOutputs;
  //std::vector< double > gradient;
  //int err;
  //err = m_Model->GetScalarAndGradientOutputs(
  // m_CurrentParameters,
  // m_ActiveParameters,
  // scalarOutputs,
  // m_OutputScalarToOptimizeIndex,
  // gradient);
  //    if (err) {
  //      // handle the error
  //    }
  //m_Trace->RecordData(m_CurrentParameters, scalarOutputs);
  //
  // // Based on:
  // //    scalarOutputs[m_OutputScalarToOptimizeIndex]
  // //    m_ActiveParameters
  // //    m_Trace
  // //    gradient
  // //    m_CurrentParameters
  // // Then we need to update
  // //    m_CurrentParameters
  //}

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
	unsigned int m_OutputScalarToOptimizeIndex;

  Optimizer() {}; // intentionally hidden

  /** Subclasses that need to reset internal state when a parameter
  value has been changed outside the operation of the optimization
  algorithm should override this method. */
  virtual void ParameterSetExternally() {};

  unsigned int GetOutputScalarIndex( const std::string & scalarName ) const;
  unsigned int GetParameterIndex( const std::string & parameterName ) const;
  
  bool IsLikeAndPrior() const;

}; // end Optimizer

} // end namespace madai

#endif // __Optimizer_h_
