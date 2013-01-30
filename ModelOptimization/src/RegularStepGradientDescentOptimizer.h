#ifndef __RegularStepGradientDescentOptimizer_h__
#define __RegularStepGradientDescentOptimizer_h__


#include "Optimizer.h"


namespace madai {

/** \class RegularStepGradientDescentOptimizer
 *
 * Straightforward implementation of a gradient descent optimizer.
 */
class RegularStepGradientDescentOptimizer : public Optimizer {
public:
  RegularStepGradientDescentOptimizer( const Model *model );
  ~RegularStepGradientDescentOptimizer();

  void NextIteration(Trace *trace);

  /** Set the step size. */
  void SetStepSize( double stepSize );

protected:
  double m_StepSize;

  RegularStepGradientDescentOptimizer() {}; // intentionally hidden
}; // end class RegularStepGradientDescentOptimizer

} // end namespace madai

#endif // __RegularStepGradientDescentOptimizer_h__
