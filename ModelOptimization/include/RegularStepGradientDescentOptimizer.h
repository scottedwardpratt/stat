#ifndef __RegularStepGradientDescentOptimizer_h_
#define __RegularStepGradientDescentOptimizer_h_

#include "Optimizer.h"

namespace madai {

class RegularStepGradientDescentOptimizer : public Optimizer {
public:
  RegularStepGradientDescentOptimizer( const Model *model );
  ~RegularStepGradientDescentOptimizer();

  void NextIteration(Trace *trace);
  void NextIteration();

  /** Set the step size. */
  void SetStepSize( double stepSize );

protected:
  double m_StepSize;

  RegularStepGradientDescentOptimizer() {}; // intentionally hidden
}; // end class RegularStepGradientDescentOptimizer

} // end namespace madai

#endif // __RegularStepGradientDescentOptimizer_h_
