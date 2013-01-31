/*********************************************************************
MADAI Model Statistical Tools
Copyright 2011-2012, The University of North Carolina at Chapel Hill.

This software was written in 2011-2012 by 
  Hal Canary <hal AT cs.unc.edu>
while working for the MADAI project <http://madai.us/>.

See copyright.txt for more information.
*********************************************************************/
#ifndef __SimpleMetropolisHastings_h__
#define __SimpleMetropolisHastings_h__

#include "Optimizer.h"

namespace madai {

class SimpleMetropolisHastings : public Optimizer {
public:
  SimpleMetropolisHastings( const Model *model );
  ~SimpleMetropolisHastings();

  void NextIteration(Trace *trace);

  /** Set the step size. */
  void SetStepSize( double stepSize );

protected:
  double m_StepSize;

  SimpleMetropolisHastings() {}; // intentionally hidden

  std::vector< bool > m_ActiveParameters;

  unsigned int m_NumberOfParameters;
  
  unsigned int m_NumberOfOutputs;

}; // end class SimpleMetropolisHastings

} // end namespace madai

#endif // __SimpleMetropolisHastings_h__
