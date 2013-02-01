#ifndef __PriorDistribution_h__
#define __PriorDistribution_h__

#include "Distribution.h"


namespace madai {

class PriorDistribution : public Distribution {
public:
  PriorDistribution();
  virtual ~PriorDistribution();
  virtual double Evaluate(std::vector<double> Theta);
};

} // end namespace madai


#endif // __PriorDistribution_h__

