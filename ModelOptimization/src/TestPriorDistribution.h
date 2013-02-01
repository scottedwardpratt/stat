#ifndef __TestPriorDistribution_h__
#define __TestPriorDistribution_h__

#include "PriorDistribution.h"

namespace madai {

class TestPriorDistribution : public PriorDistribution {
public:
  TestPriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};

} // end namespace madai

#endif // __TestPriorDistribution_h__
