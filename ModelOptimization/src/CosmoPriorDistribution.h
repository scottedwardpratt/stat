#ifndef __CosmoPriorDistribution_h__
#define __CosmoPriorDistribution_h__

#include "PriorDistribution.h"


namespace madai {

class CosmoPriorDistribution : public PriorDistribution {
public:
  CosmoPriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};


} // end namespace madai

#endif // __CosmoPriorDistribution_h__
