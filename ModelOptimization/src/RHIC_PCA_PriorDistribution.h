#ifndef __RHIC_PCA_PriorDistribution_h__
#define __RHIC_PCA_PriorDistribution_h__

#include "PriorDistribution.h"

namespace madai {

class RHIC_PCA_PriorDistribution : public PriorDistribution {
public:
  RHIC_PCA_PriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta);
};


} // end namespace madai

#endif // __RHIC_PCA_PriorDistribution_h__
