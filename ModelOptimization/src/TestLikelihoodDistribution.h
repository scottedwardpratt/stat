#ifndef __TestLikelihoodDistribution_h__
#define __TestLikelihoodDistribution_h__

#include "LikelihoodDistribution.h"

namespace madai {

class TestLikelihoodDistribution : public LikelihoodDistribution {
public:
  TestLikelihoodDistribution(Model *in_Model);
  ~TestLikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);
  
private:
  std::vector<double> GetData();
  
  std::vector<double> m_Data;
};

} // end namespace madai

#endif // __TestLikelihoodDistribution_h__
