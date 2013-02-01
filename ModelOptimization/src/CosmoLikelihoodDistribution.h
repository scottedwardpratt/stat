#ifndef __CosmoLikelihoodDistribution_h__
#define __CosmoLikelihoodDistribution_h__

#include "LikelihoodDistribution.h"

namespace madai {

class CosmoLikelihoodDistribution : public LikelihoodDistribution {
public:
  CosmoLikelihoodDistribution(Model *in_Model);
  ~CosmoLikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);

private:
  std::vector<double> GetData();
  
  std::vector<double> m_Data;
  std::vector<int>    m_intData;
};


} // end namespace madai

#endif // __CosmoLikelihoodDistribution_h__

