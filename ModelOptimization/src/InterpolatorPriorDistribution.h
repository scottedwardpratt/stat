#ifndef __InterpolatorPriorDistribution_h__
#define __InterpolatorPriorDistribution_h__

#include "PriorDistribution.h"


namespace madai {

class InterpolatorPriorDistribution : public PriorDistribution{
public:
  InterpolatorPriorDistribution(Model * in_Model);
  double Evaluate(std::vector<double> Theta);
  
private:

  std::string              m_Prior;
  bool                     m_Scaled;
  std::vector<double>      m_GaussianMeans;
  std::vector<double>      m_GaussianSTDVS;
  std::vector<double>      m_StepMeans;
  std::vector<std::string> m_StepSide;
};

} // end namespace madai

#endif // __InterpolatorPriorDistribution_h__

