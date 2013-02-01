#ifndef __RHIC_PriorDistribution_h__
#define __RHIC_PriorDistribution_h__

#include "PriorDistribution.h"

namespace madai {

class RHIC_PriorDistribution : public PriorDistribution {
public:
  RHIC_PriorDistribution(Model *in_Model);
  double Evaluate(std::vector<double> Theta); //ParameterSet Theta

private:
  std::string              m_Prior;
  bool                     m_Scaled;
  std::vector<double>      m_GaussianMeans;
  std::vector<double>      m_GaussianSTDVS;
  std::vector<double>      m_StepMeans;
  std::vector<std::string> m_StepSide;
};


} // end namespace madai

#endif // __RHIC_PriorDistribution_h__
