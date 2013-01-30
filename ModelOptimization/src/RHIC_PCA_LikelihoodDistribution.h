#ifndef __RHIC_PCA_LikelihoodDistribution_h__
#define __RHIC_PCA_LikelihoodDistribution_h__


#include "LikelihoodDistribution.h"


namespace madai {

class RHIC_PCA_LikelihoodDistribution : public LikelihoodDistribution{
public:
  RHIC_PCA_LikelihoodDistribution(Model *in_Model);
  ~RHIC_PCA_LikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);
    
private:
  std::vector<double> GetRealData(); 
  
  double*             m_DataMean;
  double*             m_DataError;
  std::vector<double> m_Data;
  std::vector<double> m_Error;
  parameterMap        m_ObservablesParamMap;
};

} // end namespace madai

#endif // __RHIC_PCA_LikelihoodDistribution_h__

