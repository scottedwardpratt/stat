#ifndef __RHIC_LikelihoodDistribution_h__
#define __RHIC_LikelihoodDistribution_h__

#include "LikelihoodDistribution.h"

namespace madai {

class RHIC_LikelihoodDistribution : public LikelihoodDistribution {
public:
  RHIC_LikelihoodDistribution(Model *in_Model);
  ~RHIC_LikelihoodDistribution();
  double Evaluate(std::vector<double> ModelMeans,
                  std::vector<double> ModelErrors);

private:

  std::vector<double>  GetFakeData();
  std::vector<double>  GetRealData();
  std::vector<double>  GetRealError();

  double*              m_DataMean;
  double*              m_DataError;
  std::vector<double>  m_Data;
  std::vector<double>  m_Error;
  bool                 m_FakeData;
  bool                 m_SuppressErrors;
  parameterMap         m_ObservablesParamMap;

  int FindParam(std::string param_name, std::vector<std::string> PNames);
};


} // end namespace madai

#endif // __RHIC_LikelihoodDistribution_h__
