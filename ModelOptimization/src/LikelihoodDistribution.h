#ifndef __LikelihoodDistribution_h__
#define __LikelihoodDistribution_h__

#include "Distribution.h"

namespace madai {

class LikelihoodDistribution : public Distribution {
public:
  LikelihoodDistribution();
  virtual ~LikelihoodDistribution();
  virtual double Evaluate(std::vector<double> ModelMeans,
                          std::vector<double> ModelErrors);
  
protected:
  virtual std::vector<double> GetData();

  std::vector<double> m_Data;
  bool                m_UseEmulator;
  bool                m_ProcessPipe;
  std::ofstream       m_EmulatorTest;
};


} // end namespace madai


#endif // __LikelihoodDistribution_h__
