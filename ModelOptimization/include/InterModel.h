#ifndef __InterModel_h__
#define __InterModel_h__

#include "MultiModel.h"
#include "rhicstat.h"

namespace madai {

class InterModel : public MultiModel {
private:
  bool       m_Verbose;
  bool       m_Timing;
  CRHICStat* m_Emulator;
public:
  InterModel();
  InterModel(std::string info_dir);
  virtual ~InterModel();
  
  virtual ErrorType GetScalarOutputs(const std::vector< double > & parameters,
                                     std::vector< double > & scalars ) const;
                     
  // Not implemented yet
  virtual ErrorType GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                                const std::vector< bool > & activeParameters,
                                                std::vector< double > & scalars,
                                                unsigned int outputIndex,
                                                std::vector< double > & gradient) const;
  
  // For interaction with the mcmc
  virtual ErrorType GetLikeAndPrior(const std::vector< double > & parameters,
                                    double & Like,
                                    double & Prior ) const;
                                    
  virtual ErrorType LoadDistributions();
}; // end class InterModel

} // end namespace madai

#endif