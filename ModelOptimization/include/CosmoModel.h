#ifndef __CosmoModel_h__
#define __CosmoModel_h__

#include "MultiModel.h"

namespace madai {
  
class CosmoModel : public MultiModel {
public:
  CosmoModel();
  CosmoModel(std::string info_dir);
  virtual ~CosmoModel();
  
  
  virtual ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                                      std::vector< double > & scalars ) const;
  
  // Not implemented yet.  Should we do this numerically?
  /** Get both scalar values and the gradient of the parameters. */
  virtual ErrorType GetScalarAndGradientOutputs(const std::vector< double > & parameters,
                                                const std::vector< bool > & activeParameters,
                                                std::vector< double > & scalars,
                                                unsigned int outputIndex, std::vector< double > & gradient) const;
  
  // Proposed function for interaction with the MCMC
  /** Get the likelihood and prior at the point parameters in parameter space. */
  virtual ErrorType GetLikeAndPrior( const std::vector< double > & parameters,
                                    double & Like,
                                    double & Prior ) const;
  
  virtual ErrorType LoadDistributions();
  
}; // end class CosmoModel
  
} // end namespace madai

#endif // end __CosmoModel_h__