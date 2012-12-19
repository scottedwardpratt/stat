#ifndef __RHIC_PCA_Model_h__
#define __RHIC_PCA_Model_h__

#include "MultiModel.h"
#include "Quad.h"

namespace madai {
  
class RHIC_PCA_Model : public MultiModel {
private:
  QuadHandler* m_Quad;
public:
  RHIC_PCA_Model();
  RHIC_PCA_Model(std::string info_dir);
  virtual ~RHIC_PCA_Model();
    
    
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
  
}; // end class RHIC_PCA_Model
  
} // end namespace madai

#endif // end __RHIC_PCA_Model_h__