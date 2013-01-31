#ifndef __TestModel_h__
#define __TestModel_h__

#include "MultiModel.h"

namespace madai {

/** \class TestModel
 *
 * Test model. */
class TestModel : public MultiModel {
public:
  TestModel();
  TestModel(std::string info_dir);
  virtual ~TestModel();
  
  virtual ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                                      std::vector< double > & scalars ) const;
  
  // Not Implemented Yet
  virtual ErrorType GetScalarAndGradientOutputs( const std::vector< double > & parameters,
                                                 const std::vector< bool > & activeParameters,
                                                 std::vector< double > & scalars,
                                                 unsigned int outputIndex, std::vector< double > & gradient) const;
  
  // Interface for MCMC
  virtual ErrorType GetLikeAndPrior( const std::vector< double > & parameters,
                                     double & Like,
                                     double & Prior ) const;
                                     
  virtual ErrorType LoadDistributions();
  
}; // end class TestModel

}// end namespace madai

#endif // end __TestModel_h__
