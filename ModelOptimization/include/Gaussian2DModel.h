#ifndef __Gaussian2DModel_h_
#define __Gaussian2DModel_h_

#include "Model.h"

namespace madai {

class Gaussian2DModel : public Model {
public:
  Gaussian2DModel();
  ~Gaussian2DModel() {};

  /** Loads a configuration from a file. For this model, the
   * configuration file is a text file containing four numbers
   * separated by whitespace (space, tab, newline). The numbers are
   * the mean in x, mean in y, standard deviation in x and standard
   * devaition in y of the 2D Gaussian function. */
  ErrorType LoadConfigurationFile( const std::string fileName );
    
  /** Get the valid range for the parameter at parameterIndex. */
  void GetRange( unsigned int parameterIndex, double range[2] );

  /** Get the scalar outputs from the model evaluated at x. */
  ErrorType GetScalarOutputs( const std::vector< double > & parameters,
                              std::vector< double > & scalars ) const;

  /** Get both scalar values and the gradient of the parameters. */
  ErrorType GetScalarAndGradientOutputs( const std::vector< double > & parameters,
                                         const std::vector< bool > & activeParameters,
                                         std::vector< double > & scalars,
                                         unsigned int outputIndex, std::vector< double > & gradient) const;
  // Not implemented yet.
  // Proposed function for interaction with the MCMC
  virtual ErrorType GetLikeAndPrior( const std::vector< double > & parameters,
                                     double & Like,
                                     double & Prior ) const;
protected:
  double m_MeanX;
  double m_MeanY;
  double m_StandardDeviationX;
  double m_StandardDeviationY;

  double PartialX( double x, double value ) const;
  double PartialY( double y, double value ) const;

}; // end class Gaussian2DModel

} // end namespace madai

#endif // __Gaussian2DModel_h_
