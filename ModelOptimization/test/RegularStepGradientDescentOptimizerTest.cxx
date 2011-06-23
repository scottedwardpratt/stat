#include <cstdlib>
#include <iostream>

#include "Gaussian2DModel.h"
#include "RegularStepGradientDescentOptimizer.h"
#include "Trace.h"


int main(int argc, char *argv[])
{
  madai::Gaussian2DModel *model =
    new madai::Gaussian2DModel();
  model->LoadConfigurationFile( "file.txt" ); // TODO - does nothing

  madai::RegularStepGradientDescentOptimizer *optimizer =
    new madai::RegularStepGradientDescentOptimizer( model );

  madai::Trace *trace = NULL; // TODO - trace is not yet defined.

  // Set the step size.
  double stepSize = 20.0;
  optimizer->SetStepSize( stepSize );

  // Pick which output scalar to optimize.
  optimizer->SetOutputScalarToOptimize( "Value " );

  // Set initial parameter values.
  optimizer->SetParameterValue( "X", 21.0 );
  optimizer->SetParameterValue( "Y", -13.5 );

  for (unsigned int i = 0; i < 50; i++)
  {
    optimizer->NextIteration( trace );

    std::vector< double > currentParameters =
      optimizer->GetCurrentParameters();
    std::cout << "[";
    unsigned int j;
    for ( j = 0; j < currentParameters.size()-1; ++j )
    {
      std::cout << currentParameters[j] << ", ";
    }
    std::cout << currentParameters[j] << "]" << std::endl;
  }

  // TODO - test for convergence

  return EXIT_SUCCESS;
}
