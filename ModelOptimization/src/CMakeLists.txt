set( SRC_FILES
  Optimizer.cxx
  MCMCRun.cxx
  Trace.cxx

  Distribution.cxx
  LikelihoodDistribution_RHIC.cxx
  PriorDistribution_RHIC.cxx
  LikelihoodDistribution_Cosmo.cxx
  PriorDistribution_Cosmo.cxx
  LikelihoodDistribution_RHIC_PCA.cxx
  PriorDistribution_RHIC_PCA.cxx
  LikelihoodDistribution_Test.cxx
  PriorDistribution_Test.cxx

  MultiModel.cxx
  RHICModel.cxx
  CosmoModel.cxx
  RHIC_PCA_Model.cxx
  TestModel.cxx

  Visualization.cxx

  ExternalModel.cxx
  Gaussian2DModel.cxx
  RegularStepGradientDescentOptimizer.cxx
  SimpleMetropolisHastings.cxx
  process_pipe.c
)

set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ModelOptimization_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
add_library( ModelOptimization ${SRC_FILES} )
