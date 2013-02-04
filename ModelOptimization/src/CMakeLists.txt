set( SRC_FILES
  Optimizer.cxx
  MCMCRun.cxx
  Trace.cxx
  TraceElement.cxx

  Distribution.cxx
  LikelihoodDistribution.cxx
  PriorDistribution.cxx
  RHIC_LikelihoodDistribution.cxx
  RHIC_PriorDistribution.cxx
  CosmoLikelihoodDistribution.cxx
  CosmoPriorDistribution.cxx
  RHIC_PCA_LikelihoodDistribution.cxx
  RHIC_PCA_PriorDistribution.cxx
  TestLikelihoodDistribution.cxx
  TestPriorDistribution.cxx

  MultiModel.cxx
  RHICModel.cxx
  CosmoModel.cxx
  RHIC_PCA_Model.cxx
  TestModel.cxx
  Quad.cxx

  DiceModel.cxx

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
