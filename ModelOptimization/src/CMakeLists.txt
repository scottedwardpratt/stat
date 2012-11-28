set( SRC_FILES
  Optimizer.cxx
  MCMCRun.cxx
  Trace.cxx

  Distributions/Distribution.cxx
  Distributions/LikelihoodRHIC.cxx
  Distributions/PriorRHIC.cxx
  Distributions/LikelihoodCosmo.cxx
  Distributions/PriorCosmo.cxx
  Distributions/LikelihoodRHIC_PCA.cxx
  Distributions/PriorRHIC_PCA.cxx
  Distributions/LikelihoodTest.cxx
  Distributions/PriorTest.cxx
  Distributions/Proposal.cxx
  Quad.cxx

  MultiModel.cxx

  Visualization.cxx

  ExternalModel.cxx
  Gaussian2DModel.cxx
  RegularStepGradientDescentOptimizer.cxx
  SimpleMetropolisHastings.cxx
  process_pipe.c
)

include_directories(
  ../include
  /opt/local/include
  ~/local/include
)

set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ModelOptimization_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
add_library( IHU ${SRC_FILES} )
