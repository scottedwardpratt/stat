set( SRC_FILES
  Optimizer.cxx
  Trace.cxx

  Gaussian2DModel.cxx
  RegularStepGradientDescentOptimizer.cxx
  ExternalModel.cxx
  process_pipe.c
  SimpleMetropolisHastings.cxx
)

include_directories(
  ../include
)

set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ModelOptimization_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
add_library( ModelOptimization ${SRC_FILES} )
