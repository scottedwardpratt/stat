set( SRC_FILES
  Gaussian2DModel.cxx
  Optimizer.cxx
  RegularStepGradientDescentOptimizer.cxx
)

include_directories(
  ../include
)

set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ModelOptimization_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} )
add_library( ModelOptimization ${SRC_FILES} )
