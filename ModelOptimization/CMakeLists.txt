project(ModelOptimization)
cmake_minimum_required(VERSION 2.8)

add_subdirectory( src )

option( BUILD_TESTING "Build the tests for ModelOptimization." ON )
if ( BUILD_TESTING )
  enable_testing()
add_subdirectory( test )
endif()

  