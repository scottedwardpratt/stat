project(ModelOptimization)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
link_directories(
  ~/local/lib
)

set(EXT_LIBS
  -lEmuPlusPlus
  -lcoralutils
)

find_package(GSL REQUIRED)
find_package(X11)
find_package(HDF5)
set(CORELIBS
  ${GSL_LIBRARIES}
  ${X11_LIBRARY}
  ${HDF5_LIBRARY}
  ${EXT_LIBS}
)
add_subdirectory( src )
add_subdirectory( main )

option( BUILD_TESTING "Build the tests for ModelOptimization." OFF )
if ( BUILD_TESTING )
  enable_testing()
add_subdirectory( test )
endif()

  
