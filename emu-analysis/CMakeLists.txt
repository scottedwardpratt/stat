cmake_minimum_required (VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
project (emu-analysis)

################### 
# the install prefix is where the bins and libs get dropped. 
message("Setting Install prefix to Home folder: $ENV{HOME}")
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local")

## just install shit for now
INSTALL_FILES("/include/emu-analysis/" FILES src/fnAnalysis.R src/fnHyperAnalysis.R 
	src/fnMcSampling.R)
