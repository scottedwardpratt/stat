#!/bin/bash

# To use, run this from the top level of the ModelOptimization
# library source directory.

files=`ls -1 include/*.h src/*.cxx src/Distributions/*.cxx`

for file in $files; do KWStyle -v $file; done
