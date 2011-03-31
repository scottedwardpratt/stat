#!/bin/sh
#
# ccs, cec24@phy.duke.edu
# computePoints: 
# input (stdin): a list of points in the paramater space for the emulator to be evaluated at
# output (stdout): for each point the ntps mean cpts and then the ntps var cpts (on a single line)
# 

#cat testPts.txt 

# return the output values to the same scale as the input ? (defaults to yes)
export rescale=0
export rescaleErr=1

cat /dev/stdin > temp
# this will read temp by default
# but we don't get output on the stdout
#R CMD BATCH --slave computePoints.R 
R --slave --no-save < computePoints.R
rm temp