#!/bin/zsh
# ccs, cec24@phy.duke.edu
# Create hyperparameters for a given data set
# wraps up the R file makeThetas.R
infile=$1
errfile=$2
paramfile=$3

export infile errfile paramfile
R --slave  < ./makeThetas.R
