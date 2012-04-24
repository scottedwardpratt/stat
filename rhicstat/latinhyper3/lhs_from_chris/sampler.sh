#!/bin/sh
# wrapper for the r sample generator
# ccs aug-10, cec24@phy.duke.edu
# RELEASED UNDER NO WARRANTY
npars=$1
nruns=$2

export npars nruns

#R CMD BATCH < ./genSamples.R
R --slave < ./genSamples.R 