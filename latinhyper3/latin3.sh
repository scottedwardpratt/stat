#!/bin/sh
# wrapper for the r sample generator
# ccs aug-10, cec24@phy.duke.edu
# RELEASED UNDER NO WARRANTY
case $# in
0|1)    echo "Usage: latin3.sh npars nruns";
        exit 1 ;;
*)
rm -f latin3.dat
npars=$1
nruns=$2
export npars nruns
#R CMD BATCH < ./genSamples.R
(R --slave < ./genSamples.R) > latin3.dat
esac