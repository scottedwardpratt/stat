#!/bin/zsh
#
# ccs, cec24@phy.duke.edu
# computePoints: 
# the argument is the folder to do the analysis on
# input (stdin): a list of points in the paramater space for the emulator to be evaluated at
# output (stdout): for each point the ntps mean cpts and then the ntps var cpts (on a single line)
# 
# ccs, may-19th: added reconstruct. If 1 we will do the project the emulated data back into the 
# functional space (using the call reconCurveAtList) otherwise we return the unprojected data
# via a call to emulateAtListNoProject
#
#
reconstruct=1 ## do we want to reconstruct the data back into the functional space? (1 = yes)
basepath=$1 # the folder where the model output is stored and the various temp files from the emulator
anadata=$2 # the R-data saved by estimate thetas

if [ $# -lt 1 ]; then
		echo "#computes predicted mean and variance"
		echo "#reads points to process from stdin"
		echo "#usage: `basename $0` {targetFolder} {Saved-R-analysisFile}"
		exit -1
fi

if [ -d $basepath -a -d $basepath/analysis -a -d $basepath/parameters ]; then
# this is most likely a good path
		echo "#computing values for: ${basepath}"
else 
		echo "#arg: ${basepath} doesn't exist"
		exit -1
fi

if [ -e $anadata ]; then
		echo "#found results from estimation"
else 
		echo "#no data files from estimation. Run make-thetas.sh on ${basepath} first"
		exit -1
fi
		
inpath=$basepath
fndata=$anadata

EMU_DIRECTORY=$(cd `dirname $0` && pwd)

export inpath fndata reconstruct EMU_DIRECTORY
cat /dev/stdin > temp
# this will read temp by default
# but we don't get output on the stdout
#R CMD BATCH --slave computePoints.R 
R --slave --no-save < $EMU_DIRECTORY/computePoints.R
rm temp