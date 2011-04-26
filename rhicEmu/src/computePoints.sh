#!/bin/zsh
#
# ccs, cec24@phy.duke.edu
# computePoints: 
# the argument is the folder to do the analysis on
# input (stdin): a list of points in the paramater space for the emulator to be evaluated at
# output (stdout): for each point the ntps mean cpts and then the ntps var cpts (on a single line)
# 

rescale=0 # rescale the output back to correct values, if you scaled the data without the errors (not sensible)
rescaleErr=1 # same but only if you scaled the data with errors (DEFAULT)
basepath=$1 # the folder where the model output is stored and the various temp files from the emulator

if [ $# -lt 1 ]; then
		echo "#computes mean and variance for a run-folder"
		echo "#reads points to process from stdin"
		echo "#usage: `basename $0` {targetFolder}"
		exit -1
fi

if [ -d $basepath -a -d $basepath/analysis -a -d $basepath/parameters ]; then
# this is most likely a good path
		echo "#computing values for: ${basepath}"
else 
		echo "#arg: ${basepath} doesn't exist"
		exit -1
fi

if [ -e $basepath/PcaTemp.dat -a -e $basepath/SampleTemp.dat -a -e $basepath/theta-table.dat ]; then
		echo "#found results from estimation"
else 
		echo "#no temp files from estimation. Run make-thetas.sh on ${basepath} first"
		exit -1
fi
		
inpath=$basepath

EMU_DIRECTORY=$(cd `dirname $0` && pwd)

export rescale rescaleErr inpath EMU_DIRECTORY
cat /dev/stdin > temp
# this will read temp by default
# but we don't get output on the stdout
#R CMD BATCH --slave computePoints.R 
R --slave --no-save < $EMY_DIRECTORY/computePoints.R
rm temp