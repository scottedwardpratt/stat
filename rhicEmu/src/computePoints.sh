#!/bin/sh
#
# ccs, cec24@phy.duke.edu
# computePoints: 
# the argument is the folder to do the analysis on
# input (stdin): a list of points in the paramater space for the emulator to be evaluated at
# output (stdout): for each point the ntps mean cpts and then the ntps var cpts (on a single line)
# 

#cat testPts.txt 

# return the output values to the same scale as the input ? (defaults to yes)
rescale=0
rescaleErr=1
basepath=$1

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
		echo "#arg: ${basepath} doesn't exiset"
		exit -1
fi

if [ -e $basepath/PcaTemp.dat -a -e $basepath/SampleTemp.dat -a -e $basepath/theta-table.dat ]; then
		echo "#found results from estimation"
else 
		echo "#no temp files from estimation. Run make-thetas.sh on ${basepath} first"
		exit -1
fi
		
inpath=$basepath

export rescale rescaleErr inpath
cat /dev/stdin > temp
# this will read temp by default
# but we don't get output on the stdout
#R CMD BATCH --slave computePoints.R 
R --slave --no-save < computePoints.R
rm temp