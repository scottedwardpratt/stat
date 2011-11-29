#!/bin/zsh
# ccs, cec24@phy.duke.edu
# Create hyperparameters for a given data set
# wraps up the R file makeThetas.R
# 
# arguments are as follows:
# 
# $1 is the folder in $MADAI_ROOT/model-data for which we're going to do analysis
# $2 is enough of the name of the output we want to process, spectra, hbt etc
# 
# the intermediate data created by the R process, i.e the PCA decomp and the processed sample etc
# get left in basepath, which stops them getting nuked when you switch between which data-set you want
# to analyze etc.
basepath=$1 
analysis=$2

if [ $# -lt 2 ]; then
		echo "estimates hyperparamters for a dataset"
		echo "arguments: "
		echo "1) folder you want to process"
		echo "2) analysis file you want to train emulator on"
		echo "usage: `basename $0` {targetFolder} {result-name}"
		exit -1
fi


if [ -d $basepath -a -d $basepath/analysis -a -d $basepath/parameters ]; then 
		echo "# estimating thetas for: ${basepath}"
else 
		echo "arg: ${basepath} doesn't exist or is bad"
		exit -1
fi

paramfile=$basepath/parameters/stats.param-combined-s
if [ -e $paramfile ]; then
		echo "# found paramfile: ${paramfile}" ## good
else 
		echo "no paramfile found in: ${basepath}"
		exit -1
fi

infile=$2

if [ -e $infile ]; then
		echo "# found infile: ${infile}"
else 
		echo "problem with infile"
		echo "infile: ${infile}"
		exit -1
fi

outpath=$basepath

export infile paramfile outpath
R --slave  < ${MADAI_ROOT}/src/makeThetas.R
