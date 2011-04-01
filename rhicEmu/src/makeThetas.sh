#!/bin/zsh
# ccs, cec24@phy.duke.edu
# Create hyperparameters for a given data set
# wraps up the R file makeThetas.R
#
# $1 is the folder in $MADAI_ROOT/model-data for which we're going to do analysis
# $2 is enough of the name of the output we want to process, spectra, hbt etc
basepath=$1 
analysis=$2

if [ $# -lt 2 ]; then
		echo "estimates hyperparamters for a run-folder"
		echo "usage: `basename $0` {targetFolder} {result-name}"
		exit -1
fi


if [ -d $basepath -a -d $basepath/analysis -a -d $basepath/parameters ]; then 
		echo "estimating thetas for: ${basepath}"
else 
		echo "arg: ${basepath} doesn't exist or is bad"
		exit -1
fi

paramfile=$basepath/parameters/stats.param-combined-s
if [ -e $paramfile ]; then
		echo "found paramfile: ${paramfile}" ## good
else 
		echo "no paramfile found in: ${basepath}"
		exit -1
fi

infile=$2
# want to take spectra-combined-s.dat
# and sub it to spectra-combined-errs-s.dat
# sub in -errs 
errfile=`echo ${infile} | awk 'gsub("-s", "-errs-s", $0)'`

if [ -e $infile -a -e $errfile ]; then
		echo "found infile: ${infile}"
		echo "found errfile: ${errfile}"
else 
		echo "problem with infile or errfile"
		echo "infile: ${infile}"
		echo "errfile: ${errfile}"
		exit -1
fi

outpath=$basepath

export infile errfile paramfile outpath
R --slave  < ${MADAI_ROOT}/src/makeThetas.R
