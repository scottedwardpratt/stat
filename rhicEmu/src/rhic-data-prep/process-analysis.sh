#!/bin/sh
# ccs, cec24@phy.duke.edu
# process the given analysis folder $1 to produce all the analysis files we need
# runs split-results.sh, combine-results.sh and then sorty.sh
# probably needs zsh :(


# check for $1
if [ $# -lt 1 ]; then 
		echo "produces {hbt|spectra|flow}-combined{-err}-s.dat files"
		echo "Usage: `basename $0` {targetFolder}"
		exit -1
fi

# some bad but basic error checking
if [ -d $1 ]; then 
		if [ -d $1/run1 ]; then
				echo "process analysis files from: $1"				
		else 
				echo "provide an analysis directory as argument"
				exit -1
		fi
else 
		echo "provide an analysis directory as argument"
		exit -1
fi

TOP=`PWD`
SCRIPT_PATH=${MADAI_ROOT}/src/rhic-data-prep
cd $1
${SCRIPT_PATH}/proc-ana/split-results.sh
${SCRIPT_PATH}/proc-ana/combine-results.sh
${SCRIPT_PATH}/proc-ana/sorty.sh
cd $TOP


