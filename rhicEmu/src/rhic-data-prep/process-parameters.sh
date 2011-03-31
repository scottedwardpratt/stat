#!/bin/sh
# ccs, cec24@phy.duke.edu
# process the given parameter folder $1 to produce a stats.param-combined-s file
# basically just run combine-params.sh on it

# check for $1
if [ $# -lt 1 ]; then 
		echo "produces stats.param-combined-s."
		echo "Usage: `basename $0` {targetFolder}"
		exit -1
fi

if [ -d $1 ]; then 
		if [ -d $1/run1 ]; then
				echo "processing params from: $1"				
		else 
				echo "provide a params directory as argument"
				echo "${1}"
				exit -1
		fi
else 
		echo "provide a params directory as argument"
		echo "${1}"
		exit -1
fi



TOP=`pwd`
SCRIPT_PATH=${MADAI_ROOT}/src/rhic-data-prep
cd $1
${SCRIPT_PATH}/proc-param/combine-params.sh
cd $TOP