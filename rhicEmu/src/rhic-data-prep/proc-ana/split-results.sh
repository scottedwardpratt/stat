#!/bin/sh
#
# ccs, cec24@phy.duke.edu (march 11)
# 
# From the root of a set  of analysis folders, loop over the run folders and for each impact 
# parameter cut up the results.dat file with the given awk scripts.
# 
# you may need to customize the awk scripts, their location and the name of the file to be cut.
# 
# on the shore of our lives
# the world blooms for the last time
FILENAME=results.dat
TOP=`pwd`
AWKROOT=${MADAI_ROOT}/src/rhic-data-prep/proc-ana
HBT=${AWKROOT}/hbt.awk
SPEC=${AWKROOT}/separate-spectra.awk
FLOW=${AWKROOT}/separate-flow.awk


for folder in $(find . -maxdepth 1 -type d -name 'run*'); do
		TEMP=`echo $folder | cut -c 6-7`
		echo "processing: ${folder}"
		cd $folder
		for impact in $(find . -maxdepth 1 -type d -name 'b*'); do
				echo $impact
				cd $impact
				IMPCUT=`echo $impact | cut -c 3-7`
				#echo "hbt-run-${TEMP}-${IMPCUT}.dat"
				awk -f $HBT results.dat > hbt-run-${TEMP}-${IMPCUT}.dat
				awk -f $SPEC results.dat > spec-run-${TEMP}-${IMPCUT}.dat
				awk -f $FLOW results.dat > flow-run-${TEMP}-${IMPCUT}.dat
				cd ..
		done
		cd ..
done
