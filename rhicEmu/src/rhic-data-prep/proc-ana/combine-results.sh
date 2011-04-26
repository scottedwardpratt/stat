#!/bin/zsh
# take the data split up by split-results and combine into 3 top level files
# each run has 5 impact paramters and will give 5 lines in the final top level file
#
# find the folders
#FILENAME=results.dat
#IMPACT=b9.87

TOP=`pwd`
FIRST=0

DATFILEPREF=("spec" "hbt" "flow")
RESFILEPREF=("spectra" "hbt" "flow")

# note that we're doing each impact parameter
INFOSTRING1="# NB each set of variables appears 5 times, once for each impact param\n#b={3.29,5.70,7.36,8.70,9.87}"
for i in $RESFILEPREF; do
		echo ${INFOSTRING1} >> ${TOP}/${i}-combined.dat
		echo ${INFOSTRING1} >> ${TOP}/${i}-combined-errs.dat
done

# loop over all the runs
for folder in $(find . -maxdepth 1 -type d -name 'run*'| sort -n); do
		#find the folders
		# good old awk, this works because we know that the folders are called 
		# run1, run2, ... if they're not then you're boned
		TEMP=`echo $folder | awk '{split($0,f,"n"); print f[2]}'`
		echo "processing: ${folder}"

		# print headers (without leading # so R can read them)
		if [[ FIRST -eq 0 ]]; then
				cd $folder/b3.29
				IMPCUT=b3.29
				COUNT=1
				for i in $RESFILEPREF; do								
						FILENAME=${DATFILEPREF[$COUNT]}-run-${TEMP}-${IMPCUT}.dat
								#echo $FILENAME
						cat $FILENAME | awk 'BEGIN{printf("#ID\tIMP\t")};{printf("%s\t", $1)};END{printf("\n")}' >> ${TOP}/${i}-combined.dat
						cat $FILENAME | awk 'BEGIN{printf("#ID\tIMP\t")};{printf("%s\t", $1)};END{printf("\n")}' >> ${TOP}/${i}-combined-errs.dat
						COUNT=$(( $COUNT + 1 ))
				done
				FIRST=1
				cd ../..
		fi

		# print the run-id
		for i in $RESFILEPREF; do
				echo -n "${TEMP}  " >> ${TOP}/${i}-combined.dat
				echo -n "${TEMP}  " >> ${TOP}/${i}-combined-errs.dat
		done

		# now loop over the impact parameters and append them one by one
		cd $folder
		for impact in $(find . -maxdepth 1 -type d -name 'b*' | sort -n); do
				#echo $impact
				cd $impact
				IMPCUT=`echo $impact | cut -c 3-7`

				# we want all the impact parameters to be on the same line
				COUNT=1
				for i in $RESFILEPREF; do
						FILENAME=${DATFILEPREF[$COUNT]}-run-${TEMP}-${IMPCUT}.dat
						#echo $FILENAME
						cat $FILENAME | awk '{printf("%lf ", $2)};' >> ${TOP}/${i}-combined.dat
						cat $FILENAME | awk '{printf("%lf ", $3)};' >> ${TOP}/${i}-combined-errs.dat
						COUNT=$(( $COUNT + 1 ))
				done
				cd ..
		done 

		# now we need to echo a newline to each file
		for i in $RESFILEPREF; do
				echo -n "\n" >> ${TOP}/${i}-combined.dat
				echo -n "\n" >> ${TOP}/${i}-combined-errs.dat
		done
		cd ..
done



