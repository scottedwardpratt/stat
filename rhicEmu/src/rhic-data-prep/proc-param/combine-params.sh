#!/bin/zsh
# grab the data
#
# find the folders
# 
# $1 is the path to process
#TARGETFOLDER=$1
for folder in $(find . -maxdepth 1 -type d -name 'run*'); do
		#find the folders
		#TEMP=`echo $folder | cut -c 6-7`
		# good old awk, this works because we know that the folders are called 
		# run1, run2, ... if they're not then you're boned
		TEMP=`echo $folder | awk '{split($0,f,"n"); print f[2]}'`
		echo $TEMP
		cd $folder
		echo -n "$TEMP\t" > temp
		cat stats.param | awk '{print $3}' | awk 'ORS=NR%7?" ":"\n"' >> temp
		cd ..
done

# now find all the temp files, cat them together and sort the results
find . -name 'temp' | xargs cat | sort -n > stats.param-combined-s
# tidy up
find . -name 'temp' -delete

