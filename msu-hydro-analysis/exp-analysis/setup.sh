#!/bin/zsh
## ccs, cec24@phy.duke.edu
## create symlinks to a given model-data analysis and parameters folder in each of the subfolders
## this allows a relatively automatic visual analysis of the default hbt, spec-and-flow and spec-only
## output
anaFolders=(hbt hydro-spec spectra)

## check that we have $1
if [ $# -lt 1 ]; then 
		echo "produces links for a visual emulator analysis"
		echo "supply the top-level model-data folder, NOT an analysis or parameters folder"
		echo "Usage: `basename $0` {targetFolder}"
		exit -1
fi



echo "creating links to: " $1

top=`pwd`
for j in $anaFolders; do
		cd $j
		if [ -L ./design ]; then 
				echo "# removing existing design link"
				rm ./design
		elif [ -e ./design -a -d ./design ]; then
				echo "a real design folder exists in: $j"
				exit -1
		fi
				
		echo "# creating new design link"
		ln -s $1/parameters ./design

		if [ -L ./model-data ]; then 
				echo "# removing existing model-data link"
				rm ./model-data
		elif [ -e ./model-data -a -d ./model-data ]; then
				echo "a real model-data folder exists in: $j"
				exit -1
		fi

		echo "# creating new model-data link"
		ln -s $1 ./model-data
		cd $top
done