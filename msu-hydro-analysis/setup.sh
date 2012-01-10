#!/bin/zsh
## ccs, cec24@phy.duke.edu
## setup.sh
##
## this script creates a folder in model-data from fresh parameters and analysis folders to hold all 
## the data for a given set of runs
## 
## data relating to the analysis of a code is stored in a model-data folder, this data includes
## the parameter values for each run (the design)
## the model output, a lot of observables (the training data)
## scratch data created by the emulator
##
## first the design data is extracted in the parameters folder, the final design file is stats.param-combined-s
## 
## model results are extracted by centrality into the analysis folder, the calls to generate-analysis-files.rb
## produce custom sets of combined data. 
## The arguments to generate-analysis-files.rb are: <s/analysis-files-path> <I:cclass> <s:output-file-path> [<s:output-error-file-path>] 
## 
## analysis-files-path: folder which contains run-1, ..., run-n folders containing runs and results.dat
## cclass: integer setting upper centrality bin: 1 -> 0 to 10, 2 -> 10 to 20 etc.
## 
## This script also takes specification files which select which observables
## to glue together, viz ./data-prep/spectra-only-spec.dat, ./data-prep/star-spectra-and-flow.dat
## 
## we need to pick which centrality classes we are going to extract, the defaults are 
## cent0to10, cent10to20, cent20to30, 
## to adjust this you'll need to change the arguments to generate-analysis-files.rb


## check that we have $1
if [ $# -lt 1 ]; then 
		echo "produces folders needed for an emulator analysis"
		echo "Usage: `basename $0` {targetFolder}"
		exit -1
fi

echo "# processing parameters"
./data-prep/process-parameters.sh $1


echo "# processing analysis files"
echo "# creating default combined analysis files"
## we'll create files for cents 1 2 3 for hbt,spec and spec and flow

## if you want to add another output type you should add its name here
outputNames=(spectra-combined hbt-combined spec-and-flow-combined)
## here you can add another centrality class, we currently have data for 30to40 and 40to50
centNames=(0to10 10to20 20to30)
## and you should add the name of the specification file here
specFiles=(spectra-only-spec.dat hbt-only-spec.dat star-spectra-and-flow.dat)

countSpec=1
for i in $outputNames; do
		countCent=1
		for j in $centNames; do
				./data-prep/generate-analysis-files.rb $1/analysis $countCent $1/$i-$j.dat $1/$i-$j-errors.dat <./data-prep/$specFiles[$countSpec] $countSpec $countCent
				countCent=$((countCent+1))
		done
		countSpec=$((countSpec+1))
done

