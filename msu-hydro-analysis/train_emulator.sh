#!/bin/sh
#
# John Novak
# June 23, 2012
#
# Given the path to the model directory, this trains an emulator using a first order regression
# and puts the state file back into the model directory named "Emulator.statefile"

if [ $# -lt 1 ]; then
                echo "Trains emulator and generates statefile"
                echo "Usage: `basename $0` {Model Folder}"
                exit -1
fi

binpath=interactive_emulator
$binpath estimate_thetas $1/DataSummary.dat $1/Emulator.statefile --regression_order=1 --pca_variance=0.99
#$binpath estimate_thetas $1/DataSummary.dat $1/Emulator.statefile --regression_order=1
#$binpath estimate_thetas $1/DataSummary.dat $1/Emulator.statefile --regression_order=1 --min_length_scale=-1
#interactive_emulator estimate_thetas model-data/September2012/DataSummary.dat model-data/September2012/Emulator.statefile --regression_order=1 --min_length_scale=0.2
