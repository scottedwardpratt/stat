#!/usr/bin/python -tt

# John Novak
# June 29, 2012

# This is a simple model which creates output for testing the MADAI MCMC code

import sys
import os
import random

# For personal refernce
# The file structure:   key: -directory, .file, @made by code
#
# ModelFile
# |.ranges.dat
# |.pcanames.dat
# |.Emulator.statefile@
# |.ObservableScales.dat@
# |.DataSummary.dat@
# |-defaultpars
# | |.mcmc.param
# | |.likelihood.param
# | |.proposal.param
# | |.theta0.param
# | |.actual.param
# |-exp_data
# | |.results.dat !This file must exist, but it is moved and replaced by the code
# | |.resutls_unscaled.dat@
# |-parameters
# | |-default
# | | |.fixed.param
# | | |.stats.param
# | |-run*
# | | |.fixed.param
# | | |.stats.param
# |-model_results
# | |-default
# | | |.results.dat
# | |-run*
# | | |.results.dat
# |-trace@
# | |.Assorted and filled in by code
# |-DensityPlots@
# | |-Assorted and filled in by code

def main():
    if len(sys.argv) < 2:
        print "This program makes a model-directory and fills it with data from a simple model for testing the MADAI emulator and MCMC code. It requires two inputs: the model-directory and the number of training points desired."
        print "Example: ./MCMCTestModel.py model-data/Test_Directory 250"
        exit(1)

    directory = sys.argv[1]
    Number_of_runs = int(sys.argv[2])
    Num_Model_inputs = 2;
    Num_Model_outputs = 1;

    SetUpFiles(directory)

    # For each run:
    for i in range(1,Number_of_runs+1):
        os.system("mkdir -p "+directory+"/parameters/run"+str(i))
        os.system("mkdir -p "+directory+"/model_results/run"+str(i))

        # Make the input params:
        x = []
        for j in range(Num_Model_inputs):
            x.append(random.random())
        # Calculate the observables:
        y = ModelFunction(x)

        # Write out the params and observables:
        params=open(directory+"/parameters/run"+str(i)+"/stats.param","w")
        for j in range(Num_Model_inputs):
            params.write("double param"+str(j)+" "+str(x[j])+" \n")
        params.close()
        results=open(directory+"/model_results/run"+str(i)+"/results.dat","w")
        for j in range(Num_Model_outputs):
            results.write("double Observable"+str(j)+" "+str(x[j])+" 0 \n")
        results.close()

    # Tie up loose ends:
    
    # Write the ranges.dat file:
    # Right now the range for everything is just 0 to 1
    ranges=open(directory+"/ranges.dat","w")
    for j in range(Num_Model_inputs):
        ranges.write("double param"+str(j)+" 0 1 \n")
    # Write the exp_data file:
    # For the moment we will just say that the "right spot" is in the middle
    x = []
    for j in range(Num_Model_inputs):
        x.append(0.5)
    y = ModelFunction(x)
    expdata=open(directory+"/exp_data/results.dat","w")
    for j in range(Num_Model_outputs):
        expdata.write("double Observable"+str(j)+" "+str(x[j])+" 0 \n")
    expdata.close()
    # Write the pcanames.dat file:
    pcanames=open(directory+"/pcanames.dat","w")
    for j in range(Num_Model_outputs):
        pcanames.write("Observable"+str(j)+"\n")
    pcanames.close()

    # Create "default" folders:
    os.system("cp -r "+directory+"/parameters/run1 "+directory+"/parameters/default")
    os.system("cp -r "+directory+"/model_results/run1 "+directory+"/parameters/default")

def SetUpFiles(directory):
#This fucntion takes a directory as a sting and makes the necessary file structure there
    os.system("mkdir -p "+directory+"/model_results")
    os.system("mkdir -p "+directory+"/parameters")
    os.system("mkdir -p "+directory+"/exp_data")
    os.system("mkdir -p "+directory+"/defaultpars")

def ModelFunction(X):
#This function is the underlying function of the fake model
#it takes a vector of inputs, then returns a vector of outputs
    Y=[]
    Y.append((X[0]*X[0])+(X[1]*X[1]))

    return Y

if __name__ == '__main__':
    main()
