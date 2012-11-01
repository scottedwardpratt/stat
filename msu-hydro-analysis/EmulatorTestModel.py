#!/usr/bin/python

# EmulatorTestModel.py

# John Novak
# Wednesday, October 31, 2012, 2:32pm

# This code is meant to be as a test for the Gaussian Process Emultator
# It generates a set of test points for the emulator to be trained on

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


import sys
import os
import numpy as np
import random
import math as m

def main():
    if len(sys.argv) < 2:
        print "This program makes the name of a model-directory and fills it with data from a simple model for testing the MADAI emulator and MCMC code. It requires two inputs: the model-directory and the number of training points desired."
        print "Example: ./EmultaorTestModel.py model-data/Test_Directory 250"
        exit(1)

    directory = sys.argv[1]
    Number_of_runs = int(sys.argv[2])

    SetUpFiles(directory)

    for i in range(1,Number_of_runs+1):
        os.system("mkdir -p "+directory+"/model_results/run"+str(i))
        os.system("mkdir -p "+directory+"/parameters/run"+str(i))
        x = [random.random(),random.random()]
        y = Function(x)
        obser = open(directory+"/model_results/run"+str(i)+"/results.dat",'w')
        param = open(directory+"/parameters/run"+str(i)+"/stats.param",'w')
        for j in range(len(x)):
            param.write("double param"+str(j)+" "+str(x[j])+" 0.1 \n")
        for j in range(len(y)):
            obser.write("double obser"+str(j)+" "+str(y[j])+" \n")
        param.close()
        obser.close()
    os.system("cp -r "+directory+"/model_results/run1 "+directory+"/model_results/default")
    os.system("cp -r "+directory+"/parameters/run1 "+directory+"/model_results/default")

    x = [0.5,0.5]
    y = Function(x)
    exp = open(directory+"/exp_data/results.dat",'w')
    pca = open(directory+"/pcanames.dat",'w')
    ranges = open(directory+"/ranges.dat",'w')
    for j in range(len(x)):
        ranges.write("double param"+str(j)+" 0 1 \n")
    for j in range(len(y)):
        exp.write("double obser"+str(j)+" "+str(y[j])+" -999 \n")
        pca.write("obser"+str(j)+" \n")
    exp.close()
    pca.close()
    ranges.close()

def SetUpFiles(directory):
#This fucntion takes a directory as a sting and makes the necessary file structure there
    os.system("mkdir -p "+directory+"/model_results")
    os.system("mkdir -p "+directory+"/parameters")
    os.system("mkdir -p "+directory+"/exp_data")
    os.system("mkdir -p "+directory+"/defaultpars")

def Functiona(x):
    """ This is the underlying function of the data which the emulator will try and emulate.
    It takes a vector x and returns a vector y."""
    m = np.array([[2,3],[4,5]])
    y = np.dot(m,x)
    return y

def Function(x):
    """ This is the underlying function of the data which the emulator will try and emulate.
    It takes a vector x and returns a vector y."""
    y = []
    y.append(m.exp(-(((x[0]-.25)*(x[0]-.25))+((x[1]-.25)*(x[1]-.25)))/(2*.25)))
    y.append(m.exp(-(((x[0]-.75)*(x[0]-.75))+((x[1]-.75)*(x[1]-.75)))/(2*.25)))
    return y

if __name__ == '__main__':
    main()
