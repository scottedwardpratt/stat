#!/usr/bin/python

# Prepare_for_Scott.py
# John Novak
# October 8, 2012 4:55 pm

# This code knocks together some '.dat' files for Scott's gunuplot script

import os
import sys
import math as m

def main():
    if not (len(sys.argv) == 3):
        print "This program requires two directories as input. The first must contain an emulator statefile (named Emulator.statefile), and the second must contain everything else: a data summary file (named DataSummary.dat), a folder containing the experimental data (named exp_data), and a pcanames file (and it should be the same as the one used for the emulator. These two directories can be the same"
        print "Example: ./Prepare_for_Scott.py model-data/trained_emulator/ model-data/emuatlor_data_file/"
        exit(0)

    Statefile = sys.argv[1] + "/Emulator.statefile"
    TestPoints = open(sys.argv[2] + "/DataSummary.dat", 'r')
    Expdata = open(sys.argv[2] + "/exp_data/results.dat", 'r')

    # Read in the names of the observables used in the data set
    obsvlist=open(sys.argv[2]+"/pcanames.dat","r")
    obsvnames_1=[]
    for line in obsvlist:
        if line[0] != "#":
            line=line[:-1]
            if line[-1] == '\r':
                line=line[:-1]
            if line[-1] == ' ':
                line=line[:-1]
            obsvnames_1.append(line)
    obsvlist.close()

    # Read in the names of the observables used in the emulator
    obsvlist_2=open(sys.argv[1]+"/pcanames.dat","r")
    obsvnames_2=[]
    for line in obsvlist_2:
        if line[0] != "#":
            line=line[:-1]
            if line[-1] == '\r':
                line=line[:-1]
            if line[-1] == ' ':
                line=line[:-1]
            obsvnames_2.append(line)
    obsvlist_2.close()

    # Compare observables lists and figure out what we are going to use
    same = True
    if len(obsvnames_1) != len(obsvnames_2):
        same = False
    else:
        for i in range(len(obsvnames_1)):
            if obsvnames_1[i] != obsvnames_2[i]:
                same = False
    if not same:
        print "The observables used in the emualtor and in the data are not the same:"
        print len(obsvnames_1),obsvnames_1
        print len(obsvnames_2),obsvnames_2
        print "We will use the intersection of the two lists:"
        obsvnames = [val for val in obsvnames_1 if val in obsvnames_2]
    else:
        obsvnames = obsvnames_1
        print "Observables:"
    print len(obsvnames),obsvnames

    # Read in the observables scales
    obsvscalesfile_1 = open(sys.argv[1]+"/ObservableScales.dat",'r')
    obsvscales_1 = []
    for line in obsvscalesfile_1:
        for i in obsvnames:
            if line.split(" ")[0] == i:
                obsvscales_1.append([float(line.split(" ")[1]),float(line[:-1].split(" ")[2])])
    obsvscalesfile_2 = open(sys.argv[2]+"/ObservableScales.dat",'r')
    obsvscales_2 = []
    for line in obsvscalesfile_2:
        for i in obsvnames:
            if line.split(" ")[0] == i:
                obsvscales_2.append([float(line.split(" ")[1]),float(line[:-1].split(" ")[2])])

    # Read in the experimental data and unscale it
    exp_data = []
    exp_file = open(sys.argv[2]+"/exp_data/results.dat",'r')
    j=0
    for line in exp_file:
        for i in obsvnames:
            if line.split(" ")[1] == i:
                exp_data.append((float(line.split(" ")[2])*obsvscales_2[j][1])+obsvscales_2[j][0])
                j += 1

    # Reading in data from DataSummary.dat
    num_obsv = int(TestPoints.readline())
    num_params = int(TestPoints.readline())
    num_points = int(TestPoints.readline())
    Xpoints = []
    Ypoints = []
    for i in range(num_points):
        line = TestPoints.readline()
        Xpoints.append(line)  # Here we are taking the points as strings, just as they are because that is how the emulator will want them
        if len(Xpoints[-1][:-2].split(" ")) != num_params:
            print "Parsing error. The number of observables stated in " + sys.argv[2] + "/DataSummary.dat is not the same as the numbe read in. Problem with point", i
            print "Expected", num_params, "observables. Got", len(Xpoints[-1][:-2].split(" ")),Xpoints[-1][:-2].split(" "), "Exiting"
            exit(1)
    for i in range(num_points):
        line = TestPoints.readline()
        Ypoints.append(map(float, line[:-2].split(" ")))  # We turn these into an array of nums so we can do math with them
        if len(Ypoints[-1]) != num_obsv:
            print "Parsing error. The number of parameters stated in " + sys.argv[2] + "/DataSummary.dat is not the same as the number read in. Problem with point", i, "Exiting"
            print "Expected", num_obsv, "observables. Got", len(Ypoints[-1]),Ypoints[-1], "Exiting"
            exit(1)

    # Print just the X points to file for the emulator
    print "There are",len(Xpoints),"points to process"
    pointstofile = open("EmulatedPoints.txt",'w')
    for i in Xpoints:
        pointstofile.write(i)
    pointstofile.close()

    # Run Emulator on points and dump output to file
    os.system("interactive_emulator interactive_mode -q "+Statefile+" < EmulatedPoints.txt > EmulatedObservables.txt")

    # Read emulator output
    means = []
    errors = []
    longdiff = []
    diffsum = 0
    index1 = 0
    observablesfromfile = open("EmulatedObservables.txt",'r')
    print "num_obsv =",num_obsv,". number of observables used=",len(obsvnames)
    for i in Xpoints:
        meanstemp=[]
        errorstemp=[]
        count = 0
        for j in range(num_obsv):
            if obsvnames[count] == obsvnames_1[j]:
                meanstemp.append((float(observablesfromfile.readline())*obsvscales_1[count][1])+obsvscales_1[count][0])
                errorstemp.append(float(observablesfromfile.readline()))
                count += 1
        means.append(meanstemp)
        errors.append(errorstemp)
    observablesfromfile.close()

    # Run through Ypoints and toss the observables we aren't using. And unscale everything
    for i in range(len(Ypoints)):
        temp = []
        count = 0
        for j in range(len(Ypoints[i])):
            if obsvnames[count] == obsvnames_1[j]:
                temp.append((Ypoints[i][j]*obsvscales_1[count][1])+obsvscales_1[count][0])
                count += 1
        Ypoints[i] = temp

    # Burn through the points and make the output
    output = open("quadcheck_num.dat",'w')
    output.write("#irun netdiff_quad netdiff_exp netdiff_quadexp\n") #runnum model-emu model-exp emu-exp
    #print exp_data
    for i in range(num_points):
        model_emu = 0
        model_exp = 0
        emu_exp   = 0
        #print Ypoints[i],means[i]
        for j in range(len(obsvnames)):
            model_emu += m.pow((means[i][j]-Ypoints[i][j])/(obsvscales_1[j][1]),2)/2
            model_exp += m.pow((means[i][j]-exp_data[j])/(obsvscales_1[j][1]),2)/2
            emu_exp   += m.pow((exp_data[j]-Ypoints[i][j])/(obsvscales_1[j][1]),2)/2
        #output.write(str(i+1)+" "+str(model_emu)+" "+str(model_exp)+" "+str(emu_exp)+"\n")
        output.write(str(i+1)+" "+str(model_emu)+" "+str(model_exp)+" "+str(emu_exp)+"\n")
    output.close()

if __name__ == '__main__':
    main()
