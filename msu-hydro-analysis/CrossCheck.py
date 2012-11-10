#!/usr/bin/python

# CrossCheck.py
# John Novak
# July 11, 2012

# This python script takes two model-data folders as arguments.
# The first model-data folder must contain a trained interactive emulator,
# the second need only contain a DataSummary file (which can be made using)
# setup.py.

# This script compares the output of the trained emulator at the points in the
# DataSummary file, then returns a measure of how well the emulator reproduced
# the model values.

import sys
import subprocess
import time
import os
import math as m


def main():
    if not (len(sys.argv) == 3):
        print "This program requires two directories as input. The first must contain a trained interactive emulator statefile, and the second must contain a DataSummary.dat file (which can be made using setup.py"
        print "Example: ./CrossCheck.py model-data/trained_emulator/ model-data/test-data/"
        exit(0)

    Statefile = sys.argv[1] + "/Emulator.statefile"
    TestPoints = open(sys.argv[2] + "/DataSummary.dat", 'r')

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
    obsvscalesfile_1.close()
    obsvscalesfile_2 = open(sys.argv[2]+"/ObservableScales.dat",'r')
    obsvscales_2 = []
    for line in obsvscalesfile_2:
        for i in obsvnames:
            if line.split(" ")[0] == i:
                obsvscales_2.append([float(line.split(" ")[1]),float(line[:-1].split(" ")[2])])
    obsvscalesfile_2.close()

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

    # Run through Ypoints and toss the observables we aren't using. And unscale everything
    for i in range(len(Ypoints)):
        temp = []
        count = 0
        for j in range(len(Ypoints[i])):
            if obsvnames[count] == obsvnames_1[j]:
                temp.append((Ypoints[i][j]*obsvscales_1[count][1])+obsvscales_1[count][0])
                count += 1
        Ypoints[i] = temp

    # Run Emulator on points and dump output to file
    os.system("interactive_emulator interactive_mode -q "+Statefile+" < EmulatedPoints.txt > EmulatedObservables.txt")

    # Read emulator output
    means = []
    errors = []
    longdiff = []
    diffsum = 0
    output1 = open("CrossCheckSummarylong.txt", 'w')
    output2 = open("CrossCheckSummarytermwise.txt", 'w')
    index1 = 0
    observablesfromfile = open("EmulatedObservables.txt",'r')
    for i in Xpoints:
        meanstemp=[]
        errorstemp=[]
        count = 0
        for j in range(num_obsv):
            if obsvnames[count] == obsvnames_1[j]:
                meanstemp.append((float(observablesfromfile.readline())*obsvscales_1[count][1])+obsvscales_1[count][0])
                errorstemp.append((float(observablesfromfile.readline())*obsvscales_1[count][1])+obsvscales_1[count][0])
                count += 1
        means.append(meanstemp)
        errors.append(errorstemp)
        diff = 0
        index2 = 0
        temp = []
        output2.write("Point #" + str(index1) + "\n")
        for j in means[-1]:
            diff += m.pow((j - Ypoints[index1][index2])/obsvscales_1[index2][1],2)/2
            temp.append(m.fabs(j - Ypoints[index1][index2])/obsvscales_1[index2][1])
            output2.write("Observable #" + str(index2) + " diff: " + str(m.fabs(j - Ypoints[index1][index2])/obsvscales_1[index2][1]) + " Model: " + str(Ypoints[index1][index2]) + " Emulator: " + str(j) + "\n")
            index2 += 1
        diffsum += m.sqrt(diff / num_obsv)
        longdiff.append(temp)
        output1.write("Point: " + str(i[:-1]) + "\nModel Values: " + str(Ypoints[index1]) + "\nEmultor Values: " + str(means[-1]) + " \nave diff: " + str(m.sqrt(diff / num_obsv)) + "\n")
        index1 += 1
    output1.close()
    output2.close()
    observablesfromfile.close()

    output = open("CrossCheckSummaryshort.txt", 'w')
    output.write("Average difference: " + str(diffsum / len(Ypoints)) + "\n")
    for i in range(len(longdiff[0])):
        termdiff = 0
        for j in range(len(longdiff)):
            termdiff += longdiff[j][i]
        termdiff = termdiff / len(longdiff)
        output.write("Observable #" + str(i) + " deviated on average " + str(termdiff) + "\n")
    output.close()

if __name__ == '__main__':
    main()
