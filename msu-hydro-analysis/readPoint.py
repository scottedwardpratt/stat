#!/usr/bin/python

# readPoint.py
# John Novak
# Saturday, November 10, 2012, 10:38am

# This code reads data points from DataSummary.dat and unscales them

import os
import sys
import math as m

def main():
    if not (len(sys.argv) == 2):
        print "This program requires one directories as input. The first must contain an emulator statefile (named Emulator.statefile), and the second must contain everything else: a data summary file (named DataSummary.dat), a folder containing the experimental data (named exp_data), and a pcanames file (and it should be the same as the one used for the emulator. These two directories can be the same"
        print "Example: ./Prepare_for_Scott.py model-data/datafile/"
        exit(0)

    TestPoints = open(sys.argv[1] + "/DataSummary.dat", 'r')

    # Read in the names of the observables used in the data set
    obsvlist=open(sys.argv[1]+"/pcanames.dat","r")
    obsvnames=[]
    for line in obsvlist:
        if line[0] != "#":
            line=line[:-1]
            if line[-1] == '\r':
                line=line[:-1]
            if line[-1] == ' ':
                line=line[:-1]
            obsvnames.append(line)
    obsvlist.close()

    # Read in the observables scales
    obsvscalesfile = open(sys.argv[1]+"/ObservableScales.dat",'r')
    obsvscales = []
    for line in obsvscalesfile:
        for i in obsvnames:
            if line.split(" ")[0] == i:
                obsvscales.append([float(line.split(" ")[1]),float(line[:-1].split(" ")[2])])

    # Read in the parameters scales
    rangefile = open(sys.argv[1]+"/ranges.dat",'r')
    ranges = []
    for line in rangefile:
        ranges.append([float(line.split(" ")[2]),float(line[:-1].split(" ")[3])])

    # Reading in data from DataSummary.dat
    num_obsv = int(TestPoints.readline())
    num_params = int(TestPoints.readline())
    num_points = int(TestPoints.readline())
    Xpoints = []
    Ypoints = []
    for i in range(num_points):
        line = TestPoints.readline()
        Xpoints.append(map(float, line[:-2].split(" ")))
        if len(Xpoints[-1]) != num_params:
            print "Parsing error. The number of observables stated in " + sys.argv[2] + "/DataSummary.dat is not the same as the number read in. Problem with point", i
            print "Expected", num_params, "observables. Got", len(Xpoints[-1][:-2].split(" ")),Xpoints[-1][:-2].split(" "), "Exiting"
            exit(1)
    for i in range(num_points):
        line = TestPoints.readline()
        Ypoints.append(map(float, line[:-2].split(" ")))
        if len(Ypoints[-1]) != num_obsv:
            print "Parsing error. The number of parameters stated in " + sys.argv[2] + "/DataSummary.dat is not the same as the number read in. Problem with point", i, "Exiting"
            print "Expected", num_obsv, "observables. Got", len(Ypoints[-1]),Ypoints[-1], "Exiting"
            exit(1)

    # Run through Xpoints and unscale everything
    Xpoints_scaled = Xpoints
    Xpoints_unscaled = []
    for i in range(len(Xpoints)):
        temp = []
        count = 0
        for j in range(len(Xpoints[i])):
            temp.append((Xpoints[i][j]*ranges[count][1])+ranges[count][0])
            count += 1
        Xpoints_unscaled.append(temp)

    # Run through Ypoints and unscale everything
    for i in range(len(Ypoints)):
        temp = []
        count = 0
        for j in range(len(Ypoints[i])):
            temp.append((Ypoints[i][j]*obsvscales[count][1])+obsvscales[count][0])
            count += 1
        Ypoints[i] = temp

    # Interact with the user so that they can select a point.
    point = int(input("Enter a point between 1 and "+str(len(Xpoints))+". Enter -1 to exit "))
    while(point != -1):
        print "scaled:"," ".join(map(str, Xpoints_scaled[point-1]))
        print "unscaled:"," ".join(map(str, Xpoints_unscaled[point-1]))
        #print " ".join(map(str, Ypoints[point]))
        for i in range(len(Ypoints[point-1])):
            print obsvnames[i], Ypoints[point-1][i]
        point = int(input("Enter a point between 1 and "+str(len(Xpoints))+". Enter -1 to exit "))

if __name__ == '__main__':
    main()
