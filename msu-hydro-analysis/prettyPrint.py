#!/usr/bin/python

# prettyPrint.py
# John F Novak (JFN)
# Saturday, November 10, 2012, 11:29am

# This python script takes up to two model-data folders as arguments.
# The first model-data folder must contain a trained interactive
# emulator and a pcanames.data file, the second need only contain a
# pcanames.dat file. If no second folder is provided, the code will
# look in the first for pcanames.dat

# This script runs the emulator at a user given point, then unscales the
# observables returned by the emulator. If a second directory is provided
# the code will only return the observables listed on both pcanames.dat

import sys
import subprocess
import time
import os
import math as m


def main():
    if (len(sys.argv) == 1):
        print "This program requires at least one directory as input. The first must contain a trained interactive emulator statefile and a pcanames.dat file. If a second directory is provided, it must contain a pcanames.dat file aswell."
        print "Example: ./pretyPrint.py model-data/trained_emulator/ model-data/test-data/"
        exit(0)

    Statefile = sys.argv[1] + "/Emulator.statefile"

    if len(sys.argv) == 3:
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
        num_obsv = len(obsvnames_2)

        # Compare observables lists and figure out what we are going to use
        same = True
        if len(obsvnames_1) != len(obsvnames_2):
            same = False
        else:
            for i in range(len(obsvnames_1)):
                if obsvnames_1[i] != obsvnames_2[i]:
                    same = False
        if not same:
            #print "The observables used in the emualtor and in the data are not the same:"
            #print len(obsvnames_1),obsvnames_1
            #print len(obsvnames_2),obsvnames_2
            #print "We will use the intersection of the two lists:"
            obsvnames = [val for val in obsvnames_1 if val in obsvnames_2]
        else:
            obsvnames = obsvnames_1
            #print "Observables:"
        #print len(obsvnames),obsvnames

    else:
        # Read in the names of the observables used in the emulator
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
        num_obsv = len(obsvnames)
        obsvnames_1 = obsvnames
        obsvlist.close()

    # Read in the observables scales
    obsvscalesfile = open(sys.argv[1]+"/ObservableScales.dat",'r')
    obsvscales = []
    for line in obsvscalesfile:
        for i in obsvnames:
            if line.split(" ")[0] == i:
                obsvscales.append([float(line.split(" ")[1]),float(line[:-1].split(" ")[2])])
    obsvscalesfile.close()

    point = str(input("Enter a point (-1 to exit): "))
    while point != "-1":

        # Make file containing point to be passed to emulator
        emupoint = open("EmulatedPoints.txt",'w')
        emupoint.write(point+"\n")
        emupoint.close()

        # Run Emulator on points and dump output to file
        os.system("interactive_emulator interactive_mode -q "+Statefile+" < EmulatedPoints.txt > EmulatedObservables.txt")

        # Read emulator output
        observablesfromfile = open("EmulatedObservables.txt",'r')
        count = 0
        for j in range(num_obsv):
            if obsvnames[count] == obsvnames_1[j]:
                mean=((float(observablesfromfile.readline())*obsvscales[count][1])+obsvscales[count][0])
                #error=float(observablesfromfile.readline())
                #error=((float(observablesfromfile.readline())*obsvscales[count][1])+obsvscales[count][0])
                error=obsvscales[count][1]
                print obsvnames[count],mean,error
                count += 1
        observablesfromfile.close()

        # Do it again
        point = str(input("Enter a point (-1 to exit): "))

if __name__ == '__main__':
    main()
