#!/usr/bin/python -tt

# John Novak
# June 22, 2012

# This program prepares a .dat file to be used to train the Interactive Emulator

import sys
import os
import math as m

def main():
    if len(sys.argv) != 2:
        print "This program prepares a .dat file to be used to train the Interactive Emulator"
        print "It names the output DataSummary.dat and puts it in the model-data folder"
        print "This is meant to be run from msu-hyrdro-analysis"
        print "usage: python setup.py <model-data folder>"
        print "   or: ./setup.py <model-data folder>"
        print "example: python setup.py model-data/June2012"
        exit(1)

    # Read in the names of the observables
    obsvlist=open(sys.argv[1]+"/pcanames.dat","r")
    obsvnames=[]
    for line in obsvlist:
        if line[0] != "#":
            obsvnames.append(line[:-1])
    print obsvnames
    obsvlist.close()

    # Read in the names and ranges of the parameters
    paramlist=open(sys.argv[1]+"/ranges.dat","r")
    paramnames=[]
    ranges=[]
    for line in paramlist:
        paramnames.append(line[:-1].split(" ")[1])
        ranges.append([float(line[:-2].split(" ")[2]),float(line[:-2].split(" ")[3])])
    print paramnames
    print ranges
    paramlist.close()

    paramvals=[]
    obsvvals=[]

    # Read in the parameter values
    os.system("ls -1 "+sys.argv[1]+"/parameters > dirlist")
    dirlist=open("dirlist","r")
    for run in dirlist:
        if run[:3]=="run":
            temp=[]
            paramfile=open(sys.argv[1]+"/parameters/"+run[:-1]+"/stats.param")
            for line in paramfile:
                for i in paramnames:
                      if line.split(" ")[1] == i:
                          temp.append(line.split(" ")[2][:-1])
            paramvals.append(temp)
            paramfile.close()
    print paramvals[0]
    dirlist.close()

    # Read in the observables values
    os.system("ls -1 "+sys.argv[1]+"/model_results > dirlist")
    dirlist=open("dirlist","r")
    for run in dirlist:
        if run[:3]=="run":
            temp=[]
            obsvfile=open(sys.argv[1]+"/model_results/"+run[:-1]+"/results.dat")
            for line in obsvfile:
                for i in obsvnames:
                      if line.split(" ")[1] == i:
                          temp.append(line.split(" ")[2][:-1])
            obsvvals.append(temp)
            obsvfile.close()
    dirlist.close()
    os.system("rm dirlist")

    # Determine the mean and standard deviation of the observables
    means=[0]*len(obsvvals[0])
    stdv=[0]*len(obsvvals[0])
    for i in obsvvals:
        for j in range(len(i)):
            means[j]+=float(i[j]);
    for i in range(len(means)):
        means[i]=means[i]/len(obsvvals)
    for i in obsvvals:
        for j in range(len(i)):
            stdv[j]+=(float(i[j])-means[j])*(float(i[j])-means[j]);
    for i in range(len(stdv)):
        stdv[i]=m.sqrt(stdv[i]/len(obsvvals))

    # Write out the .dat file
    #if sys.argv[1][-1] == "/":
    #    outputname = sys.argv[1][:-1].split("/")[-1]+".dat"
    #else:
    #    outputname = sys.argv[1].split("/")[-1]+".dat"
    outputname = sys.argv[1]+"/DataSummary.dat"
    output=open(outputname,"w")
    output.write(str(len(obsvvals[0]))+"\n")
    output.write(str(len(paramvals[0]))+"\n")
    output.write(str(len(obsvvals))+"\n")
    for i in paramvals:
        temp=""
        c=0
        for j in i:
            k = (float(j)-ranges[c][0])/(ranges[c][1]-ranges[c][0])
            c+=1
            temp+=str(k)+" "
        output.write(temp+"\n")
    for i in obsvvals:
        temp=""
        c=0
        for j in i:
            k = (float(j)-means[c])/(stdv[c])
            c+=1
            temp+=str(k)+" "
        output.write(temp+"\n")
    output.close()
    print "Output written to",outputname

if __name__ == '__main__':
    main()
