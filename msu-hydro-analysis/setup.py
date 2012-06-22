#!/usr/bin/python -tt

# John Novak
# June 22, 2012

# This program prepares a .dat file to be used to train the Interactive Emulator

import sys
import os

def main():
    if len(sys.argv) != 2:
        print "This program prepares a .dat file to be used to train the Interactive Emulator"
        print "This is meant to be run from msu-hyrdro-analysis"
        print "usage: python setup.py <model-data folder>"
        print "example: python setup.py model-data/June2012"
        exit(1)

    obsvlist=open("data-prep/standardlong.dat","r")
    obsvnames=[]
    for line in obsvlist:
        obsvnames.append(line[:-1])
    print obsvnames
    obsvlist.close()

    paramlist=open(sys.argv[1]+"/ranges.dat","r")
    paramnames=[]
    for line in paramlist:
        paramnames.append(line[:-1].split(" ")[1])
    print paramnames
    paramlist.close()

    paramvals=[]
    obsvvals=[]

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

    output=open(sys.argv[1].split("/")[-1]+".dat","w")
    output.write(str(len(obsvvals[0]))+"\n")
    output.write(str(len(paramvals[0]))+"\n")
    output.write(str(len(obsvvals))+"\n")
    for i in paramvals:
        temp=""
        for j in i:
            temp+=j+" "
        output.write(temp+"\n")
    for i in obsvvals:
        temp=""
        for j in i:
            temp+=j+" "
        output.write(temp+"\n")
    output.close()
    print "Output written to",sys.argv[1].split("/")[-1]+".dat"

if __name__ == '__main__':
    main()
