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
            line=line[:-1]
            if line[-1] == '\r':
                line=line[:-1]
            if line[-1] == ' ':
                line=line[:-1]
            obsvnames.append(line)
    print "Observables:",obsvnames
    obsvlist.close()

    # Read in the names and ranges of the parameters
    paramlist=open(sys.argv[1]+"/ranges.dat","r")
    paramnames=[]
    ranges=[]
    for line in paramlist:
        #print line.split(" ")
        line=line[:-1]
        if line[-1] == '\r':
            line=line[:-1]
        if line[-1] == ' ':
            line=line[:-1]
        paramnames.append(line.split(" ")[1])
        ranges.append([float(line.split(" ")[2]),float(line.split(" ")[3])])
    print "Paramteres:",paramnames
    print "Parameter ranges:",ranges
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
                line=line[:-1]
                #print line
                for i in paramnames:
                      if line.split(" ")[1] == i:
                          temp.append(line.split(" ")[2])
            paramvals.append(temp)
            paramfile.close()
    print "A sample set of parameter values:",paramvals[0]
    dirlist.close()

    # Read in the observables values
    os.system("ls -1 "+sys.argv[1]+"/model_results > dirlist")
    dirlist=open("dirlist","r")
    for run in dirlist:
        if run[:3]=="run":
            temp=[]
            obsvfile=open(sys.argv[1]+"/model_results/"+run[:-1]+"/results.dat")
            #print "Reading in",sys.argv[1]+"/model_results/"+run[:-1]+"/results.dat"
            for line in obsvfile:
                line=line[:-1]
                #print line
                for i in obsvnames:
                      if line.split(" ")[1] == i:
                          #print i,line.split(" ")[2],line
                          temp.append(line.split(" ")[2:])
            obsvvals.append(temp)
            obsvfile.close()
    if len(paramvals) != len(obsvvals):
        print "You have a differnt number of parameter sets than observable sets! ABORT! ABORT! OH, THE HUMANITY!!!"
        print "Number of parameter sets:",len(paramvals)
        print "Number of observable sets:",len(obsvvals)
        exit(1)
    print "A sample set of observable values:",obsvvals[0]
    dirlist.close()
    os.system("rm dirlist")

    # Determine the mean and standard deviation of the observables
    means=[0]*len(obsvvals[0])
    stdv=[0]*len(obsvvals[0])

    #print "Means and sigmas:",means,stdv
    for i in obsvvals:
        for j in range(len(i)):
            means[j]+=float(i[j][0]);
    for i in range(len(means)):
        means[i]=means[i]/len(obsvvals)
    for i in obsvvals:
        for j in range(len(i)):
            stdv[j]+=float(i[j][1]);
    for i in range(len(stdv)):
        stdv[i]=stdv[i]/len(obsvvals)

    print "Means and sigmas:",means,stdv
    outputname = sys.argv[1]+"/DataSummary.dat"
    output=open(outputname,"w")
    #output.write("#Number of outputs\n")
    output.write(str(len(obsvvals[0]))+"\n")
    #output.write("#Number of inputs\n")
    output.write(str(len(paramvals[0]))+"\n")
    #output.write("#Number of training points\n")
    output.write(str(len(obsvvals))+"\n")
    #output.write("#")
    #for i in paramnames:
        #output.write(i+" ")
    #output.write("\n")
    for i in paramvals:
        temp=""
        c=0
        for j in i:
            k = (float(j)-ranges[c][0])/(ranges[c][1]-ranges[c][0])
            c+=1
            temp+=str(k)+" "
        output.write(temp+"\n")
    #output.write("#")
    #for i in obsvnames:
        #output.write(i+" ")
    #output.write("\n")
    for i in obsvvals:
        temp=""
        c=0
        for j in i:
            k = (float(j[0])-means[c])/(stdv[c])
            #spectra[c][int(float(j)] += 1
            c+=1
            temp+=str(k)+" "
        output.write(temp+"\n")
    output.close()
    print "Output written to",outputname

    outputname=sys.argv[1]+"/ObservableScales.dat"
    output=open(outputname,"w")
    output.write("#These are the values which the observables were scaled with.\n")
    output.write("#They were scaled as: y'=(y-<y>)/y_error.\n")
    output.write("#The values below are: Name, <y>, y_error\n")
    #print len(obsvnames)
    #print obsvnames[0],means[0],stdv[0]
    for i in range(len(obsvnames)):
        output.write(str(obsvnames[i])+" "+str(means[i])+" "+str(stdv[i])+"\n")
    output.close()
    print "Observable's means and errors printed to:",outputname

    if not os.path.exists(sys.argv[1]+"/exp_data/results_unscaled.dat"):
        os.system("mv "+sys.argv[1]+"/exp_data/results.dat "+sys.argv[1]+"/exp_data/results_unscaled.dat")
        print "Moved experimental data to",sys.argv[1]+"/exp_data/results_unscaled.dat"

    outputname=sys.argv[1]+"/exp_data/results.dat"
    exp_data=open(sys.argv[1]+"/exp_data/results_unscaled.dat","r")
    output=open(outputname,"w")
    j=0
    for line in exp_data:
        for i in obsvnames:
            if line.split(" ")[1] == i:
                #print line.split(" ")[1],i,j,float(line.split(" ")[2]),means[j],stdv[j]
                output.write("double "+i+" "+str((float(line.split(" ")[2])-means[j])/stdv[j])+" -999\n")
                j+=1
    output.close()
    print "Scaled experimental data printed to ",outputname

if __name__ == '__main__':
    main()
