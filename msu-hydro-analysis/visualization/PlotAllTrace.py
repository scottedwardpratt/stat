#!/usr/bin/python -tt

# This script is to be run inside a trace folder (model-data/*/trace/)
# It creates trace density plots from the trace which is the aggregate of all
# the traces in the various default folders.

import sys
import os
import string
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def main():
    print("And we are off"),
    os.system("ls -1 > temp") # a list of all the directorys

    count=0

    data = np.array([[[[0]*51]*51]*6]*6)

    # First, open the list directories
    list1=open("temp","r")
    for line1 in list1:
        line1 = line1[:-1]
        if line1[:3]=="def": # if its a trace, use it
            #print line1
            os.system("ls -1 "+line1+"/ > temp"+line1) # make a list of the outputs in it
            list2=open("temp"+line1) # open the list
            for line2 in list2: # for each file in the directory
                if line2[:4] == "outp": # if the file is an output use it
                    #print line2
                    trace=open(line1+"/"+line2[:-1],"r") # open the file
                    names=trace.readline() # the first line is the param names
                    names=names[1:-1].split(",")
                    vals=map(int,map(lambda x: x*50,map(float,names[8:])))
                    #print names[8:]
                    #print vals
                    names = names[1:7]
                    for i in range(6):
                        for j in range(i+1,6):
                            if vals[i] < 50 and vals[j] < 50:
                                data[i,j,vals[i],vals[j]] += 1
                            elif vals[i] >= 50:
                                print line1+"/"+line2[:-1],":",names[i],"outside of range:",vals[i]/50
                    for line3 in trace: # for each line
                        vals=map(int,map(lambda x: x*50,map(float,line3.split(",")[1:]))) # First, turn the strings to nums, then mult by 50, then take the floor
                        #print vals
                        for i in range(6):
                            for j in range(i+1,6):
                                if vals[i] < 50 and vals[j] < 50:
                                    data[i,j,vals[i],vals[j]] += 1
                                elif vals[i] >= 50:
                                    print line1+"/"+line2[:-1],":",names[i],"outside of range:",vals[i]/50
                    trace.close()
                    count += 1
                    if count == 5000:
                        writeout( names, data )
                        data = np.array([[[[0]*51]*51]*6]*6)
            list2.close()
            os.system("rm temp"+line1)
    list1.close()
    output.close()
    os.system("rm temp")

    # Now, we have all the data loaded, print it all out
    writeout( names, data )

def writeout( names, data ):
    print("."),
    for i in range(6):
        for j in range(i+1,5):
            fileout = open(names[i]+names[j]+".txt",'w')
            for k in range(50):
                for l in range(50):
                    fileout.write(str(k)+" "+str(l)+" "+str(data[i,j,k,l].sum())+"\n")
            fileout.close()

if __name__ == '__main__':
    main()
