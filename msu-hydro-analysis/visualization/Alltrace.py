#!/usr/bin/python -tt

# This script is to be run inside a trace folder (model-data/*/trace/)
# It creates one file named fulltrace.txt which is an aggregate of all
# the traces in the various default folders.

import sys
import os
import string
import numpy
import matplotlib
import matplotlib.pyplot as plt

def main():
  os.system("ls -1 > temp")

  count=0

  output=open("fulltrace.txt","w")

  list1=open("temp","r")
  for line1 in list1:
    line1 = line1[:-1]
    if line1[:3]=="def":
      #print line1
      os.system("ls -1 "+line1+"/ > temp"+line1)
      list2=open("temp"+line1)
      for line2 in list2:
        if line2[:4] == "outp":
          #print line2
          trace=open(line1+"/"+line2[:-1],"r")
          trash=trace.readline()
          for line3 in trace:
            vals=line3.split(",")[1:]
            count+=1
            toprint=str(count)
            for entry in vals:
              toprint=toprint+","+entry
            output.write(toprint)
          trace.close()
      list2.close()
      os.system("rm temp"+line1)
  list1.close()
  output.close()
  os.system("rm temp")

if __name__ == '__main__':
  main()
