#!/usr/bin/python -tt

# This script is meant to be run in a sub-folder of trace (mode-data/*/trace/default*/)
# This takes all the output*.dat files and compiles an aggreate trace file named fulltrace.csv

import sys
import os
import string
import numpy
import matplotlib
import matplotlib.pyplot as plt

def main():
  os.system("ls -1 > temp")

  count=0

  output=open("fulltrace.csv","w")

  list1=open("temp","r")
  for line1 in list1:
    line1 = line1[:-1]
    if line1[:3]=="out":
      trace=open(line1,"r")
      trash=trace.readline()
      for line3 in trace:
        vals=line3.split(",")[1:]
        count+=1
        toprint=str(count)
        for entry in vals:
          toprint=toprint+","+entry
        output.write(toprint)
      trace.close()
  list1.close()
  output.close()
  os.system("rm temp")

if __name__ == '__main__':
  main()
