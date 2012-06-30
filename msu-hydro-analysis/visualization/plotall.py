#!/usr/bin/python -tt

# John Novak
# June 22, 2012

# This code is a more robust version of the old plotall gnuplot script

import sys
import os
import subprocess

def main():
    command = "set term postscript enhanced color; set view map; set palette rgbformula -23,-28,-3; set pointsize; set points -1;\n"
    #files=subprocess.call(["ls", "-l"])
    ranges=open("../ranges.dat","r")
    params=[]
    for line in ranges:
        params.append(line.split()[1])

    for i in range(len(params)):
        for j in range(i+1,len(params)):
            command+='set ylabel "'+params[i].replace("_"," ")+'";\n'
            command+='set xlabel "'+params[j].replace("_"," ")+'";\n'
            command+='set output "'+params[i]+"_"+params[j]+'.ps";\n'
            command+='splot "'+params[i]+"_"+params[j]+'.txt" matrix with image;\n'

    command+='set output "multi.ps"; set multiplot layout 4,4 rowsfirst; set bmargin 0.001; set tmargin 0.001; set lmargin 0.001; set rmargin 0.001; unset key; unset colorbox; unset xtics; unset ytics; set ylabel font "Courier,11"; set xlabel font "Courier,11" offset 0,2;\n'

    for i in range(len(params)):
        for j in range(i+1,len(params)):
            command+='set ylabel "'+params[i].replace("_"," ")+'";\n'
            command+='set xlabel "'+params[j].replace("_"," ")+'";\n'
            command+='splot "'+params[i]+"_"+params[j]+'.txt" matrix with image;\n'

    #talking to gnuplot
    plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
    plot.communicate(command)

if __name__ == '__main__':
    main()
