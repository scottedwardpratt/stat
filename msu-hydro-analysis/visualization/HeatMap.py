#!/usr/bin/python

# HeatMap.py
# This is a general plotting script which takes a text file and generates a heat map from it using gnuplot
# WARNING: This code is not nearly as robust as some of the other plotting scripts, so it probably wont take
# much creativity to crash.

# John Novak
# July 1, 2012

import os
import sys
import subprocess

def main():
    arguments=['filename','xlabel','ylabel','title']
    state = 0
    for flag in sys.argv:
        if flag == "-xl":
            state = 1
        elif flag == "-yl":
            state = 2
        elif flag == "-t":
            state = 3
        else:
            arguments[state]=flag

    command = "set term postscript enhanced color; set view map; set palette rgbformula -23,-28,-3; set pointsize; set points -1;\n"

    command+='set xlabel "'+arguments[1]+'";\n'
    command+='set ylabel "'+arguments[2]+'";\n'
    command+='set title "'+arguments[3]+'";\n'
    command+='set output "'+arguments[0].split('.')[0]+'.ps";\n'
    command+='splot "'+arguments[0]+'" matrix with image;\n'

    plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
    plot.communicate(command)

if __name__ == '__main__':
    main()
