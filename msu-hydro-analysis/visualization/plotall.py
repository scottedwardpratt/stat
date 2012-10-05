#!/usr/bin/python -tt

# John Novak
# June 22, 2012

# This code is a more robust version of the old plotall gnuplot script

# This script is meant to be run from within the DensityPlots folder (model-data/*/DensityPlots)
# This looks at the ranges.dat file, anticipates what files should be in DensityPlots, then
# creates gnuplot visualizations of them.

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

    OneD=[]
    for i in range(len(params)):
        for j in range(i+1,len(params)):
            command+='set ylabel "'+params[i].replace("_"," ")+'";\n'
            command+='set xlabel "'+params[j].replace("_"," ")+'";\n'
            command+='set output "'+params[i]+"_"+params[j]+'.ps";\n'
            command+='splot "'+params[i]+"_"+params[j]+'.txt" matrix with image;\n'
            matrix=open(params[i]+"_"+params[j]+'.txt','r')
            data=[]
            for line in matrix:
                data.append(line[:-1].split(" "))
            matrix.close()
            output=open(params[j]+"_"+params[i]+'.txt','w')
            for k in range(len(data[0])):
                for l in range(len(data)):
                    output.write(data[l][k]+" ")
                output.write("\n")
            output.close()
        temp=[]
        temp2=[0]*100
        if i+1 != len(params):
            matrix2=open(params[i+1]+"_"+params[i]+'.txt','r')
        if i+1 == len(params):
            matrix2=open(params[i-1]+"_"+params[i]+'.txt','r')
        for line in matrix2:
            if line[:-2].split(" ")[0] != '':
                temp.append(map(float,line[:-2].split(" ")))
        matrix2.close()
        for thing in temp:
            for k in range(len(thing)):
                temp2[k]+=thing[k]
        OneD.append(temp2)
                
    command+='set output "multi.ps"; set multiplot layout 6,6 rowsfirst; set bmargin 0.001; set tmargin 0.001; set lmargin 0.001; set rmargin 0.001; unset key; unset colorbox; unset xtics; unset ytics; set ylabel font "Courier,6"; set xlabel font "Courier,6" offset 0,2;\n'

    for i in range(len(params)):
        for j in range(i+1,len(params)):
            command+='set origin '+str(1-(j-i)*.166666666)+','+str(1-((i+1)*.166666666))+';\n'
            command+='set ylabel "'+params[i].replace("_"," ")+'";\n'
            command+='set xlabel "'+params[len(params)-j+i].replace("_"," ")+'";\n'
            command+='splot "'+params[i]+"_"+params[len(params)-j+i]+'.txt" matrix with image;\n'

    for i in range(len(params)):
        for j in range(i+1,len(params)):
            command+='set origin '+str(i*.166666666)+','+str((j-i-1)*.166666666)+';\n'
            command+='set ylabel "'+params[len(params)-j+i].replace("_"," ")+'";\n'
            command+='set xlabel "'+params[i].replace("_"," ")+'";\n'
            command+='splot "'+params[len(params)-j+i]+"_"+params[i]+'.txt" matrix with image;\n'
    
    command+='set bmargin 1.5; set tmargin 1.5; set lmargin 2; set rmargin 3; unset key; unset colorbox; unset xtics; unset ytics; set ylabel font "Courier,6"; set xlabel font "Courier,6" offset 0,0; set style data points;\n'

    for i in range(len(params)):
        command+='set origin '+str((i*.166666666)-.01)+','+str((1-((i+1)*.166666666))-.02)+';\n'
        command+='set ylabel "Density";\n'
        command+='set xlabel "'+params[i].replace("_"," ")+'";\n'
        command+="plot '-' pt 7;\n" 
        for j in OneD[i]:
            command+=str(j)+';\n'
        command+='e;\n'

    #talking to gnuplot
    plot = subprocess.Popen(['gnuplot'], stdin=subprocess.PIPE)
    plot.communicate(command)

if __name__ == '__main__':
    main()
