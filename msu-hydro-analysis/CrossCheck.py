#!/usr/bin/python

# CrossCheck.py
# John Novak
# July 11, 2012

# This python script takes two model-data folders as arguments.
# The first model-data folder must contain a trained interactive emulator,
# the second need only contain a DataSummary file (which can be made using)
# setup.py.

# This script compares the output of the trained emulator at the points in the
# DataSummary file, then returns a measure of how well the emulator reproduced
# the model values.

import sys
import subprocess
import time
import math as m


def main():
    if not (len(sys.argv) == 3):
        print "This program requires two directories as input. The first must contain a trained interactive emulator statefile, and the second must contain a DataSummary.dat file (which can be made using setup.py"
        print "Example: ./CrossCheck.py model-data/trained_emulator/ model-data/test-data/"
        exit(0)

    Statefile = sys.argv[1] + "/Emulator.statefile"
    TestPoints = open(sys.argv[2] + "/DataSummary.dat", 'r')

    Emulator = subprocess.Popen(['interactive_emulator', 'interactive_mode', '-q', Statefile], stdin=subprocess.PIPE, shell=False, stdout=subprocess.PIPE)

    time.sleep(2)

    num_params = int(TestPoints.readline())
    num_obsv = int(TestPoints.readline())
    num_points = int(TestPoints.readline())
    Xpoints = []
    Ypoints = []
    for i in range(num_points):
        line = TestPoints.readline()
        Xpoints.append(line)  # Here we are taking the points as strings, just as they are because that is how the emulator will want them
        if len(Xpoints[-1][:-2].split(" ")) != num_obsv:
            print "Parsing error. The number of observables expected in " + sys.argv[2] + "/DataSummary.dat is not the same as read in. Point", i
            print "Expected", num_obsv, "observables. Got", Xpoints[-1][:-2].split(" "), "Exiting"
            exit(1)
    for i in range(num_points):
        line = TestPoints.readline()
        Ypoints.append(map(float, line[:-2].split(" ")))  # We turn these into an array of nums so we can do math with them
        if len(Ypoints[-1]) != num_params:
            print "Parsing error. The number of parameters expected in " + sys.argv[2] + "/DataSummary.dat is not the same as read in. Point", i, "Exiting"
            exit(1)

    means = []
    errors = []
    longdiff = []
    diffsum = 0
    output1 = open("CrossCheckSummarylong.txt", 'w')
    output2 = open("CrossCheckSummarytermwise.txt", 'w')
    index1 = 0
    for i in Xpoints:
        Emulator.stdin.write(i)
        time.sleep(1)
        meanstemp = []
        errorstemp = []
        for j in range(num_obsv):
            meanstemp.append(float(Emulator.stdout.readline()))
            errorstemp.append(float(Emulator.stdout.readline()))
        means.append(meanstemp)
        errors.append(errorstemp)
        diff = 0
        index2 = 0
        temp = []
        output2.write("Point #" + str(index1) + "\n")
        for j in means[-1]:
            diff += (j - Ypoints[index1][index2]) * (j - Ypoints[index1][index2])
            temp.append(m.sqrt((j - Ypoints[index1][index2]) * (j - Ypoints[index1][index2])))
            output2.write("Observable #" + str(index2) + " diff: " + str(m.sqrt((j - Ypoints[index1][index2]) * (j - Ypoints[index1][index2]))) + " Model: " + str(Ypoints[index1][index2]) + " Emulator: " + str(j) + "\n")
            index2 += 1
        diffsum += m.sqrt(diff / num_obsv)
        longdiff.append(temp)
        output1.write("Point: " + str(i[:-1]) + "\nModel Values: " + str(Ypoints[index1]) + "\nEmultor Values: " + str(means[-1]) + " \nave diff: " + str(m.sqrt(diff / num_obsv)) + "\n")
        index1 += 1
        #if index1 == 1:
            #break
    output1.close()
    output2.close()

    output = open("CrossCheckSummaryshort.txt", 'w')
    output.write("Average difference: " + str(diffsum / len(Ypoints)) + "\n")
    for i in range(len(longdiff[0])):
        termdiff = 0
        for j in range(len(longdiff)):
            termdiff += longdiff[j][i]
        termdiff = termdiff / len(longdiff)
        output.write("Observable #" + str(i) + " deviated on average " + str(termdiff) + "\n")
    output.close()

    Emulator.terminate()

if __name__ == '__main__':
    main()
