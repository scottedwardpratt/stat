#!/usr/bin/env python
"""
Dice Model
Copyright 2013, The University of North Carolina at Chapel Hill.

This software was written in 2012 by Hal Canary <cs.unc.edu/~hal/>
while working for the MADAI project <http://madai.us/>.

See <https://madai-public.cs.unc.edu/software/license/> for license
information.

This is intended as an example of using the dice model to create a
emulator.

	$ run_interactive_emulator estimate_thetas \
	    train.dat model_snapshot_file.txt
"""

import sys, os, math, numpy, subprocess
import dice_model_example

output_file = 'train.dat'
data_file_seperator = ' '
model = dice_model_example.dice_model2

if __name__ == '__main__':
	# To specify a particular set of input points, use a input file.
	# data_file = 'test_points.dat'
	# X = numpy.genfromtxt(data_file, delimiter=data_file_seperator)

	sub = subprocess.Popen(
		['./latin-hyper-sample','6','60'], stdout=subprocess.PIPE)

	""" X is a matrix of points in parameter space """
	X = numpy.genfromtxt(sub.stdout, delimiter=data_file_seperator)

	""" Y is a matrix of points in the model's output space """
	Y = [ model(list(x)) for x in X ]

	o = open(output_file,'w') # o = sys.stdout

	"""
	Creates an input file suitable for interactive_emulator's
	estimate_thetas function.
	"""
	print >>o, len(Y[0])
	print >>o, len(X[0])
	print >>o, len(X)
	for x in X:
		print >>o, ' '.join(map(str,x))
	for y in Y:
		print >>o, ' '.join(map(str,y))
