#!/usr/bin/python -u

# Emulator Emulator (used for testing emulator interface)
# MADAI Model Statistical Tools
# Copyright 2011-2012, The University of North Carolina at Chapel Hill.
#
# This software was written in 2011-2012 by 
# 	Hal Canary <hal AT cs.unc.edu>
# while working for the MADAI project <http://madai.us/>.
#
# See copyright.txt for more information.
#
# Be sure to call with '-u' flag.


import sys
from math import *
def tok(readable):
	while True:
		line = readable.readline()
		if line == '':
			break
		for token in line.split():
			yield token
def f(x,y):
	scale = 0.2 * exp(-16 * ((x - 0.5)**2 + (y - 0.5)**2))
	return scale * (3.0 + sin(25 * x) + sin(25 * y))
o = sys.stdout
it = tok(sys.stdin)
o.write('2\nparam_0\nparam_1\n2\nmean_0\nvariance_0\n')
o.flush()
while True:
	try:
		x = float(it.next())
		y = float(it.next())
		o.write('%r\n%r\n' % (f(x,y), 0.0))
		o.flush()
	except (StopIteration,):
		break
