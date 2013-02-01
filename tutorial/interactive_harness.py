#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test and use interactive_emulator from the MADAI Emulator
Copyright 2013, The University of North Carolina at Chapel Hill.

This software was written in 2012 by Hal Canary <cs.unc.edu/~hal/>
while working for the MADAI project <http://madai.us/>.

See <https://madai-public.cs.unc.edu/software/license/> for license
information.

1.	Create a model.  This is a Python function that takes a list or
  	tuple as an argument and returns an list of outputs

2.	Generate a latin hypercube of training points using either of the
 	maximum_latin_hypercube() or random_latin_hypercube() functions.

3.	Evaluate the model at the training points.

5.	Train the emulator (interactive_emulator estimate_thetas).

6.	Test emulator at each of the training points.

7.	 Test

"""

import subprocess
import os
import numpy
import tempfile
from dice_model import dice_model

executable = '/home/hal/local/bin/run_interactive_emulator'

class Emulator:
	"""
	The Emulator class wraps an external MADAIEmulator in a function
	"""
	def __init__(self, executable, model_snapshot_file):
		self.model_snapshot_file = model_snapshot_file
		self.process = subprocess.Popen(
			[ executable, 'interactive_mode', model_snapshot_file ],
			stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		l = self.process.stdout.readline()
		self.nparams = int(l)
		self.parameter_names = [ 			self.process.stdout.readline().strip() 
			for i in xrange(self.nparams) ]
		self.nreturns = int(self.process.stdout.readline())
		self.returns_names = [ self.process.stdout.readline().strip()
			for i in xrange(self.nreturns) ]

	def __call__(self, x):
		assert len(x) == self.nparams
		for xi in x:
			print >> self.process.stdin, xi
		self.process.stdin.flush()
		return [ float(self.process.stdout.readline()) 
			for i in xrange(self.nreturns) ]

	@staticmethod
	def from_points(X, Y, executable, quiet=False):
		"""
		Create and return a new Emulator object.  Access the filename
		of the trained emulator via:
				emulator = Emulator.from_points(X,Y,executable)
				print emulator.model_snapshot_file
		"""
		training_file = Emulator.create_training_pts_file(X,Y)
		f = tempfile.NamedTemporaryFile(
			suffix='.dat', prefix='madai_model_snapshot_', delete=False)
		model_snapshot_file = f.name
		f.close()

		output = open(os.devnull,"w") if quiet else None
		proc = subprocess.Popen(
			[ executable, 'estimate_thetas', training_file,
			  model_snapshot_file ],
			stdout=output, stderr=output )
		proc.wait()
		return Emulator(executable, model_snapshot_file)

	@staticmethod
	def create_training_pts_file(X,Y):
		"""
		Creates an input file suitable for interactive_emulator's
		estimate_thetas function.
		"""
		o = tempfile.NamedTemporaryFile(
			suffix='.dat', prefix='madai_training_pts_', delete=False)
		filename = o.name
		print >>o, len(Y[0])
		print >>o, len(X[0])
		print >>o, len(X)
		for x in X:
			print >>o, ' '.join(map(str,x))
		for y in Y:
			print >>o, ' '.join(map(str,y))
		o.close()
		return filename


def maximum_latin_hypercube(nparams,nsamples):
	"""
	Make use of a R script to call the maximinLHS() function.
	"""
	exe = './latin-hyper-sample'
	proc = subprocess.Popen([ exe, str(nparams), str(nsamples) ],
		stdout=subprocess.PIPE)
	return numpy.genfromtxt(proc.stdout, delimiter=' ')

def random_latin_hypercube(nparams,nsamples, seed=None):
	"""
	If maximum_latin_hypercube is unavailible, use this.
	"""
	if seed is not None:
		state = numpy.random.get_state()
		numpy.random.seed(seed)
	X = numpy.empty((nsamples,nparams))
	linspace = numpy.linspace(0, 1, nsamples)
	for i in xrange(nparams):
		X[:,i] = numpy.random.permutation(linspace)
	if seed is not None:
		numpy.random.set_state(state)
	return X

def test_the_emulator(nparams, nsamples, model):
	"""
	model must take in a list (or iterable castable to a list) and
	return a list.
	"""
	#X = maximum_latin_hypercube(nparams,nsamples)
	testingseed = 1234
	X = random_latin_hypercube(nparams, nsamples, testingseed)
	Y = [ model(list(x)) for x in X ]
	emulator = Emulator.from_points(X,Y,executable, quiet=True)
	print emulator.model_snapshot_file,
	print max(map(abs, (Y[i][0] - emulator(X[i])[0] for i in xrange(nsamples))))


if __name__ == '__main__':
	nparams = 6
	nsamples = nparams * 10

	# test_model1 = lambda x: [ sum(xi**2 for xi in x), ]
	# test_the_emulator(nparams, nsamples, test_model1)

	# test_model2 = lambda x: [
	# 	(1e-3 * numpy.random.random_sample()) + 
	# 	sum(xi**2 for xi in x), ]
	# test_the_emulator(nparams, nsamples, test_model2)

	test_model3 = lambda x: [
		(1e-2 * numpy.random.random_sample()) + 
		sum(xi**2 for xi in x), ]
	test_the_emulator(nparams, nsamples, test_model3)

	test_the_emulator(nparams, nsamples, dice_model)
