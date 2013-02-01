#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
statistics.py
Copyright 2012-2013, The University of North Carolina at Chapel Hill.

This software was written in 2012 by Hal Canary <cs.unc.edu/~hal/>
while working for the MADAI project <http://madai.us/>.

See <https://madai-public.cs.unc.edu/software/license/> for license
information.

This is intended as an example of:

-  A comparison of a parametraized model of a physical process to an
   expiriment to produce a liklihood estimate.
-  The Metropolis-Hastings MCMC algorithm,

"""
import numpy
import random
import math, sys
from dice_model import dice_model

def log_liklihood(y, expected, inverse_covariance):
	diff = y - expected
	return -0.5 * numpy.dot( diff,numpy.dot(inverse_covariance, diff).flatten())

def MCMC(log_function, number_of_parameters, N, length_scale=0.1, init_x=None):
	"""
	an implementation of the Metropolis-Hastings algorithm.

	log_function should return the log of the distribution to be sampled.

	number_of_parameters is the number of arguments to pass to log_function

	N is the number of samples to return.
	"""
	trace = numpy.empty((number_of_parameters,N))
	count = 0
	x = [0.5 for i in xrange(number_of_parameters)] if init_x is None else init_x
	val = log_function(*x)
	while count < N:
		x2 = [ random.gauss(xi,length_scale) for xi in x ]
		val2 = log_function(*x2)
		a = val2 - val
		if (a > 0.0) or (math.exp(a) > numpy.random.uniform(0.0,1.0)):
			val,x = val2,x2
			trace[:,count] = x
			count += 1
	return trace

def main():
	"""
	Yf is the matrix of field data (real-world observations) The
	experiment was repeated eight times.  To remove linear dependance, the
	final column was deleted.  Otherwise the covariace calulation will
	give wierd answers.
	[[Here's a pit to fall into.]]

	Yf was created with this:
		from dice_model import dice_model
		GROUND_TRUTH = [8.0/15,2.0/3,8.0/15,8.0/15,8.0/15,0.0]
		print repr([dice_model(GROUND_TRUTH) for i in xrange(8)])
	"""

	Yf = numpy.array(
		[[0.153, 0.224, 0.166, 0.158, 0.141],
		 [0.148, 0.208, 0.171, 0.16, 0.164],
		 [0.141, 0.222, 0.151, 0.172, 0.16],
		 [0.156, 0.208, 0.161, 0.176, 0.147],
		 [0.16, 0.202, 0.15, 0.148, 0.171],
		 [0.141, 0.221, 0.163, 0.156, 0.147],
		 [0.152, 0.198, 0.151, 0.156, 0.162],
		 [0.164, 0.199, 0.164, 0.158, 0.164]]
		).transpose() # N = 1000

	"""
	EYf = E[Y_f] = expectation value of the observations.
	"""
	EYf = numpy.mean(Yf, 1)
	"""
	CovYf = covariace of the observations.  If only one set of
	observations, is just variance.
	"""
	CovYf = numpy.cov(Yf, bias=1)
	inverse_CovYf = numpy.linalg.inv(CovYf)

	"""
	To calulate likeyhood ratio, compare the output of the dice_model at a
	point x to the expected value EYf, using the Covariance of the output
	to scale.
	"""
	dice_log_likelihood = lambda *x : log_liklihood(
		numpy.array(dice_model(x),dtype=float), EYf, inverse_CovYf)

	"""
	I took some time to discover this:
	"""
	best_length_scale = 1e-2

	trace = MCMC(dice_log_likelihood, 6, 5000, length_scale=best_length_scale)

	numpy.savetxt(sys.stdout,trace.transpose(),'%r',delimiter=",")

if __name__ == '__main__':
	main()
