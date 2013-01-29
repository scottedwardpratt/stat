#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Dice Model
Copyright 2012-2013, The University of North Carolina at Chapel Hill.

This software was written in 2012 by Hal Canary <cs.unc.edu/~hal/>
while working for the MADAI project <http://madai.us/>.

See <https://madai-public.cs.unc.edu/software/license/> for license
information.

This is intended as an example of:

-   A parametraized model of a physical process

-   A comparison of the model to an expiriment to produce a liklihood estimate.

-   The Metropolis-Hastings MCMC algorithm,

"""
import numpy
import random
import math, sys

def dice_model(x,N=96):
	"""
	This function takes in a point in parameter space (x) and models N
	roles of a weighted dice.
	Returns a histogram of the results as a list.  Will sum to N.
	"""
	(x1,x2,x3,x4,x5,x6) = x
	P1, P2, P3, P4 = x1, x2, x3, x4
	P5 = x5 + x6
	P6 = 1.0 - x1 - x2 - x3 - x4 - x5 - x6
	count = [0 for i in xrange(6)]
	for i in xrange(N):
		r = random.random()
		if r < P1:
			count[0] += 1
		elif r < P1 + P2:
			count[1] += 1
		elif r < P1 + P2 + P3:
			count[2] += 1
		elif r < P1 + P2 + P3 + P4:
			count[3] += 1
		elif r < P1 + P2 + P3 + P4 + P5:
			count[4] += 1
		else:
			count[5] += 1
	assert sum(count) == N
	return count[:5]

def dice_model2(x,N=96):
	"""
	reimplementation of dice_model().
	"""
	(x1,x2,x3,x4,x5,x6) = x
	P = [x1 * 0.3, x2 * 0.3, x3 * 0.3, x4 * 0.3, 0, 0]
	P[4] = (x5 + x6) * 0.3
	P[5] = 1.0 - x1 - x2 - x3 - x4 - x5 - x6
	count = [0 for i in xrange(6)]
	#count = [int(math.floor(P[i]*N)) for i in xrange(6)]
	#for i in xrange(N - sum(count)):
	for i in xrange(N):
		r = random.random()
		if r < sum(P[:1]):
			count[0] += 1
		elif r < sum(P[:2]):
			count[1] += 1
		elif r < sum(P[:3]):
			count[2] += 1
		elif r < sum(P[:4]):
			count[3] += 1
		elif r < sum(P[:5]):
			count[4] += 1
		else:
			count[5] += 1
	assert sum(count) == N
	return count[:5]

def liklihood(y, expected, covariance):
	return numpy.exp(-0.5 * numpy.dot(
		(y - expected),
		numpy.dot( 
			numpy.linalg.inv(covariance),
			(y - expected)).flatten()))

def get_likelihood(model, x, expected, covariance):
	"""
	Given a model function and a point x, calculate the liklihood of
	that point.
	"""
	return liklihood(numpy.array(model(x),dtype=float), expected, covariance)
	# invcov = numpy.linalg.inv(covariance),
	# #expected = numpy.array(expected, dtype=float)
	# #variance = numpy.array(variance, dtype=float)
	# J2 = numpy.dot(
	# 	(y - expected),
	# 	numpy.dot( invcov, (y - expected)).flatten())
	# return numpy.exp(-J2/2.0)

def MCMC(func, N, skip=0, length_scale=0.1, init_x=None):
	"""
	Metropolis-Hastings algorithm
	func should be a tuple including the number of arguments.
	returns a trace
	"""
	argnum, f = func
	assert skip >= 0
	skip += 1
	trace = numpy.empty((argnum,N))
	count = total = 0
	x = tuple(0.0 for i in xrange(argnum)) if init_x is None else init_x
	val = f(*x)
	while total < N:
		r = numpy.random.uniform(-length_scale, length_scale, argnum)
		x2 = tuple(a + b for a,b in zip(x,r))
		val2 = f(*x2)
		if val2 > 0.0:
			a = float(val2) / float(val) if (float(val) > 0) else 1.0
			if (a >= 1.0) or (a > numpy.random.uniform(0.0,1.0)):
				val,x = val2,x2
				count += 1
				if count%skip == 0:
					trace[:,total] = x
					total +=1
	return trace

def print_matrix(matrix, transpose=False, file_object=sys.stdout):
	"""
	I don't like Numpy's default print methods.
	"""
	view = matrix.transpose() if transpose else matrix
	for i in xrange(view.shape[0]):
		print >>sys.stdout, ','.join(map(str,view[i,:]))


if __name__ == '__main__':
	# x = [0.16,0.20,0.16,0.16,0.16,0.0]
	# Yf = [dice_model(x, 96) for i in xrange(8)]
	# print '\n'.join(map(repr, Yf))
	"""
	Yf is the matrix of field data (real-world observations) The
	experiment was repeated eight times.  To remove linear dependance, the
	final column was deleted.  Otherwise the covariace calulation will
	give wierd answers.
	[[Here's a pit to fall into.]]
	"""
	Yf = numpy.array(          # Yf=numpy.array(
		[[16, 15, 16, 19, 16], # 	[[16, 15, 16, 19, 16, 14],
		 [12, 21, 21, 15, 18], # 	 [12, 21, 21, 15, 18,  9],
		 [11, 17, 15, 21, 12], # 	 [11, 17, 15, 21, 12, 20],
		 [12, 20, 21, 13, 16], # 	 [12, 20, 21, 13, 16, 14],
		 [14, 19, 18, 14, 19], # 	 [14, 19, 18, 14, 19, 12],
		 [18, 19, 15, 14, 17], # 	 [18, 19, 15, 14, 17, 13],
		 [19, 18, 14, 16, 15], # 	 [19, 18, 14, 16, 15, 14],
		 [17, 15, 19, 18, 13]] # 	 [17, 15, 19, 18, 13, 14]]
		).transpose()		   # 	).transpose()

	"""
	EYf = E[Y_f] = expectation value of the observations.
	"""
	EYf = numpy.mean(Yf, 1)
	"""
	CovYf = covariace of the observations.  If only one set of
	observations, is just variance.
	"""
	CovYf = numpy.cov(Yf, bias=1)

	"""
	To calulate likeyhood ratio, compare the output of the dice_model at a
	point x to the expected value EYf, using the Covariance of the output
	to scale.
	"""
	# dice_likelihood = lambda *x : get_likelihood(dice_model2, x, EYf, CovYf)
	dice_likelihood = lambda *x : liklihood(
		numpy.array(dice_model2(x),dtype=float), EYf, CovYf)

	trace =  MCMC((6, dice_likelihood), 5000)

	print_matrix(trace, True)



#VarYf = numpy.var(Yf,1)
#CovYf = numpy.var(Yf)
#y = numpy.array(dice_model([0.12,0.1,0.1,0.1,0.1,0.1]),dtype=float)
#y = numpy.array(dice_model2(x),dtype=float)
# x = [0.16,0.20,0.16,0.16,0.16,0.0]
# y = numpy.array(dice_model(x),dtype=float)
# J2 = numpy.dot((y-EYf),numpy.dot(numpy.linalg.inv(numpy.diag(VarYf)),(y-EYf)))
# print math.sqrt(J2), '\t', numpy.exp(-J2/2.0)
# x = [0.12,0.1,0.1,0.1,0.1,0.1]
# y = numpy.array(dice_model(x),dtype=float)
# J2 = numpy.dot((y-EYf),numpy.dot(numpy.linalg.inv(numpy.diag(VarYf)),(y-EYf)))
# print math.sqrt(J2), '\t', numpy.exp(-J2/2.0)
# def get_likelihood(model, x, expected, variance):
# 	y = numpy.array(model(x),dtype=float)
# 	expected = numpy.array(expected, dtype=float)
# 	variance = numpy.array(variance, dtype=float)
# 	J2 = numpy.dot(
# 		(y-expected),
# 		numpy.dot(
# 			numpy.linalg.inv(numpy.diag(variance)),
# 			(y-expected)))
# 	return numpy.exp(-J2/2.0)
#print EYf/96
#print numpy.mean(trace,1)
