#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import numpy as np
# import matplotlib.pyplot as plt
# import cmath
from math import sin, cos, atan, sqrt
from sympy.mpmath import acot

mu = 47.
d1 = 0.289
f = 21.195
phi = 0.00869
N = 69

twodelta = 2.359E-6
R = 5.E-5
atanterm = atan(d1 / (f * phi))

pythterm = d1**2 + f**2 * phi**2
Mx = 16.27

sa = \
	np.sqrt(R / (mu * pythterm *
	(N + 1 - (1 / phi) * sin((N + 1) * phi) *
	cos((N - 1) * phi + 2 * atanterm))
	))

saappr = \
	(twodelta / 2) * Mx / (Mx + 1) * sqrt((2 * N) / (mu * R))

print '\nsigma_a =  {:08.5}\nsigma_a (approximation) = {:08.5}'.format(sa, saappr)
print 'Difference is: {:04.2} %\n'.format(100 * (saappr - sa) / saappr)

gammaappr = 1 / sqrt(d1**2 + f**2 * phi**2)
#
# gamtop = (f + 2 * d1 * N) * sqrt(pythterm) * cos(2 * N * phi + atanterm)
# gambot = 2 * (d1 * f + N * (d1**2 + f**2 * phi**2)) + (1 / phi) * pythterm * sin(2 * N * phi)
#
# gammaappr = gamtop / gambot

AN = (mu / (2 * R)) * pythterm * (N - 1 + (1 / phi) * sin((N + 1) * phi) * cos((N - 1) * phi - 2 * acot(d1 / (f * phi))))
BN = (2 * mu / R) * ((d1 / 2) * (N - 1) + sqrt(pythterm) * sin((N + 1) * phi) / (2 * phi) * cos((N - 1) * phi - acot(d1 / (f * phi))))

gamma = BN / (2 * AN)

print 'gamma = {:08.5}\ngamma (approximation)  = {:08.5}'.format(gamma, gammaappr)
print 'Difference is: {:04.2} %\n'.format(100 * (gammaappr - gamma) / gammaappr)
