#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math

# Cylinder radius:
R = 30.

x1array = np.linspace(0., 100, 1000)
x0 = 50.
x = np.linspace(-R, R, 1000 * (2 * R) / 500)

cx1 = np.zeros((1000), dtype=np.complex)
# cx2 = np.zeros((1000), dtype=np.complex)

r0 = 40000000.
r1 = 40000.

lamb = 41.33E-6
# delta = -0.0001

rho = 10E12
scatt_ampl = 2.82E-9

delta = ((lamb ** 2) * rho * scatt_ampl) / (2 * np.pi)
print delta

r = r0 + r1

j = np.complex(0, 1)

prefix = 1 / cmath.sqrt(j * lamb * r1)
prefix2 = cmath.sqrt(r / (j * lamb * r0 * r1))

phi = ((4 * np.pi) / lamb) * delta * R * np.sqrt(1 - (x / R) ** 2)

for i, x1 in enumerate(x1array):
	# print (np.exp(j * phi) - 1)
	# print ((2 * j * np.pi) / lamb)
	# print ((x - x1) ** 2) * 2 * r1
	# print np.sum(np.exp(((2 * j * np.pi) / lamb) * ((x - x1) ** 2) * 2 * r1))

	# cx1[i] = prefix * np.sum(np.exp(((2 * j * np.pi) / lamb) * ((x - x1) ** 2) * 2 * r1) * (np.exp(j * phi) - 1))
	cx1[i] = prefix2 * np.exp((2 * np.pi * j) / lamb * (-(x1 - x0)**2) / (2 * r)) * \
		np.sum(np.exp(((2 * j * np.pi) / lamb) * ((x - x0) ** 2) / (2 * r0) + ((x1 - x) ** 2) / (2 * r1)) *
		(np.exp(j * phi) - 1))

realcx1 = np.real(cx1)

cx1[cx1 < 0] = cx1[cx1 < 0] * -1

I1 = 1 + 2 * realcx1 + cx1 ** 2

# cx2 = 2 * np.pi * 1 * 3.05E14 * lamb * 10
#
# realcx2 = np.real(cx2)
# cx2[cx2 < 0] = cx2[cx2 < 0] * -1
# I2 = 1 + 2 * realcx2 + cx2 ** 2

plt.plot(x1array, I1)
plt.show()
