#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""

import numpy as np
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
from scipy.optimize import curve_fit


def sinus(omega, ampl, phase, offset):
	return ampl * np.sin(np.radians(omega - phase)) + offset


def fitsinus(omega, ydata, parm, mib, mab):
	yd_mean = np.mean(ydata) - 20
	yd_len = len(ydata)

	try:
		popt, pcov = curve_fit(
			sinus, omega, ydata,
			p0=parm,
			bounds=(mib, mab),
			method='trf', max_nfev=30000)  # kwargs={"max_nfev": 300})
		return popt, pcov
	except ValueError:
		print "Value error"

path = '/Users/andcj/bitsync/Documents/PostDoc/scripts/output/061016-1115'

peakarray = np.load(path + '/peakarray_smooth.npy')

line = 300
npeaks = np.shape(peakarray)[2]
print npeaks
peakarray[:50, :, :] = 0
peakarray[305:401, :, :] = 0
peakarray[674:, :, :] = 0


print np.shape(peakarray)

diffrx = np.arange(len(peakarray[:, 0, 0]))

comb_array = np.zeros((721 * npeaks * 2, 2))
peakarray_4peak = np.delete(peakarray, [0, 1, 2, 3, 8, 9, 10, 11, 12, 13], 2)

for l in range(len(peakarray_4peak[0, 0, :])):
	comb_array[l * 721:(l + 1) * 721, 0] = diffrx
	comb_array[l * 721:(l + 1) * 721, 1] = peakarray_4peak[:, line, l]

print np.shape(comb_array)

comb_array = np.delete(comb_array, np.where(comb_array[:, 1] == 0), 0)

# peakarray = np.delete(peakarray, np.where(peakarray[:, line, :] == 0), 0)

plt.plot(comb_array[:, 0] / 2, comb_array[:, 1], '.')

mib = [190., 330, 450.]
mab = [200., 360., 490.]
popt1, pcov = fitsinus(comb_array[:, 0], comb_array[:, 1], [200, 330, 470], mib, mab)
sine_ydata1 = sinus(diffrx / 2, popt1[0], popt1[1], popt1[2])
plt.plot(diffrx / 2, sine_ydata1)

mib = [320., 200, 450.]
mab = [np.inf, 260., 520.]
popt2, pcov = fitsinus(comb_array[:, 0], comb_array[:, 1], [360, 200, 470], mib, mab)
sine_ydata2 = sinus(diffrx / 2, popt2[0], popt2[1], popt2[2])
plt.plot(diffrx / 2, sine_ydata2)

mib = [190., 0, 450.]
mab = [200., 100., 490.]
popt3, pcov = fitsinus(comb_array[:, 0], comb_array[:, 1], [200, 50, 470], mib, mab)
sine_ydata3 = sinus(diffrx / 2, popt3[0], popt3[1], popt3[2])
plt.plot(diffrx / 2, sine_ydata3)

print popt1
print popt2
print popt3

fig0, ax0 = plt.subplots()
for j in range(npeaks):
	ax0.scatter(range(1000), peakarray[50, :, j])

fig1, ax1 = plt.subplots()

for i in range(16):
	if i < npeaks / 2:
		ax1.plot(peakarray[:, line, i], '.', color='black')
	else:
		ax1.plot(peakarray[:, line, i], '.', color='red')

plt.show()
