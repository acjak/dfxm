# !/bin/python
"""blah."""
import numpy as np
import scipy.interpolate as inter
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import seaborn as sns

# import scipy.ndimage
import itertools
import sys
import warnings


def gaus(x, a, x0, sigma):
	return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def fitGaussian(x, y):
	from scipy.optimize import curve_fit

	try:
		popt, pcov = curve_fit(gaus, x, y, p0=[max(y), x[np.argmax(y)], .003], maxfev=10000)
		return popt, pcov
	except RuntimeError:
		pass


def fitMyShitUp(center_int, xr):
	popt, pcov = fitGaussian(xr, center_int)
	return popt[0], popt[1], popt[2]


def makeGaussian(xr, ampl, midp, sigma):
	xr_fine = np.arange(min(xr), max(xr), abs(xr[1] - xr[0]) / 10)
	yr = []
	for i in xr_fine:
		yr.append(gaus(i, ampl, midp, sigma))
	return xr_fine, yr


BG = 102.
scan1 = np.loadtxt('input/scan524.txt')

diffrz = scan1[:, 0]
ccdavg = scan1[:, 7]

pop1, pop2, pop3 = fitMyShitUp(ccdavg - BG, diffrz)

xr_fine, yr = makeGaussian(diffrz, pop1, pop2, pop3)


# mpy = len(detnorm[0,:])/2
# mpx = len(detnorm[:,0])/2
fwhm = 2 * math.sqrt(2 * math.log(2)) * pop3
ampl = (pop3 * math.sqrt(2 * math.pi)) / pop1


gaussoutput = [pop1, pop2, pop3]

np.savetxt('output/gaussoutput.txt', gaussoutput)
np.savetxt('output/gaussfit.txt', np.vstack((xr_fine, yr)).T)
np.savetxt('output/rocking.txt', np.vstack((diffrz, ccdavg - BG)).T)

# sns.set_style('white')
# fig, ax = plt.subplots(figsize=(7, 7))
# # fig2, ax2 = plt.subplots(figsize=(7,7))
# # ax2 = ax.twiny()
#
#
# legend_tt = 'Ampl:' + str(pop1) + '\n Midp:' + str(pop2) + '\n FWHM:' + str(fwhm)

# ax.plot(diffrz, ccdavg - BG)
# ax.plot(xr_fine, yr, label=legend_tt)
# ax.legend(prop={'size': 10})
# ax.ticklabel_format(useOffset=False, axis='x')
# ax.set_title('Rocking curve, diamond, resolution experiment')
# ax.set_xlabel('diffrz [degrees]')
# ax.set_ylabel('counts (BG subtracted)')
#
# fig.savefig('rockingcurve.png')

# legend_roll = 'Roll Ampl:' + str(ampl_roll) + '\n Roll Midp:' + str(midp_roll) + '\n Roll FWHM:' + str(fwhm_roll)
#
# ln1 = ax.plot(xr_tt, yr_tt, label=legend_tt, color='red')
# ln1_data = ax.plot(b_d, tt_stdlist, color='black', linestyle=':', label="2T STD")
# ax.set_xlabel('2theta [degrees]')
#
# ln2 = ax2.plot(xr_roll, yr_roll, label=legend_roll)
# ln2_data = ax2.plot(a, roll_stdlist, color='black', linestyle='-.', label="Roll STD")
# ax2.set_xlabel('samrz [degrees]')
#
# lns = ln1+ln1_data+ln2+ln2_data
# labs = [l.get_label() for l in lns]
# ax.legend(lns, labs, loc=0, prop={'size':8})
# ax.ticklabel_format(useOffset=False, axis='x')

# ax.legend(loc='upper right', prop={'size':8})
# ax2.legend(loc='upper right', prop={'size':8})


# fig.savefig(data.directory + '/' + sampletitle + '_tt_onaxis_' + str(roi) + '.pdf')
