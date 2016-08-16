# !/bin/python
"""blah."""
import EdfFile
from lib.getedfdata import *
from lib.dfxm import *

import numpy as np
import scipy.interpolate as inter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
#  import seaborn as sns

# import scipy.ndimage
import itertools
import sys
import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

print rank, mpisize

if rank == 0:
	start = time.time()

# for i in range(len(sys.argv)):
# 	print i, sys.argv[i]

path = '/u/data/andcj/Resolution_march_2016/rollscan_big'
bg_path = '/u/data/andcj/Resolution_march_2016/rollscan_top_far'

filename = 'scan1_'
bg_filename = 'scan1_0960.edf'
# filename2 = 'ff2_'
sampletitle = 'rollscan_big'
datatype = 'res_paper'

test_switch = True

poi = [1024, 1024+500]
# size = [2, 20]
size = [40, 40]

pixelsize = 180.  # 180.

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)

try:
	directory = data.directory
except AttributeError:
	directory = 0
tools = DFXM(path, data.data_files, directory, roi, datatype, data.dirhash, data.meta, test_switch)

data.setTest(True)
data.adjustOffset(False)

meta = data.getMetaArray()
a, b, c = data.getMetaValues()


def plotImageArray(sampletitle):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025,  hspace=0.03)

	for i in range(len(data.imgarray[:, 0, 0])):
		img = data.imgarray[i, :, :]
		axarr = plt.subplot(gs1[i])

		axarr.imshow(img, cmap="Greens")
		axarr.set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	plt.savefig(data.directory + '/%s_array.pdf' % (sampletitle))


def makeIndexList_ROLL(a, b, c):
	index_list = []

	for i in a:
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(i, b[len(b)/2], c[0])
		index_list.append(index[0])

	return index_list

def makeIndexList_TT(a, b, c):
	index_list = []

	for i in b:
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(a[len(a)/2], i, c[0])
		index_list.append(index[0])

	return index_list

def getIntList(index_list):
	center_int = np.zeros((len(index_list)))

	srcur = meta[index_list, 3]
	srcur_max = max(srcur)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')

		for i in range(len(data.imgarray)):
			center_int[i] = srcur_max*np.sum(data.imgarray[i])/srcur[i]
		return center_int

def fitMyShitUp(center_int, xr):
	popt, pcov = tools.fitGaussian(xr, center_int)
	return popt[0], popt[1], popt[2]

def makeGaussian(xr, ampl, midp, sigma):
	xr_fine = np.arange(min(xr), max(xr), abs(xr[1]-xr[0])/10)
	yr = []
	for i in xr_fine:
		yr.append(tools.gaus(i, ampl, midp, sigma))
	return xr_fine, yr

def convertToDegrees(b):
	detx_specular = data.beta0
	# tt0 = 10.89*2
	sample_det = 5720
	b_degrees = []
	for i, detx in enumerate(b):
		b_degrees.append(math.degrees(math.atan(detx/sample_det)))
		print i, detx, math.degrees(math.atan(detx/sample_det))
	return b_degrees

index_list_tt = makeIndexList_TT(a, b, c)

index_list_roll = makeIndexList_ROLL(a, b, c)


b_d = convertToDegrees(b)
print b[len(b)/2], a[len(a)/2]

# xr_common = b_d

center_int_tt = getIntList(index_list_tt)
center_int_roll = getIntList(index_list_roll)

ampl_tt, midp_tt, sigma_tt = fitMyShitUp(center_int_tt, b_d)
ampl_roll, midp_roll, sigma_roll = fitMyShitUp(center_int_roll, a)

xr_tt, yr_tt = makeGaussian(b_d, ampl_tt, midp_tt, sigma_tt)
xr_roll, yr_roll = makeGaussian(a, ampl_roll, midp_roll, sigma_roll)

fwhm_tt = 2*math.sqrt(2*math.log(2))*sigma_tt
fwhm_roll = 2*math.sqrt(2*math.log(2))*sigma_roll

fig, ax = plt.subplots(figsize=(7,7))
# fig2, ax2 = plt.subplots(figsize=(7,7))
ax2 = ax.twiny()


legend_tt = '2T Ampl:' + str(ampl_tt) + '\n 2T Midp:' + str(midp_tt) + '\n 2T FWHM:' + str(fwhm_tt)
legend_roll = 'Roll Ampl:' + str(ampl_roll) + '\n Roll Midp:' + str(midp_roll) + '\n Roll FWHM:' + str(fwhm_roll)

ln1 = ax.plot(xr_tt, yr_tt, label=legend_tt, color='red')
ax.plot(b_d, center_int_tt, color='black', linestyle=':')
ax.set_xlabel('2theta [degrees]')

ln2 = ax2.plot(xr_roll, yr_roll, label=legend_roll)
ax2.plot(a, center_int_roll, color='black', linestyle=':')
ax2.set_xlabel('samrz [degrees]')

lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0, prop={'size':8})
ax.ticklabel_format(useOffset=False, axis='x')

# ax.legend(loc='upper right', prop={'size':8})
# ax2.legend(loc='upper right', prop={'size':8})


fig.savefig(data.directory + '/' + sampletitle + '_tt_offaxis_bot500.pdf')
# fig2.savefig(data.directory + '/' + sampletitle + '_roll_onaxis.pdf')

# plotImageArray(sampletitle)
# tools.plotPPOS(gaussarray)
# tools.plotPPOSBig(gaussarray)
# saveGaussArray(sampletitle, gaussarray)
# tools.plotStrain(strainpic, data.imgarray)
