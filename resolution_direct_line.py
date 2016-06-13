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
# import seaborn as sns

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

# path = '/data/hxrm/Resolution_march_2016/rollscan_big'
# bg_path = '/data/hxrm/Resolution_march_2016/rollscan_top_far'

filename = 'scan1_'
bg_filename = 'scan1_0960.edf'
# filename2 = 'ff2_'
sampletitle = 'sharpness_line'
datatype = 'res_paper'

test_switch = True

poi = [1134, 1224]  # -float(sys.argv[1])]
# size = [2, 20]
size = [50, 50]

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

def plotImageArray(imgarray, sampletitle, detnorm):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025,  hspace=0.03)

	for i in range(len(imgarray[:, 0, 0])):
		img = imgarray[i, :, :]*detnorm
		axarr = plt.subplot(gs1[i])

		axarr.imshow(img, cmap="Greens")
		axarr.set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())


	plt.savefig(data.directory + '/%s_array.pdf' % (sampletitle))

def plotPlots(linearray, sampletitle, roi):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.25,  hspace=0.3)
	from scipy.signal import argrelextrema

	for i, li in enumerate(linearray):  # [0, :]):

		axarr = plt.subplot(gs1[i])

		lin = linearray[:, i]

		axarr.plot(lin)

		try:
			ma = np.mean(lin[argrelextrema(lin, np.greater)[0]])
			mi = np.mean(lin[argrelextrema(lin, np.less)[0]])
			axarr.plot([0,len(lin)],[ma,ma])
			axarr.plot([0,len(lin)],[mi,mi])
			axarr.set_title('{:.2%}'.format((ma-mi)/ma))
		except IndexError:
			print "Index error."

		axarr.set_ylim(0.6,1)
		# set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	print data.directory, sampletitle, str(roi)
	# plt.savefig('{}/{}_{}_lines.pdf'.format(data.directory, sampletitle, str(roi)))
	plt.savefig(data.directory + '/%s_%s_lines.pdf' % (sampletitle, str(roi)))

def plotHists(linearray, sampletitle, roi):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(8, 8)
	gs1.update(wspace=0.025,  hspace=0.03)

	for i in range(len(linearray[0, :])):

		hist,values = np.histogram(linearray[:, i],bins=10)
		# av_line.append(np.std(hist))

		axarr = plt.subplot(gs1[i])

		axarr.plot(hist)
		axarr.set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	plt.savefig(data.directory + '/%s_%s_hists.pdf' % (sampletitle, str(roi)))

def getSTD(imgarray, sampletitle, detnorm):
	img_std = np.zeros((len(imgarray[:, 0, 0])))
	for i in range(len(imgarray[:, 0, 0])):
		img_std[i] = np.std(imgarray[i, :, :]*detnorm)

	return img_std

def getLine(imgarray, detnorm):
	line_array = np.zeros((len(imgarray[0, 0, :]),len(imgarray[:, 0, 0])))
	for i, li in enumerate(line_array):
		line = np.sum(imgarray[i, :, :]*detnorm, axis=0)
		# line[line<100] = 1
		line_array[:, i] = line/max(line)
	return line_array

def convertToDegrees(b):
	detx_specular = data.beta0
	# tt0 = 10.89*2
	sample_det = 5720
	b_degrees = []
	for i in b:
		b_degrees.append(math.degrees(math.tan(i/sample_det)))
	return b_degrees

def fitMyShitUp(center_int, xr):
	popt, pcov = tools.fitGaussian(xr, center_int)
	return popt[0], popt[1], popt[2]

def makeGaussian(xr, ampl, midp, sigma):
	xr_fine = np.arange(min(xr), max(xr), abs(xr[1]-xr[0])/10)
	yr = []
	for i in xr_fine:
		yr.append(tools.gaus(i, ampl, midp, sigma))
	return xr_fine, yr

def getIntList(index_list):
	center_int = np.zeros((len(index_list)))

	srcur = meta[index_list, 3]
	srcur_max = max(srcur)
	print srcur
	print srcur_max

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')




# index = data.getIndex(a[25], b[25], c[0])
# img = data.getImage(index[0], False)

detnorm = np.load('tmp/detector_norm.npy')[roi[2]:roi[3],roi[0]:roi[1]]


# mpy = len(detnorm[0,:])/2
# mpx = len(detnorm[:,0])/2

index_list_tt = makeIndexList_TT(a, b, c)

print a
print b
print c

# index_list_roll = makeIndexList_ROLL(a, b, c)

b_d = convertToDegrees(b)

# xr_common = b_d

getIntList(index_list_tt)
tt_array = data.imgarray

# getIntList(index_list_roll)
# roll_array = data.imgarray

plotImageArray(tt_array, sampletitle + '_tt', detnorm)
# plotImageArray(roll_array, sampletitle + '_roll', detnorm)

tt_line = getLine(tt_array, detnorm)
# roll_line = getLine(roll_array, detnorm)


# tt_stdlist = getSTD(tt_array, sampletitle + '_tt', detnorm)
# roll_stdlist = getSTD(roll_array, sampletitle + '_roll', detnorm)
#
# ampl_tt, midp_tt, sigma_tt = fitMyShitUp(tt_stdlist, b_d)
# ampl_roll, midp_roll, sigma_roll = fitMyShitUp(roll_stdlist, a)
#
# xr_tt, yr_tt = makeGaussian(b_d, ampl_tt, midp_tt, sigma_tt)
# xr_roll, yr_roll = makeGaussian(a, ampl_roll, midp_roll, sigma_roll)
#
# fwhm_tt = 2*math.sqrt(2*math.log(2))*sigma_tt
# fwhm_roll = 2*math.sqrt(2*math.log(2))*sigma_roll

fig, ax = plt.subplots(figsize=(7,7))
# fig2, ax2 = plt.subplots(figsize=(7,7))

for i in range(len(tt_line[0, :])):
	ax.plot(tt_line[:,i])

figcenter, axcenter = plt.subplots(figsize=(7,7))

axcenter.plot(tt_line[:,25])
np.savetxt(data.directory + '/' + sampletitle + '_' + str(roi) + '.txt', tt_line[:,25])

fig2, ax2 = plt.subplots(figsize=(7,7))
# fig2, ax2 = plt.subplots(figsize=(7,7))

av_line = []

# for i in range(len(tt_line[0, :])):
# 	hist,values = np.histogram(tt_line[:,i])
# 	av_line.append(np.std(hist))
#
# ax2.plot(b_d, av_line)

plotPlots(tt_line, sampletitle, roi)
# plotHists(tt_line, sampletitle, roi)

# legend_tt = '2T Ampl:' + str(ampl_tt) + '\n 2T Midp:' + str(midp_tt) + '\n 2T FWHM:' + str(fwhm_tt)
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


fig.savefig(data.directory + '/' + sampletitle + '_' + str(roi) + '.pdf')
fig2.savefig(data.directory + '/' + sampletitle + '_average' + str(roi) + '.pdf')
figcenter.savefig(data.directory + '/' + sampletitle + '_center' + str(roi) + '.pdf')
