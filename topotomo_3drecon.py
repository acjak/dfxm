#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import time
from lib.getedfdata import *
# from lib.gauss import *
from lib.dfxm import *
import numpy as np
import itertools
from mpi4py import MPI
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
from scipy.signal import argrelextrema, savgol_filter
# import seaborn as sns
# import scipy.ndimage

# import warnings


def line(x, a, b):
	return a * x + b


def fitLine(x, y):
	from scipy.optimize import curve_fit

	try:
		popt, pcov = curve_fit(line, x, y, p0=[0, 30], maxfev=10000)
		return popt, pcov
	except RuntimeError:
		pass


def plot(img0, img1, s):
	ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4), dtype=np.uint8) * 0

	ta[:, :, 3] = 255
	ta[:, :, 0] = 255 * img0 / np.max(img0)
	ta[:, :, 2] = 255 * img0 / np.max(img0)
	ta[:, :, 1] = 255 * img1 / np.max(img1)

	if np.mean(img0) < 0.01 and np.mean(img1) < 0.01:
		ta[:, :, 3] = 255
		ta[:, :, 0] = 0
		ta[:, :, 2] = 0
		ta[:, :, 1] = 0

	fig, ax = plt.subplots()

	ax.imshow(ta[:, 1::2, :], interpolation="none")
	ax.autoscale(enable=False)
	# ax.text(50, 50, str(local_c[i]), color='white')
	ax.plot(
		[s - 40 - 2 * 56, s - 40],
		[s - 65, s - 65],
		linewidth=3,
		color='white')

	# ax.plot(
	# 	[30, 30],
	# 	[0, s],
	# 	color='blue')

	ax.text(s - 170, s - 25, u'10 μm', color='white', fontsize=16)

	fig1, ax1 = plt.subplots()
	ax1.plot(img0[:, 30], color=(1, 0, 1, 1))
	ax1.plot(img1[:, 30], color=(0, 1, 0, 1))

	fig1.savefig(
		data.directory +
		'/topo_' +
		str('%04d' % (i + rank * local_n)) +
		'_line.png')

	fig.savefig(
		data.directory +
		'/topo_' +
		str('%04d' % (i + rank * local_n)) +
		'_im.png')
	fig.clf()


def findpeaks(x):
	xsmooth = savgol_filter(x, 5, 2)
	peaks = np.zeros((npeaks))
	maxima = argrelextrema(xsmooth, np.greater, order=5, mode='wrap')
	maxvals = np.sort(xsmooth[maxima])[-npeaks:]
	for i, val in enumerate(maxvals):
		peaks[i] = np.where(xsmooth == val)[0]
	return peaks

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

# rank = 0
# mpisize = 1

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/ff_topo_2'
bg_path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff1_'
filename2 = 'ff2_'
sampletitle = filename
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'topotomo'

test_switch = False

# poi = [500, 500]
poi = [500, 500]
size = [1000, 1000]
s = 1000

npeaks = 8

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]


data = GetEdfData(
	path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

data2 = GetEdfData(
	path, filename2, bg_path, bg_filename, roi, datatype, test_switch)
data2.setTest(True)
data2.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0

tools = DFXM(
	path,
	data.data_files,
	directory,
	roi,
	datatype,
	data.dirhash,
	data.meta,
	test_switch)

a, b, f = data.getMetaValues()
c, d, g = data2.getMetaValues()
ab_vals = list(itertools.product(a, b))
cd_vals = list(itertools.product(c, d))

local_n = len(a) / mpisize
istart = rank * local_n + 50
istop = (rank + 1) * local_n
local_a = a[istart:istop]
local_c = c[istart:istop]
local_f = f[istart:istop]
local_g = g[istart:istop]

fig = plt.figure(frameon=False)
# fig.set_size_inches(5, 5)
# ax = plt.Axes(fig, [0., 0., 1., 1.])
# ax.set_axis_off()
# fig.add_axes(ax)

peak_array = np.zeros((len(local_a), size[0], npeaks * 2))

plt.figure(figsize=(20, 16))
gs1 = matplotlib.gridspec.GridSpec(8, 8)
gs1.update(wspace=0.025, hspace=0.03)

for (i, a) in enumerate(local_a):
	print i
	if i == 1:
		break
	print "Working on angle ", str(a), "."
	# fig.set_size_inches(11, 11)
	# ax = plt.Axes(fig, [0., 0., 1., 1.])
	# ax.set_axis_off()
	# fig.add_axes(ax)

	index = data.getIndex(a, -10000, -10000)
	img0 = data.getImage(index[0], False)

	index2 = data2.getIndex(a, -10000, -10000)
	img1 = data2.getImage(index2[0], False)

	imgsum1 = np.sum(img1, 1) / len(img1[0, :])
	ran = np.array(range(len(imgsum1)))
	popt, pcov = fitLine(ran, imgsum1)
	fittedline1 = ran * popt[0] + popt[1]
	fittedline1 = fittedline1 - fittedline1[len(fittedline1) / 2]
	gradient1 = np.tile(fittedline1, (len(img1[:, 0]), 1)).transpose()
	img1 -= gradient1

	imgsum0 = np.sum(img0, 1) / len(img0[0, :])
	ran = np.array(range(len(imgsum0)))
	popt, pcov = fitLine(ran, imgsum0)
	fittedline0 = ran * popt[0] + popt[1]
	fittedline0 = fittedline0 - fittedline0[len(fittedline0) / 2]
	gradient0 = np.tile(fittedline0, (len(img0[:, 0]), 1)).transpose()
	img0 -= gradient0

	img0[img0 < 0] = 0
	img1[img1 < 0] = 0

	img0 = img0 / np.max(img0)
	img1 = img1 / np.max(img1)

	plot(img0, img1, s)

	for j, x in enumerate(img0.transpose()):
		try:
			peak_array[i, j, :npeaks] = findpeaks(x)
		except ValueError:
			pass
			# print "Value error 0."

	for k, x in enumerate(img1.transpose()):
		try:
			peak_array[i, k, npeaks:] = findpeaks(x)
		except ValueError:
			pass
			# print "Value error 1."

	# for line in range(s / 50):
	# 	axarr = plt.subplot(gs1[line])
	# 	xsmooth = savgol_filter(img0[:, line * 50], 5, 2)
	# 	axarr.plot(xsmooth, color=(1, 0, 1, 1))
	#
	# 	axarr.plot(img1[:, line * 50], color=(0, 1, 0, 1))
	# 	print line * 50
	# 	for peak in range(npeaks * 2):
	# 		if peak < npeaks:
	# 			xval = peak_array[i, line * 50, peak]
	# 			axarr.plot([xval, xval], [0, 1], color=(1, 0, 1, 1))
	# 		else:
	# 			xval = peak_array[i, line * 50, peak]
	# 			axarr.plot([xval, xval], [0, 1], color=(0, 1, 0, 1))

	# plt.savefig(data.directory + '/peakfind_array.png')
	# plot(img0, img1, s)

# np.save(data.directory + '/peakarray_smooth.npy', peak_array)
# if rank == 0:
# 	end = time.time()
# 	print "Time:", end - start
