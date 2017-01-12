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


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

# rank = 0
# mpisize = 1

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/nf_topo_2'
bg_path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/nf_topo_2'
# bg_path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'nf2_'
filename2 = 'nf3_'
sampletitle = filename
bg_filename = 'nf2_0001'

datatype = 'topotomo'

test_switch = True

# poi = [500, 500]
poi = [1000, 1000]
size = [2000, 2000]
s = 1500

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
istart = rank * local_n
istop = (rank + 1) * local_n
local_a = a[istart:istop]
local_c = c[istart:istop]
local_f = f[istart:istop]
local_g = g[istart:istop]

fig = plt.figure(frameon=False)
# fig.set_size_inches(5, 10)
# ax = plt.Axes(fig, [0., 0., 1., 1.])
#
# fig.add_axes(ax)

for (i, a) in enumerate(local_a):
	print "Working on angle ", str(a), "."
	fig.set_size_inches(10, 10)
	fig.set_size_inches(size[0], size[1])
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	# ax.set_axis_off()
	fig.add_axes(ax)

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

	ax.imshow(ta, interpolation="none")
	ax.autoscale(enable=False)
	# ax.text(50, 50, a, color='white')
	# ax.plot(
	# 	[s - 40 - 2 * 80, s - 40],
	# 	[s - 55, s - 55],
	# 	linewidth=3,
	# 	color='white')

	# ax.plot(
	# 	[30, 30],
	# 	[0, s],
	# 	color='blue')

	# ax.text(s - 180, s - 25, u'100 μm', color='white', fontsize=16)

	ax.set_axis_off()
	extent = ax.get_window_extent().transformed(
		plt.gcf().dpi_scale_trans.inverted())

	# fig1, ax1 = plt.subplots()
	# ax1.plot(img0[300, :], color='Greens')
	# ax1.plot(img1[300, :], color='Pink')
	#
	# fig1.savefig(
	# 	data.directory +
	# 	'/topo_im_' +
	# 	str('%04d' % (i + rank * local_n)) +
	# 	'.png')

	fig.savefig(
		data.directory +
		'/topo_im_' +
		str('%04d' % (i + rank * local_n)) +
		'.png',
		bbox_inches=extent,
		dpi=1)
	# pad_inches=0)

	# if i == 1:
	# 	break

# if rank == 0:
# 	end = time.time()
# 	print "Time:", end - start
