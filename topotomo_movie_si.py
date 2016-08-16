#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import time
from lib.getedfdata import *
# from lib.gauss import *
from lib.dfxm import *
import numpy as np
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
# import seaborn as sns
# import scipy.ndimage
import itertools
# import warnings
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

# rank = 0
# mpisize = 1

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/u/data/andcj/hxrm/Dislocation_may_2015/dislocations/mapping'
bg_path = '/u/data/andcj/hxrm/Dislocation_may_2015/dislocations/mapping'

filename = 'map1_'
# filename2 = 'ff2_'
sampletitle = filename
bg_filename = 'bg1_'

datatype = 'topotomo'

test_switch = True

# poi = [500, 500]
poi = [500, 500]
size = [1000, 1000]
s = 1000

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]


data = GetEdfData(
	path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

# data2 = GetEdfData(
# 	path, filename2, bg_path, bg_filename, roi, datatype, test_switch)
# data2.setTest(True)
# data2.adjustOffset(False)

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
# c, d, g = data2.getMetaValues()
ab_vals = list(itertools.product(a, b))
# cd_vals = list(itertools.product(c, d))

print a
print b
print f

local_n = len(a) / mpisize
istart = rank * local_n
istop = (rank + 1) * local_n
local_a = a[istart:istop]
local_c = c[istart:istop]
local_f = f[istart:istop]
local_g = g[istart:istop]

# fig = plt.figure(frameon=False)
# # fig.set_size_inches(5, 5)
# # ax = plt.Axes(fig, [0., 0., 1., 1.])
# # ax.set_axis_off()
# # fig.add_axes(ax)
#
# for (i, a) in enumerate(local_a):
# 	print "Working on angle ", str(a), "."
# 	fig.set_size_inches(11, 11)
# 	ax = plt.Axes(fig, [0., 0., 1., 1.])
# 	ax.set_axis_off()
# 	fig.add_axes(ax)
#
# 	index = data.getIndex(a, -10000, -10000)
# 	img0 = data.getImage(index[0], False)
#
# 	index2 = data2.getIndex(a, -10000, -10000)
# 	img1 = data2.getImage(index2[0], False)
#
# 	ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4), dtype=np.uint8) * 0
#
# 	ta[:, :, 3] = 255
# 	ta[:, :, 0] = 255 * img0 / np.max(img0)
# 	ta[:, :, 2] = 255 * img0 / np.max(img0)
# 	ta[:, :, 1] = 255 * img1 / np.max(img1)
#
# 	if np.mean(img0) < 0.01 and np.mean(img1) < 0.01:
# 		ta[:, :, 3] = 255
# 		ta[:, :, 0] = 0
# 		ta[:, :, 2] = 0
# 		ta[:, :, 1] = 0
#
# 	ax.imshow(ta, interpolation="none", cmap="Greens")
# 	ax.autoscale(enable=False)
# 	# ax.text(50, 50, str(local_c[i]), color='white')
# 	ax.plot(
# 		[s - 40 - 4 * 56, s - 40],
# 		[s - 60, s - 60],
# 		linewidth=10,
# 		color='white')
#
# 	ax.text(s - 190, s - 25, u'20 μm', color='white', fontsize=20)
#
# 	fig.savefig(
# 		data.directory +
# 		'/topo_im_' +
# 		str('%04d' % (i + rank * local_n)) +
# 		'.png')
# 	fig.clf()

# if rank == 0:
# 	end = time.time()
# 	print "Time:", end - start
