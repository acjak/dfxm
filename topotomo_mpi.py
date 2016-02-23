#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""

from lib.getedfdata import *
from lib.gauss import *
import numpy as np

import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
# import seaborn as sns

# import scipy.ndimage
import itertools

# import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

# rank = 0
# mpisize = 1

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/data/hxrm/Dislocation_november_2015/diamond/ff_topo_2'
bg_path = '/data/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff1_'
filename2 = 'ff2_'
sampletitle = filename
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'topotomo'

test_switch = True

# poi = [500, 500]
poi = [500, 650]
size = [400, 400]
s = 400

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]


data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

data2 = GetEdfData(path, filename2, bg_path, bg_filename, roi, datatype, test_switch)
data2.setTest(True)
data2.adjustOffset(False)


a, b, f = data.getMetaValues()
c, d, g = data2.getMetaValues()
ab_vals = list(itertools.product(a, b))
cd_vals = list(itertools.product(c, d))

print len(a), len(b), len(f)
print len(c), len(d), len(g)

local_n = len(a)/mpisize
istart = rank*local_n
istop = (rank+1)*local_n
local_a = a[istart:istop]
local_c = c[istart:istop]
local_f = f[istart:istop]
local_g = g[istart:istop]

fig = plt.figure(frameon=False)
# fig.set_size_inches(5, 5)
# ax = plt.Axes(fig, [0., 0., 1., 1.])
# ax.set_axis_off()
# fig.add_axes(ax)

# print b, f, d, g

for j in [1]:  # range(len(local_a)):
	i = 436
	fig.set_size_inches(6, 6)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)
	index = data.getIndex(float(local_a[i]), float(b[0]), float(local_f[i]))
	print local_a[i], f[i], index
	img0 = data.getImage(index[0], False)

	index2 = data2.getIndex(float(local_c[i]), float(d[0]), float(local_g[i]))
	print local_c[i], g[i], index2
	img1 = data2.getImage(index2[0], False)

	ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

	ta[:, :, 3] = 255
	ta[:, :, 0] = 255*img0/np.max(img0)
	ta[:, :, 2] = 255*img0/np.max(img0)
	ta[:, :, 1] = 255*img1/np.max(img1)

	if np.mean(img0) < 0.01 and np.mean(img1) < 0.01:
		ta[:, :, 3] = 255
		ta[:, :, 0] = 0
		ta[:, :, 2] = 0
		ta[:, :, 1] = 0

	ax.imshow(ta, interpolation="none", cmap="Greens")
	ax.autoscale(enable=False)
	# ax.text(50, 50, str(local_c[i]), color='white')
	ax.plot([s-20-2*56, s-20], [s-40, s-40], linewidth=5, color='white')
	ax.text(s-100, s-20, u'20 μm', color='white')

	# ax.clf()
	# plt.colorbar()
	fig.savefig(data.directory + '/topo_im_' + str('%04d' % (i+rank*local_n)) + '.png')
	fig.clf()

if rank == 0:

	# hist, datx, daty = data.makeMeanGrid()
	# hist_array = np.zeros((len(hist[:, 0]), len(hist[0, :]), 1))
	# hist_array[:, :, 0] = hist
	# data.makeHistogram(hist_array, a, b, 'ff_topo')

	end = time.time()
	print "Time:", end-start
