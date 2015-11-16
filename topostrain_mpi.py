#!/bin/python
"""blah."""

from lib.getedfdata import *
from lib.gauss import *
import numpy as np

# import matplotlib
import matplotlib.pylab as plt
# import seaborn as sns

import scipy.ndimage
import itertools

import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/data/id06/inhouse/2015/run5/diamond/ff_strain'

filename = 'ff2_'
#filename2 = 'ff2_'
sampletitle = filename
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'strain_tt'

poi = [500, 500]
size = [1000, 1000]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]


data = GetEdfData(path, filename, bg_filename, roi, datatype)
data.setTest(True)
data.adjustOffset(False)

#data2 = GetEdfData(path, filename2, bg_filename, roi, datatype)
#data2.setTest(True)
#data2.adjustOffset(False)


a, b = data.getMetaValues()
#c, d = data2.getMetaValues()
ab_vals = list(itertools.product(a, b))
#cd_vals = list(itertools.product(c, d))

# print a

local_n = len(a)/mpisize
istart = rank*local_n
istop = (rank+1)*local_n
local_a = a[istart:istop]
#local_c = c[istart:istop]


fig = plt.figure(frameon=False)
fig.set_size_inches(5,5)
ax = plt.Axes(fig,[0.,0.,1.,1.])
ax.set_axis_off()
fig.add_axes(ax)

#for i in range(len(local_a)):
	#print local_a[i], local_c[i]
	#index = data.getIndex(float(local_a[i]), float(b[0]))
	#img0 = data.getImage(index[0], False)
	
	#print np.mean(img0)
	
	#index = data2.getIndex(float(local_c[i]), float(d[0]))
	#img1 = data2.getImage(index[0], False)
	
	#ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0
	
	#ta[:, :, 3] = 255
	#ta[:, :, 0] = 255*img0/np.max(img0)
	#ta[:, :, 2] = 255*img0/np.max(img0)
	#ta[:, :, 1] = 255*img1/np.max(img1)
	
	#if np.mean(img0) < 0.01 and np.mean(img1) < 0.01:
		#ta[:, :, 3] = 255
		#ta[:, :, 0] = 0
		#ta[:, :, 2] = 0
		#ta[:, :, 1] = 0
		
	
	#ax.imshow(ta,interpolation="none", cmap = "Greens")
	##plt.colorbar()
	#fig.savefig(data.directory + '/topo_im_' + str('%04d' % (i+rank*local_n)) + '.png')
	##fig.clf()

if rank == 0:

	hist, datx, daty = data.makeMeanGrid()
	hist_array = np.zeros((len(hist[:, 0]), len(hist[0, :]), 1))
	hist_array[:, :, 0] = hist
	data.makeHistogram(hist_array, a, b, 'ff_topo')

	

	end = time.time()
	print "Time:", end-start
