#!/bin/python
"""blah."""

from lib.getedfdata import *
from lib.dfxm import *
from lib.gauss import *
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
# import seaborn as sns

# import scipy.ndimage
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

# ESRF path:
# /data/id06/inhouse/2015/run5/

path = '/data/hxrm/Dislocation_november_2015/diamond/ff_strain'
bg_path = '/data/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff2_'
# filename2 = 'ff2_'
sampletitle = 'topo_strain'
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'strain_tt'

poi = [500, 500]
size = [500, 500]

diffrz_pos1 = 13
diffrz_pos2 = 17

test_switch = True

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)


a, b, c = data.getMetaValues()
ab_vals = list(itertools.product(a, b))


def plotImageArray(imgarray1, imgarray2):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(7, 7)
	gs1.update(wspace=0.025,  hspace=0.3)

	for i in range(len(imgarray1[:, 0, 0])):
		ta = np.ones((len(imgarray1[0, :, 0]), len(imgarray1[0, 0, :]), 4),  dtype=np.uint8)*0
		ta[:, :, 3] = 255
		ta[:, :, 0] = 255*imgarray1[i, :, :]/np.max(imgarray1[i, :, :])
		ta[:, :, 2] = 255*imgarray1[i, :, :]/np.max(imgarray1[i, :, :])
		ta[:, :, 1] = 255*imgarray2[i, :, :]/np.max(imgarray2[i, :, :])

		# print 255*imgarray1[i, :, :]/np.max(imgarray1[i, :, :])
		axarr = plt.subplot(gs1[i])
		axarr.imshow(ta)
		axarr.set_title('%.4f' % (float(c[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	plt.savefig(data.directory + '/%s_array.pdf' % ('overview'))


def makeImgList(a, b, c, diffrz_pos):
	index_list = []

	for i in range(len(c)):
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(float(a[diffrz_pos]), float(b[diffrz_pos]), c[i])
		index_list.append(index[0])

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')

		return data.imgarray


imgarray1 = makeImgList(a, b, c, diffrz_pos1)
imgarray2 = makeImgList(a, b, c, diffrz_pos2)

if rank == 0:
	plotImageArray(imgarray1, imgarray2)

	end = time.time()
	print "Time:", end-start, "seconds."
