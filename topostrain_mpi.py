#!/bin/python
"""blah."""

from lib.getedfdata import *
from lib.dfxm import *
from lib.gauss import *
# import numpy as np

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

poi = [400, 675]
size = [10, 10]

diffrx_pos = 31

test_switch = False

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0

tools = DFXM(directory, roi, datatype, test_switch)

a, b, c = data.getMetaValues()
ab_vals = list(itertools.product(a, b))


def plotImageArray(diffrx_pos):
	plt.figure(figsize=(14, 14))
	gs1 = matplotlib.gridspec.GridSpec(6, 6)
	gs1.update(wspace=0.025,  hspace=0.03)

	for i in range(len(data.imgarray[:, 0, 0])):
		img = data.imgarray[i, :, :]
		axarr = plt.subplot(gs1[i])

		axarr.imshow(img, cmap="Greens")
		axarr.set_title('%.4f, %.4f' % (float(a[i]), float(b[i])))
		axarr.xaxis.set_major_formatter(plt.NullFormatter())
		axarr.yaxis.set_major_formatter(plt.NullFormatter())

	plt.savefig(data.directory + '/%s_%s_array.pdf' % ('topostrain', str(c[diffrx_pos])))


def makeIndexList(a, b, c, diffrx_pos):
	index_list = []

	for i in range(len(a)):
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(float(a[i]), float(b[i]), c[diffrx_pos])
		index_list.append(index[0])

	xr = b-data.beta0

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')
		return tools.makeStrainArrayMPI(data.imgarray, 1, xr, data.beta0), index_list


strainpic, index_list = makeIndexList(a, b, c, diffrx_pos)

if rank == 0:
	if test_switch:
		plotImageArray(diffrx_pos)
		tools.plotStrain(strainpic, data.imgarray)

	end = time.time()
	print "Time:", end-start, "seconds."
