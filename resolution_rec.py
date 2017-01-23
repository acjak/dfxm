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

path = sys.argv[1]  # '/u/data/andcj/Resolution_march_2016/rollscan_center'
bg_path = '/u/data/andcj/Resolution_march_2016/rollscan_top_far'

filename = 'scan1_'
bg_filename = 'scan1_0960.edf'
# filename2 = 'ff2_'
sampletitle = sys.argv[2]  # 'rollscan_center_roll'
datatype = 'res_paper'

test_switch = True

poi = [1024, 1024]
# size = [2, 20]
size = [2048, 2048]

pixelsize = 180.  # 180.

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]

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
	gs1.update(wspace=0.025, hspace=0.03)

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
		index = data.getIndex(i, b[len(b) / 2], c[0])
		index_list.append(index[0])

	xr = a - data.alpha0

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')
		return tools.makeNewGaussArrayMPI(data.imgarray, 1, xr), index_list


def makeIndexList_TT(a, b, c):
	index_list = []

	for i in b:
		# print rank, float(a[i]-data.alpha0), float(b[i]-data.beta0)#, i
		index = data.getIndex(a[len(a) / 2], i, c[0])
		index_list.append(index[0])

	xr = b - data.beta0

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		data.makeImgArray(index_list, 50, 'linetrace')
		return tools.makeNewGaussArrayMPI(data.imgarray, 1, xr), index_list


def saveGaussArray(sampletitle, gaussarray):
	np.save(data.directory + '/' + sampletitle + '_data.npy', gaussarray)

if sys.argv[3] == 'TT':
	gaussarray, index_list = makeIndexList_TT(a, b, c)
if sys.argv[3] == 'ROLL':
	gaussarray, index_list = makeIndexList_ROLL(a, b, c)

if rank == 0:
	print np.shape(data.imgarray)
	if test_switch:
		plotImageArray(sampletitle)
		# tools.plotPPOS(gaussarray)
		tools.plotPPOSBig(gaussarray)
		saveGaussArray(sampletitle, gaussarray)
		# tools.plotStrain(strainpic, data.imgarray)

	end = time.time()
	print "Time:", end - start, "seconds."
