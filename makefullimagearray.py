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

path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/ff_topo_2'
bg_path = '/u/data/andcj/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff3_'
# filename2 = 'ff2_'
sampletitle = 'fulltopo'
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'topotomo'

# poi = [500, 500]
# size = [1000, 1000]

poi = [512, 512]
size = [200, 200]


test_switch = True

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0

# tools = DFXM(path, data.data_files, directory, roi, datatype, data.dirhash, data.meta, test_switch)

# tools = DFXM(directory, roi, datatype, test_switch)

a, b, c = data.getMetaValues()
# ab_vals = list(itertools.product(a, b))
def line(x, a, b):
	return a*x+b

def fitLine(x, y):
	from scipy.optimize import curve_fit

	try:
		popt, pcov = curve_fit(line, x, y, p0=[0, 30], maxfev=10000)
		return popt, pcov
	except RuntimeError:
		pass

def allFiles(a, b):
	index_list = range(len(data.meta))

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		imgarray = data.makeImgArray(index_list, 50, 'linetrace')

	if rank == 0:
		np.save('/u/data/andcj/tmp/alpha.npy', a)
		np.save('/u/data/andcj/tmp/beta.npy', b)
		print np.shape(imgarray)
		reshapedarray = np.reshape(imgarray,[41, 181, np.shape(imgarray)[1], np.shape(imgarray)[2]])

		print np.shape(reshapedarray)
		# reshapedarray[18:22,:,:,:] = 0
		np.save('/u/data/andcj/tmp/largetest_200x200.npy', reshapedarray)

if test_switch:
	allFiles(a, b)
	# strainpic, index_list = makeIndexList(a, b, c, diffrx_pos)
else:
	makeFullTomoArray(a, b, c, 12)

#
# if rank == 0:
# 	if test_switch:
# 		plotImageArray(diffrx_pos)
# 		tools.plotStrain(strainpic, data.imgarray, diffrx_pos)
#
# 	end = time.time()
# 	print "Time:", end-start, "seconds."
