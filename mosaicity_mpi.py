#!/bin/python
"""blah."""

from lib.getedfdata import GetEdfData
from lib.dfxm import DFXM
# import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
# import seaborn as sns
from scipy import ndimage
import itertools
import warnings
import numpy as np

from mpi4py import MPI
import time
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

print rank, mpisize

if rank == 0:
	start = time.time()

path = '/u/data/andcj/hxrm/Al_march_2017/al/mapping/run0/layer0'
bg_path = '/u/data/andcj/hxrm/Al_march_2017/al/mapping/run0/layer0'

filename = 'mosa_frelon_far_'

sampletitle = 'mosa_r0_l0_'
bg_filename = 'mosa_frelon_far_0001_0001_0001'

datatype = 'mosaicity'

poi = [1024, 1024]
size = [600, 600]

test_switch = True

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

motors = ["diffry", "phi", "diffrz"]

runvars = [
	path,
	filename,
	bg_path,
	bg_filename,
	datatype,
	roi,
	test_switch,
	motors]

data = GetEdfData(runvars)
data.setTest(True)
data.adjustOffset(False)

alist, blist, clist = data.getMetaValues()


def mosaicity_part(istart, istop, local_n):
	fname = "tmp/imagepart_{}".format(rank)
	if os.path.isfile(fname + '.npy'):
		print "Loading data from file. " + str(rank)
		image_part = np.load(fname + '.npy')

	else:
		print "Reading data from data files. " + str(rank)
		image_part = np.zeros((local_n, size[1], 273, len(blist)), np.dtype(np.float16))
		data.roi[2] = data.roi[2] + istart
		data.roi[3] = data.roi[2] + (istop-istart)

		for i, b in enumerate(blist):
			alistx = data.meta[np.where(data.meta[:, 1] == b), 0]
			for j, a in enumerate(alistx[0]):
				ind = data.getIndex(a, b, -10000)
				image_part[:, :, j, i] = data.getImage(ind[0], False)

			if rank == 0:
				print i
			if rank == 0 and float(i) / len(blist) % 1 == 0.0:
				done = 100 * (float(i) / len(blist))
				print "Calculation is %g perc. complete..." % done
		np.save(fname, image_part)

		print "Image saved on cpu" + str(rank)
	mosaic_part = np.zeros((local_n, size[1], 2), np.dtype(np.float16))

	if rank == 0:
		print "\nStarting mosaicity fitting. \n"

	image_part[image_part < 2] = 0

	for k in range(local_n):
		for l in range(size[1]):
			if np.count_nonzero(image_part[k, l, :, :]) > 1000:
				try:
					com = ndimage.measurements.center_of_mass(image_part[k, l, :, :])
					# print com, rank
				except KeyError:
					print "KeyError", str(rank)
					com = [6., 6.]
				mosaic_part[k, l, :] = com
				# print mosaic_part[k, l, :]

		if rank == 0 and 10 * float(k) / local_n % 1 == 0.0:
			done = 100 * (float(k) / local_n)
			print "Calculation is %g perc. complete..." % done

	print np.mean(mosaic_part), np.max(mosaic_part), rank
	np.save(data.directory + '/mosapart_{}'.format(rank), mosaic_part)
	mosaic_part[0, 0, 0] = rank
	return mosaic_part


def mosaicityMPI():
	print rank, mpisize

	ypix = (data.roi[1] - data.roi[0])
	# Chose part of data set for a specific CPU (rank).
	local_n = ypix / mpisize
	istart = rank * local_n
	istop = (rank + 1) * local_n
	# local_data = blist[istart:istop]

	# Calculate strain on part of data set.
	mosaic_part = mosaicity_part(istart, istop, local_n)

	# CPU 0 (rank 0) combines data parts from other CPUs.
	if rank == 0:
		# Make empty arrays to fill in data from other cores.
		recv_buffer = np.zeros((local_n, size[1], 2), np.dtype(np.float16))
		mosaic_array = np.zeros((size[0], size[1], 2), np.dtype(np.float16))
		recv_buffer = [0., 0.]

		datarank = mosaic_part[0, 0, 0]
		mosaic_part[0, 0, 0] = 0
		mosaic_array[istart:istop, :, :] = mosaic_part
		for i in range(1, mpisize):
			try:
				comm.Recv(recv_buffer, MPI.ANY_SOURCE)
				datarank = int(recv_buffer[0, 0, 0])
				recv_buffer[0, 0, 0] = 0
				mosaic_array[datarank * local_n:(datarank + 1) * local_n, :, :] = recv_buffer
			except Exception:
				print "MPI error."

	else:
		# all other process send their results
		a = [54389.54, 5435.9545]
		# comm.Send(a, dest=0)
		# comm.Send(mosaic_part, dest=0)

	# root process prints results
	if comm.rank == 0:
		return mosaic_array


mosaic_array = mosaicityMPI()
if rank == 0:
	print np.shape(mosaic_array)
	np.save(data.directory + '/mosaic_array.npy', mosaic_array)
