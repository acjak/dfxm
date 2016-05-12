# !/bin/python
"""blah."""

from lib.getedfdata import *

import numpy as np

# import matplotlib
import matplotlib.pylab as plt
# import seaborn as sns

# import scipy.ndimage
import itertools

# import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print rank, size

if rank == 0:
	start = time.time()

test_switch = True

path = '/data/hxrm/Dislocation_november_2015/diamond/ff_strain'
bg_path = '/data/hxrm/Dislocation_november_2015/diamond/bg_ff'

filename = 'ff2_'
# filename2 = 'ff2_'
sampletitle = 'topo_fancy_txt_'
bg_filename = 'bg_ff_2x2_0p5s_'

datatype = 'strain_tt'

# poi = [1023, 1023]
# size = [600, 300]
#
# roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

roi = [200, 1800, 200, 1800]

data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)

meta = data.getMetaArray()
# hist, datx, daty = data.makeMeanGrid()
a, b, c = data.getMetaValues()

fulldata = [a,b,c]
print np.shape(fulldata)

data.setTest(True)
data.adjustOffset(False)

ab_vals = list(itertools.product(a, b))

# Chose part of data set for a specific CPU (rank).
local_n = len(ab_vals)/size
istart = rank*local_n
istop = (rank+1)*local_n
local_data = b[istart:istop]
# local_data = fulldata[istart:istop]

print len(a), len(b), len(c)

# if rank == 0:
# 	end = time.time()
# 	print "Init time: ", end-start

fig = plt.Figure()

# for i in range(len(local_data)):

for i,beta in enumerate(local_data):
	ax = fig.add_subplot(111)
	#print rank, local_data[i][0], local_data[i][1], local_data[1]
	print rank, a[len(a)/2], beta, c[len(c)/2]
	#index = data.getIndex(float(local_data[i][0]), float(local_data[i][1]), float(local_data[1]))
	index = data.getIndex(a[len(a)/2], beta, c[len(c)/2])
	print index
	img = data.getImage(index[0], False)
	ax.imshow(img, cmap='Greens', interpolation='none')
	ax.text(5, 2, beta)
	# plt.show()
	# plt.savefig('output/topofancy_diamond2/topo_ff_angletext' + str("%04d" % (i+rank*local_n))+'.png')  #
	plt.savefig(data.directory + '/topo_ff_angletext' + str("%04d" % (i+rank*local_n))+'.png')  # str(i+80+rank*local_n))
	fig.clf()
