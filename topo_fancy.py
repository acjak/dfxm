# !/bin/python
"""blah."""

from lib.getedfdata import *
from lib.gauss import *
import numpy as np

# import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

import scipy.ndimage
import itertools

import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print rank, size

if rank == 0:
	start = time.time()


path = '/data/hxrm/Dislocations_november_2015/diamond/nearfield'

filename = 'topotomo_nf1_'
sampletitle = filename
bg_filename = 'bg_1x1_1s_'

datatype = 'topotomo'

# poi = [1023, 1023]
# size = [600, 300]
#
# roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

roi = [200, 1800, 200, 1800]

data = GetEdfData(path, filename, bg_filename, roi, datatype)

meta = data.getMetaArray()
hist, datx, daty = data.makeMeanGrid()
a, b = data.getMetaValues()

data.setTest(True)
data.adjustOffset(False)

ab_vals = list(itertools.product(a, b))

print a
print b

# Chose part of data set for a specific CPU (rank).
local_n = len(ab_vals)/size
istart = rank*local_n
istop = (rank+1)*local_n
local_data = ab_vals[istart:istop]

# if rank == 0:
# 	end = time.time()
# 	print "Init time: ", end-start

for i in range(len(local_data)):
    index = data.getIndex(float(local_data[i][0]), float(local_data[i][1]))
    print index, rank, local_data[i][0], local_data[i][1]
    img = data.getImage(index[0], False)
    plt.imshow(img, cmap='Greens', interpolation='none')
    # plt.show()
    plt.savefig('output/topofancy_diamond/topo_nf1_' + str("%04d" % (i+rank*local_n))+'.png') # str(i+80+rank*local_n))
