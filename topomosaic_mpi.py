# !/bin/python
"""blah."""

from lib.getedfdata import *
from lib.dfxm import *
from lib.gauss import *
import numpy as np

#  import matplotlib
import matplotlib.pylab as plt
#  import seaborn as sns

# import scipy.ndimage
import itertools

# import warnings

from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
mpisize = comm.Get_size()

print rank, mpisize

if rank == 0:
	start = time.time()

# path = '/data/id06/inhouse/2015/run5/diamond/ff_mosaicity_tomo_2'
# bg_path = '/data/id06/inhouse/2015/run5/diamond/bg_ff'
#
# filename = 'ff1_'
# # filename2 = 'ff2_'
# sampletitle = 'mosaicity_1'  # filename
# bg_filename = 'bg_ff_2x2_0p5s_'
#
# datatype = 'mosaicity'

path = '/data/hxrm/BaTiO3_january_2016/pole1'
bg_path = '/data/hxrm/BaTiO3_january_2016/bg'

filename = 'pole1_zap_v2_'
# filename2 = 'ff2_'
sampletitle = 'BaTiO3_fancy'
bg_filename = 'bg_frelon_far_zap_3000ms_frelon_far_0001_'

path2 = '/data/hxrm/BaTiO3_january_2016/pole2'
filename2 = 'pole2_zap_v1_'

datatype = 'topotomo'

test_switch = True

poi = [512, 500]
poi2 = [500, 500]
size = [500, 500]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]
roi2 = [poi2[0]-size[0]/2, poi2[0]+size[0]/2, poi2[1]-size[1]/2, poi2[1]+size[1]/2]


data = GetEdfData(path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

data2 = GetEdfData(path2, filename2, bg_path, bg_filename, roi2, datatype, test_switch)
data2.setTest(True)
data2.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0
tools = DFXM(path, data.data_files, directory, roi, datatype, data.dirhash, data.meta, test_switch)

diffrz_pos1 = 3
diffrz_pos2 = 5


a, b, c = data.getMetaValues()
# c, d = data2.getMetaValues()
ab_vals = list(itertools.product(a, b))
# cd_vals = list(itertools.product(c, d))

local_n = len(a)/mpisize
istart = rank*local_n
istop = (rank+1)*local_n
local_a = a[istart:istop]
# local_c = c[istart:istop]

# print local_a

fig = plt.figure(frameon=False)
fig.set_size_inches(10, 10)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)


for i in range(len(local_a)):

	index1 = data.getIndex(float(local_a[i]), -10000, -10000)
	index2 = data2.getIndex(float(local_a[i]), -10000, -10000)

	# index1 = data.getIndex(float(local_a[i]), float(b[diffrz_pos1]), -10000)
	# index2 = data.getIndex(float(local_a[i]), float(b[diffrz_pos2]), -10000)

	img1 = data.getImage(index1[0], False)
	img2 = data2.getImage(index2[0], False)

	if np.mean(img1) < 10.:
		img1 = np.zeros((size[0],size[1]))

	if np.mean(img2) < 10.:
		img2 = np.zeros((size[0],size[1]))
	# print index1[0], local_a[i], i+rank*local_n, img1.mean()
	# img2 = data.getImage(index2[0], False)

	img1[img1 < 0] = 0
	img2[img2 < 0] = 0

	ta = np.ones((len(img1[:, 0]), len(img1[0, :]), 4),  dtype=np.uint8)*0

	ta[:, :, 3] = 255
	ta[:, :, 0] = 255*img1/np.max(img1)
	ta[:, :, 2] = 255*img2/np.max(img2)
	ta[:, :, 1] = 255*img2/np.max(img2)

	# if np.mean(img1) < 0.01 and np.mean(img2) < 0.01:
	# 	ta[:, :, 3] = 255
	# 	ta[:, :, 0] = 0
	# 	ta[:, :, 2] = 0
	# 	ta[:, :, 1] = 0

	ax.imshow(ta,interpolation="none", origin="lower")
	# ax.imshow(img1,interpolation="none", origin="lower", cmap='Greens')
	# plt.colorbar()
	fig.savefig(data.directory + '/topo_im_' + str('%04d' % (i+rank*local_n)) + '.png')
	# fig.clf()

if rank == 0:
	meta = data.getMetaArray()
	data.getMeanData()
	hist, datx, daty = tools.makeMeanGrid(data.data_mean)
	hist_array = np.zeros((len(hist[:, 0]), len(hist[0, :]), 1))
	hist_array[:, :, 0] = hist
	tools.makeHistogram(hist_array, a, b, 'ff_topo')

	end = time.time()
	print "Time:", end-start
