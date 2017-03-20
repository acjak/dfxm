#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
# import time
from lib.getedfdata import GetEdfData
from lib.dfxm import DFXM
import numpy as np
# import itertools
# from mpi4py import MPI
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
# from matplotlib.colors import colorConverter
# from scipy import ndimage
# from mpl_toolkits.axes_grid.inset_locator import inset_axes


def getimages(alpha):
	print alpha
	fig, ax = plt.subplots(1, 2, figsize=(14, 7))
	for i, oby in enumerate(alpha):
		img = data.getImage(data.getIndex(oby, -10000, -10000)[0], False)
		ff = data2.getImage(data.getIndex(oby, -10000, -10000)[0], False)
		print "Index: {}, Oby: {}, Im. mean: {}, Flatfield mean: {}".format(i, oby, np.mean(img), np.mean(ff))
		ff[ff <= 0] = 1
		cr = img * ff.mean() / ff
		cr[cr > 600] = 600
		cr[cr < 0] = 0
		# cr_filter = data.rfilter(cr, 18, 3)

		z = np.zeros((100))
		# length = 0

		proj = [90, 150, 160, 25]
		pp = np.zeros((4, 20))

		for j in range(20):
			pp[:, j] = proj[0] - j, proj[1] - j * 0.65, proj[2] - j, proj[3] - j * 0.65
			z += tools.getProjection(cr, pp[0, j], pp[1, j], pp[2, j], pp[3, j], 100)

		plotimages(cr, z, pp, oby, fig, ax)


def plotimages(im, projline, pp, oby, fig, ax):

	# xy = np.zeros(8, 2)
	# for k in range(8):
	# 	xy[k, :] = pp
	# 	x = np.linspace(pp[0, k], pp[2, k], 100)
	# 	y = np.linspace(pp[1, k], pp[3, k], 100)
	# 	ax[0].scatter(x, y, s=1)

	ax[0].plot(pp[::2, 0], pp[1::2, 0], 'b')
	ax[0].plot(pp[::2, -1], pp[1::2, -1], 'b')
	ax[0].plot([pp[::2, 0][0], pp[::2, -1][0]], [pp[1::2, 0][0], pp[1::2, -1][0]], 'b')
	ax[0].plot([pp[::2, 0][1], pp[::2, -1][1]], [pp[1::2, 0][1], pp[1::2, -1][1]], 'b')

	ax[0].imshow(im, cmap='Greys')
	ax[0].autoscale(False)
	ax[1].plot(projline)
	ax[1].set_ylim(5500, 8500)

	# fig.title('OBY pos. = {} mm'.format(oby))
	fig.savefig(directory + '/rawimage_{:06.3f}.png'.format(oby))
	ax[0].cla()
	ax[1].cla()
	# plt.close(fig)


def fastplot(alpha):

	fig, ax = plt.subplots(figsize=(10, 10))
	extent = ax.get_window_extent().transformed(
		plt.gcf().dpi_scale_trans.inverted())
	ax.set_axis_off()
	for i, oby in enumerate(alpha):
		img = data.getImage(data.getIndex(oby, -10000, -10000)[0], False)
		ff = data2.getImage(data.getIndex(oby, -10000, -10000)[0], False)
		ff[ff <= 0] = 1
		cr = img * ff.mean() / ff
		cr[cr > 600] = 600
		cr[cr < 0] = 0
		ax.imshow(cr, cmap='Greys')

		fig.savefig(directory + '/rawimage_{:06.3f}.png'.format(oby), bbox_inches=extent)
		ax.cla()


path = '/u/data/andcj/hxrm/MLL_december_2016/nanorod/part2'
bg_path = '/u/data/andcj/hxrm/MLL_december_2016/alignment'

filename = 'focusscan2_wide_10s_3_0'
filename2 = 'focusscan2_wide_10s_3_flatfield'

sampletitle = filename
bg_filename = 'bg_10s_'

datatype = 'strain_tt'
test_switch = True


poi = [602, 512]
# size = [300, 300]
size = [500, 500]

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]

print roi

motors = ['obx', 'diffrz', 'chi']

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

runvars[1] = filename2

data2 = GetEdfData(runvars)

data2.setTest(True)
data2.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0

tools = DFXM(
	path,
	data.data_files,
	directory,
	roi,
	datatype,
	data.dirhash,
	data.meta,
	test_switch)

a, b, c = data.getMetaValues()

# getimages(a)
fastplot(a)
