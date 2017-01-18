#!/bin/python
# -*- coding: utf-8 -*-
"""blah."""
import time
from lib.getedfdata import *
from lib.dfxm import *
import numpy as np
# import itertools
from mpi4py import MPI
import matplotlib
matplotlib.rc('font', family='DejaVu Sans')
import matplotlib.pylab as plt
# from matplotlib.colors import colorConverter
# from scipy import ndimage
# from mpl_toolkits.axes_grid.inset_locator import inset_axes


def saveimages(alpha):
	print alpha
	for i, oby in enumerate(alpha):
		img = data.getImage(data.getIndex(oby, -10000, -10000)[0], False)
		ff = data2.getImage(data.getIndex(oby, -10000, -10000)[0], False)

		cr = img * ff.mean() / ff
		cr[cr > 600] = 600
		cr[cr < 0] = 0

		cr_filter = data.rfilter(cr, 18, 3)

		plt.imshow(cr, cmap='Greys')
		plt.savefig(directory + '/rawimage_{:06.3f}.png'.format(oby))


path = '/u/data/andcj/hxrm/MLL_december_2016/nanorod/part2'
bg_path = '/u/data/andcj/hxrm/MLL_december_2016/alignment'

filename = 'focusscan2_wide_10s_3_0'
filename2 = 'focusscan2_wide_10s_3_flatfield'

sampletitle = filename
bg_filename = 'bg_10s_'

datatype = 'strain_tt'
test_switch = True


poi = [602, 512]
size = [300, 300]

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

a, b, c = data.getMetaValues()

saveimages(a)
