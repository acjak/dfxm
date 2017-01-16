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

path = '/u/data/andcj/hxrm/MLL_december_2016/BaTiO3'
bg_path = '/u/data/andcj/hxrm/MLL_december_2016/BaTiO3'

filename = 'strainmesh_5s_1_'
sampletitle = filename
bg_filename = 'bg1_5s_'

datatype = 'strain_tt'

theta = 10.25

test_switch = True


poi = [600, 500]
size = [50, 50]
s = 40

roi = [
	poi[0] - size[0] / 2,
	poi[0] + size[0] / 2,
	poi[1] - size[1] / 2,
	poi[1] + size[1] / 2]

print roi

data = GetEdfData(
	path, filename, bg_path, bg_filename, roi, datatype, test_switch)
data.setTest(True)
data.adjustOffset(False)

try:
	directory = data.directory
except AttributeError:
	directory = 0
