#!python

import EdfFile

import sys
import os

# import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

import scipy.ndimage
import itertools

path = os.getcwd()

filename = sys.argv[0]
sampletitle = filename

bg_filename = 'bg1_5s_'

datatype = 'strain_tt'

poi = [1250, 1150]
# size = [100, 100]
size = [600, 300]
# poi = [1000, 1000]
# size = [1000, 1000]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

img = EdfFile.EdfFile(path + '/' + filename)



# eta = -0.01
# cols = 8

# data = GetEdfData(path, filename, bg_filename, roi, datatype)
