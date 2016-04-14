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

path = '/data/hxrm/Resolution_march_2016/vignet_mesh'
bg_path = '/data/hxrm/Resolution_march_2016/rollscan_top_far'

filename = 'mesh1_'
bg_filename = 'scan1_0960.edf'
# filename2 = 'ff2_'
sampletitle = 'ffdet_normalization'
datatype = 'ffdet_normalization'

test_switch = True

poi = [1024, 1024]
# size = [2, 20]
size = [64, 64]

pixelsize = 180.  # 180.

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

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

print a, b

cal_image = np.zeros((2048, 2048))
a0 = a[len(a)/2]
b0 = b[len(b)/2]
lc = len(cal_image)
xlen = size[0]/2
ylen = size[1]/2

maxcurrent = max(data.meta[:, 3])

# pixel per mm on detector:
ppm = 630

fig, ax = plt.subplots()

ax.imshow(cal_image, cmap='Greens', interpolation=None)
ax.autoscale(False)

n = 0
for i in a:
    for j in b:
        xpos = (a0-i)*ppm
        ypos = -(b0-j)*ppm
        print ypos, j, xpos, i
        print maxcurrent, data.meta[n, 3]
        ind = data.getIndex(i, j, c[0])
        data.roi = [lc/2-xpos-xlen, lc/2-xpos+xlen, lc/2-ypos-ylen, lc/2-ypos+ylen]

        # fig2, ax2 = plt.subplots()
        # ax2.imshow(img, cmap='Greens', interpolation=None)
        # fig2.savefig(data.directory + '/full_img_%g.png' % n)
        # fig2.clf()

        img = data.getImage(ind[0], False)
        img = maxcurrent*img/data.meta[n, 3]
        # Add ring current normalization and check the detx/detz directions of the result.
        n += 1
        print "Image #:", n

        roi = data.roi

        # ax.plot([roi[0], roi[0]], [roi[2], roi[3]], color='black')
        # ax.plot([roi[1], roi[1]], [roi[2], roi[3]], color='black')
        # ax.plot([roi[0], roi[1]], [roi[3], roi[3]], color='black')
        # ax.plot([roi[0], roi[1]], [roi[2], roi[2]], color='black')

        cal_image[lc/2-ypos-ylen:lc/2-ypos+ylen, lc/2-xpos-xlen:lc/2-xpos+xlen] = img

ax.imshow(cal_image, cmap='Greens', interpolation=None)
fig.savefig(data.directory + '/cal_image.png')
