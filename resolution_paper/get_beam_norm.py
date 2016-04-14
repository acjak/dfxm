# !/bin/python
"""blah."""
import EdfFile

import numpy as np
import scipy.interpolate as inter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
#  import seaborn as sns

# import scipy.ndimage
import itertools

path = '/data/hxrm/Resolution_march_2016/nf_images'

filename = 'nf_s5_0p5x0p5_0001.edf'
# filename2 = 'ff2_'
sampletitle = 'beam_normalization'

poi = [1160, 1050]
size = [800, 800]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

file_with_path = path + '/' + filename
img = EdfFile.EdfFile(file_with_path)
im = img.GetData(0).astype(np.int64)[roi[2]:roi[3], roi[0]:roi[1]]

xr = np.arange(-1024., 1024., 2048./800)
yr = np.arange(-1024., 1024., 2048./800)

im_interp = inter.interp2d(xr,yr,im, kind='cubic')

poi_interp = []

img_zoom = im_interp(np.arange(-100.,100,1.), np.arange(-100.,100,1.))

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

ax2.imshow(im, cmap="Greens")
ax.imshow(img_zoom, cmap="Greens")

ax.set_xlabel('[um]')
ax.set_ylabel('[um]')
fig.savefig('output/beam_norm.png')
fig2.savefig('output/beam_norm_full.png')
