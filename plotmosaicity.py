#!/bin/python
"""blah."""

import math
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.colors import hsv_to_rgb
from math import atan2, degrees, pi

mosaic_unfinished = np.load('tmp/mosaic_array.npy')
mosaic_data = np.zeros_like(mosaic_unfinished)

fdir = '/Users/andcj/bitsync/Documents/PostDoc/scripts/output/210317-0334/'

for i in range(16):
	fname = fdir + 'mosapart_{}.npy'.format(i)
	mosaic_part = np.load(fname)

	mosaic_data[i*37:(i+1)*37, :, :] = mosaic_part


def getangle(x2, y2):
	rads = atan2(-y2, x2)
	rads %= 2*pi
	return rads


le = 100
leh = le/2

hsva = np.zeros((le, le, 3))

for x in range(-leh, leh):
	for y in range(-leh, leh):

		dist = math.sqrt((x)**2 + (y)**2)
		ang = getangle(x, y)

		hsva[x-leh, y-leh, 0] = abs(ang) / (2 * np.pi)
		hsva[x-leh, y-leh, 1] = abs(dist / math.sqrt((leh)**2 + (leh)**2))
		hsva[x-leh, y-leh, 2] = 1


HSV = np.dstack((hsva[:, :, 0], hsva[:, :, 1], hsva[:, :, 2]))
RGB = hsv_to_rgb(HSV)

outputpic = np.zeros((np.shape(mosaic_data)[0], np.shape(mosaic_data)[1], 3))
outputpic[:, :, 2] = 1.

for i in range(600):
	for j in range(600):
		theta_phi = mosaic_data[i, j, :]
		if theta_phi[0] != [0.] and theta_phi[1] != [0.]:
			try:
				u = int(100 * (theta_phi[0]/273 - 0.5))
				v = int(100 * (theta_phi[1]/42 - 0.5))
			except ValueError:
				print "ValueError", u, v
			try:
				outputpic[i, j, :] = hsva[u+50, v+50, :]
			except IndexError:
				print "IndexError", u, v
		# else:
		# 	print "theta_phi is zero."

print outputpic[450:460, 290:291, :]
outputpic[:, :, 1] *= 20
print "Making plot"
fig, ax = plt.subplots(1, 2)
output_rgb = hsv_to_rgb(outputpic)
ax[0].imshow(RGB, extent=[0, 100, 0, 100], aspect=1)
ax[1].imshow(output_rgb)
# ax[1].imshow(mosaic_data[:, :, 0], origin="lower")
plt.show()
