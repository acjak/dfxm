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


def makeGvector(bv, xylin, ang):
	gv = np.ones((3, len(xylin), len(xylin)))  # G-vector

	xcpos = len(gv[0, :, 0])/2
	ycpos = len(gv[0, 0, :])/2

	prefix = 2*math.pi/bv
	# prefix = 2*math.sin(ang)

	for i in range(len(xylin)):
		for j in range(len(xylin)):
			if i != xcpos and j != ycpos:
				gv[0, i, j] = prefix*bv/(2*math.pi)*(xylin[j])/((xylin[i])**2 + (xylin[j])**2)
				gv[1, i, j] = prefix*bv/(2*math.pi)*(xylin[i])/((xylin[i])**2 + (xylin[j])**2)
				gv[2, i, j] = prefix

	return gv

def makeQvector(wl, omega, xi):
	q = [0, \
		(2*math.pi/wl)*2*np.sin(xi/2)*np.sin(omega-xi/2), \
		(2*math.pi/wl)*2*np.sin(xi/2)*np.cos(omega-xi/2)]

	return q


def fullArray(gv, wl, xylin, omega0):
	# linearray = np.zeros((len(data.alphavals),len(data.betavals)))

	a, b = data.getMetaValues()

	dist = np.zeros((3))

	sigmax = 1E-5
	sigmay = 5E-4
	sigmaz = 1E-4

	om_start = 18
	xi_start = 10

	with open(data.infFile,  "a") as myfile:
		myfile.write("sigmax: %s \n" % (str(sigmax)))
		myfile.write("sigmay: %s \n" % (str(sigmay)))
		myfile.write("sigmaz: %s \n" % (str(sigmaz)))
		myfile.write("\n")
		myfile.write("om_start: %s \n" % (str(om_start)))
		myfile.write("xi_start: %s \n" % (str(xi_start)))
		myfile.write("\n\n\n")

	k = 0

	fig = plt.figure(figsize=(18, 4))
	gs1 = matplotlib.gridspec.GridSpec(1,10)
	gs1.update(wspace=0.02,  hspace=0.2)

	sns.set_style("white")

	for i in range(1):  # len(a)/2):
		for j in range(10):
			# Empty array for model calculations.
			matrix_part = np.zeros((len(xylin), len(xylin)))

			# Omega and xi defined from the peak position in degrees.
			omega = b[j+om_start]-10.995
			xi = a[i+xi_start]-22.086

			# Import image from data set.
			ind = data.getIndex(float(a[i+xi_start]),float(b[j+om_start]))[0]
			img = data.getImage(ind,False)

			# Creat subplot in plot array.
			axarr = plt.subplot(gs1[k])
			# axarr.set_title('%.3f %.3f' % (omega, xi),fontsize=8)
			axarr.set_title('%.4f %.4f' % (b[j+om_start], a[i+xi_start]),fontsize=8)

			# Get an average across the dislocation line.
			line = np.sum(img, axis=1)/len(img[0, :])

			# Make Q vector from the same angle as data.
			q = makeQvector(wl, np.radians(omega0+omega), np.radians(2*omega0+xi))

			# Length of q.
			ql = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)
			print q[1]**2, q[2]**2
			print ql

			sigma_a = 7.698e-5
			sigma_a = 2.8e-4

			sigmax = 2*ql*sigma_a
			sigmay = ql*1e-5
			sigmaz = ql*(1/math.tan(np.radians(omega0+omega)))*sigma_a


			# Loop over all pixels in the cross-sectional model image.
			for c in range(len(xylin)):
				for d in range(len(xylin)):

					# Find the distance from Q to G vector in x, y and z.
					for m in range(len(q)):
						dist[m] = abs(gv[m, c, d] - q[m])

					# Use the distance to fold with a Gaussian.
					matrix_part[c, d] += math.exp(-(dist[0]**2)/(2*sigmax**2))
					matrix_part[c, d] += math.exp(-(dist[1]**2)/(2*sigmay**2))
					matrix_part[c, d] *= math.exp(-(dist[2]**2)/(2*sigmaz**2))

			# Plot data.
			axarr.plot(xylin, line)
			axarr.plot([xylin[np.argmax(line)], xylin[np.argmax(line)]], [0, np.max(line)])

			# Make a line from the model and plot it on the same subplot.
			model_line = np.sum(matrix_part[:,:], axis=1)
			axarr.plot(xylin, np.max(line)*model_line/np.max(model_line))
			# axarr.plot(xylin, 10*model_line)

			# axarr.imshow(matrix_part, cmap='Greens')
			# axarr.set_ylim(0,2000)
			axarr.xaxis.set_major_formatter(plt.NullFormatter())
			axarr.yaxis.set_major_formatter(plt.NullFormatter())

			k += 1

	fig.savefig(data.directory + '/comp_screwdisl.pdf')
	plt.show()


path = '/Users/andcj/hxrm_data/disl_may_2015/dislocations/strain'

filename = 'strainmap_tt_2'
sampletitle = filename
bg_filename = 'bg1_5s_'

datatype = 'strain_tt'

#poi = [750, 750]
# size = [800, 200]

poi = [1220, 1150]
# size = [100, 100]
size = [30, 80]
# poi = [1000, 1000]
# size = [1000, 1000]

roi = [poi[0]-size[0]/2, poi[0]+size[0]/2, poi[1]-size[1]/2, poi[1]+size[1]/2]

eta = -0.01
cols = 8

data = GetEdfData(path, filename, bg_filename, roi, datatype)
# data.printMeta()
meta = data.getMetaArray()
hist, datx, daty = data.makeMeanGrid()
a, b = data.getMetaValues()

data.setTest(True)
data.adjustOffset(True)
# alignToCenter(meta, a, b, hist)
# overlapImg(meta, a, b, hist, [10, 26], [11, 26])
# compareTheta(meta, eta, cols)

omega0 = 10.94822
ang = np.radians(omega0)

bv = 0.192  # Burger's vector
wl = 0.07293  # Wavelength

ln = size[1]*80.

xylin = np.ogrid[-ln/2:ln/2:80]

gv = makeGvector(bv, xylin, ang)
print xylin
fullArray(gv, wl, xylin, omega0)

# bins, length, xr = diagonalPlot(meta, a, b, hist, 1, 0.53, datatype)
